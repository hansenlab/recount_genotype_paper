library(optparse)

option_list <- list(
  make_option(c("-s", "--metadata"), type = "character",
              help = "Input: Study metadata in .csv format."),
  make_option(c("-r", "--ref"), type = "character",
              help = "Input: Reference matrix as a data.table in .rds format"),
  make_option(c("-a", "--alt"), type = "character",
              help = "Input: Alternative matrix as a data.table in .rds format"),
  make_option(c("-g", "--trueGenotype"), type = "character",
              help = "Input: Ground truth genotype matrix as a data.table in .rds format"),
  make_option(c("-p", "--prior"), type = "character",
              help = "Output: Prior based on empirical distribution of training data, in .rds format."),
  make_option(c("-m", "--model"), type = "character",
              help = "Output: Trained model, in .rds format"),
  make_option(c("-l", "--modelLattice"), type = "character",
              help = "Output: Trained model lattice, in .rds format"),
  make_option(c("-c", "--coverageThreshold"), type = "numeric", default = 3,
              help = "Option: Coverage threshold to filter. Default: 3"),
  make_option(c("-t", "--trainSampling"), type = "numeric", default = 1,
              help = "Option: Sampling proportion for training (per genotype). Default: 1 (no sampling)."))

opt <- parse_args(OptionParser(option_list = option_list))

if (length(opt$metadata) == 0 | length(opt$ref) == 0 |
    length(opt$alt) == 0 | length(opt$trueGenotype) == 0 |
    length(opt$prior) == 0 | length(opt$model) == 0 |
    length(opt$modelLattice) == 0) {
  stop("Not all arguments provided. Check --help for description.")
}

library(tidyverse)
library(stats4)
library(broom)
library(caret)
library(data.table)
library(mgcv)

report_time_and_mem <- function(msg) {
  cat(paste0(msg, "\n"))
  print(Sys.time())
  print(gc())
}

filter_by_index <- function(input_matrix, idx) {
  input_matrix$ref<- input_matrix$ref[idx, ]
  input_matrix$alt <- input_matrix$alt[idx, ]
  input_matrix$geno <- input_matrix$geno[idx, ]
  return(input_matrix)
}

filter_by_mean_coverage <- function(input_matrix, coverageThreshold) {
  #Keep SNPs that have an average coverage (across samples) >= coverageThreshold. 
  total_matrix <- input_matrix$ref + input_matrix$alt
  avgCov <- apply(total_matrix, 1, mean)
  rm(total_matrix)
  keep_idx <- which(avgCov >= coverageThreshold)
  return(filter_by_index(input_matrix, keep_idx))
}

recode_true_genotype <- function(true_genotype) {
  #Transform true_genotype into numerical values.
  true_genotype <- as.matrix(true_genotype)
  true_genotype[true_genotype == "0/0"] <- "0" #Homozygous reference
  true_genotype[true_genotype == "0/1" | true_genotype == "1/0" ] <- "1" #Heterozygous
  true_genotype[true_genotype == "1/1"] <- "2" #Homozygous alternative
  true_genotype[true_genotype == "./." | true_genotype == "."] <- NA # Missing
  return(true_genotype)
}

raw_to_training_data <- function(metadata, ref_matrix, alt_matrix, true_genotype, coverageThreshold, trainSampling) {

  input_matrix <- list(ref = ref_matrix, 
                     alt = alt_matrix,
                     geno = recode_true_genotype(true_genotype))

  #Filter by mean coverage specified in `coverageThreshold`. 
  input_matrix <- filter_by_mean_coverage(input_matrix, coverageThreshold)
  
  #Downsample training data by as specified in `trainSampling`.
  downsample_index <- sample(1:nrow(input_matrix$geno), as.numeric(trainSampling) * nrow(input_matrix$geno))
  input_matrix <- filter_by_index(input_matrix, downsample_index)
  
  #Transform coverage ref and alt counts to M, S matrix.
  M <- as.matrix(log2((input_matrix$ref + 1) / (input_matrix$alt + 1)))
  S <- as.matrix(log2(sqrt((input_matrix$ref + 1) * (input_matrix$alt + 1))))
  true_genotype <- input_matrix$geno
  rm(input_matrix)

  #Calculate mean(M), mean(S), sd(SD) values, where for each SNP we average across samples from the same tissues.
  training_data <- vector("list", nrow(true_genotype)) #create a list with preset size.
  cat("Total number of SNPs to iterate: ", nrow(true_genotype), "\n")
  for(snp_i in 1:nrow(true_genotype)) {
    if(snp_i %% 10000 == 0) {
      cat(snp_i, "\n")
    }
    dt = data.table(M = M[snp_i ,],
                    S = S[snp_i ,],
                    genotype = true_genotype[snp_i ,],
                    tissue = metadata$tissue)
    dt = dt[!is.na(dt$genotype) ,]
    dt = dt[, .(mu_M = mean(M), sd_M = sd(M), mu_S = mean(S), n = .N),
              by = .(genotype, tissue)] 
    training_data[[snp_i]] <- dt
  }
  training_data <- do.call(rbind, training_data)
  training_data <- training_data[!is.na(training_data$sd_M) ,]

  #Estimate the true genotype distribution for the prior.
  true_genotype_distribution <- sapply(c("0", "1", "2"), function(g){
    length(which(true_genotype == g))/length(which(!is.na(true_genotype)))
  })

  return(list(training = training_data, 
              prior = true_genotype_distribution))
}

train_model <- function(training_data) {
  model <- list()
  pdf("test.pdf", width = 8, height = 8)
  #Iterate through each genotype 0, 1, 2. 
  for(my_genotype in c("0", "1", "2")) {
    cat("training genotype: ", my_genotype, "\n")
    training_data_g <- training_data %>% filter(genotype == my_genotype)
    #Mu_M linear model
    model_mu_M <- train(mu_M ~ mu_S, 
                        data = training_data_g, 
                        method = "lm",
                        trControl = trainControl(method = "none",
                                                 returnData = F,
                                                 returnResamp = "none"))  
    #sd(M) GAM model
    model_sd_M <- train(sd_M ~ mu_S, 
                        data = training_data_g, 
                        method = "gam", 
                        trControl = trainControl(method = "none",
                                                 returnData = F,
                                                 returnResamp = "none"))
    model[[my_genotype]] = list(model_mu_M, model_sd_M)
    print(ggplot(training_data_g[sample(1:nrow(training_data_g), 5000) ,], aes(x = mu_S, y = mu_M)) + geom_point())
  }
  dev.off()
  return(model)
}

generate_lattice <- function(model) {
  #IMPORTANT CONSTANTS
  max_S <- 8.5 #max value for S in generating our model lattice.
  max_ref_times_alt <- (2^max_S)^2  #max value for ref*alt in generating our model lattice.
  
  #Create lattice by transforming integer values to log2(sqrt()) space. 
  lattice_vals <- log2(sqrt(seq(1, max_ref_times_alt)))
  lattice <- data.frame(value = lattice_vals)
  #make predictions for lattice. 
  lattice$mu_M_prediction_1 <- predict(model[[1]][[1]], data.frame(mu_S = lattice$value))
  lattice$sd_M_prediction_1 <- predict(model[[1]][[2]], data.frame(mu_S = lattice$value))
  lattice$mu_M_prediction_2 <- predict(model[[2]][[1]], data.frame(mu_S = lattice$value))
  lattice$sd_M_prediction_2 <- predict(model[[2]][[2]], data.frame(mu_S = lattice$value))
  lattice$mu_M_prediction_3 <- predict(model[[3]][[1]], data.frame(mu_S = lattice$value))
  lattice$sd_M_prediction_3 <- predict(model[[3]][[2]], data.frame(mu_S = lattice$value))

  lattice <- as.data.table(lattice)
  setkey(lattice, value)

  return(lattice)
}


#######################
#Main code starts here.

report_time_and_mem("Loading and transforming data for training.")
raw_to_training_data_result  <- raw_to_training_data(metadata = read.csv(opt$metadata),
                                                     ref_matrix = readRDS(opt$ref),
                                                     alt_matrix = readRDS(opt$alt),
                                                     true_genotype = readRDS(opt$trueGenotype),
                                                     coverageThreshold = as.numeric(opt$coverageThreshold),
                                                     trainSampling = as.numeric(opt$trainSampling))

report_time_and_mem("Training model.")
model_MS <- train_model(raw_to_training_data_result[["training"]])
saveRDS(raw_to_training_data_result[["prior"]], file = opt$prior, compress = TRUE)
saveRDS(model_MS, file = opt$model, compress = T)

report_time_and_mem("Generating model lattice.")
lattice <- generate_lattice(model_MS)
saveRDS(lattice, file = opt$modelLattice, compress = T)

report_time_and_mem("Finished.")