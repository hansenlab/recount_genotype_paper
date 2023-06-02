library(optparse)

option_list <- list(
  make_option(c("-f", "--filteredSNPs"), type = "character",
              help = "Input: Filtered biallelic SNPs in .rds format."),
  make_option(c("-r", "--ref"), type = "character",
              help = "Input: Reference matrix as a data.table in .rds format"),
  make_option(c("-a", "--alt"), type = "character",
              help = "Input: Alternative matrix as a data.table in .rds format"),
  make_option(c("-g", "--trueGenotype"), type = "character",
              help = "Input: Ground truth genotype matrix as a data.table in .rds format"),
  make_option(c("-p", "--predictedGenotype"), type = "character",
              help = "Input: Predicted model genotype matrix as a data.table in .rds format"),
  make_option(c("-m", "--model"), type = "character",
              help = "Output: Accuracy model in .rds format."),
    make_option(c("-l", "--lattice"), type = "character",
              help = "Output: Accuracy model lattice in .rds format."),
  make_option(c("-c", "--coverageThreshold"), type = "numeric", default = 3,
              help = "Option: Coverage threshold to filter. Default: 3"),
  make_option(c("-t", "--trainSampling"), type = "numeric", default = 1,
              help = "Option: Sampling proportion for training. Default: 1 (no sampling)."))

opt <- parse_args(OptionParser(option_list = option_list))

if(length(opt$filteredSNPs) == 0 | length(opt$ref) == 0 |
    length(opt$alt) == 0 | length(opt$trueGenotype) == 0 |
    length(opt$predictedGenotype) == 0 | length(opt$coverageThreshold) == 0 |
    length(opt$model) == 0 | length(opt$lattice) == 0) {
  stop("Not all arguments provided. Check --help for description.")
}

library(tidyverse)
library(data.table)
library(mgcv) #GAMs
library(rms) #restricted cubic splines
library(GenomicRanges)


report_time_and_mem <- function(msg) {
  cat(paste0(msg, "\n"))
  print(Sys.time())
  print(gc())
}

expit <- function(x) {
  return(1/(1+exp(-x)))
}

merge_raw_data_to_long <- function(snps, ref, alt, pred_genotype, true_genotype) {
  stopifnot(identical(dim(ref), dim(alt)))
  stopifnot(identical(dim(ref), dim(pred_genotype)))
  stopifnot(identical(dim(ref), dim(true_genotype)))

  coverage <- ref + alt
  rm(ref, alt)

  #create ID based on chr and start position
  chr_start <- paste0(as.character(seqnames(snps)), "_", start(snps))
  true_genotype$chr_start <- chr_start
  pred_genotype$chr_start <- chr_start
  coverage$chr_start <- chr_start

  #convert from wide to long
  true_genotype <- true_genotype %>% pivot_longer(!chr_start, names_to = "sample_id_rep", values_to = "true_genotype")
  pred_genotype <- pred_genotype %>% pivot_longer(!chr_start, names_to = "sample_id_rep", values_to = "pred_genotype")
  coverage <- coverage %>% pivot_longer(!chr_start, names_to = "sample_id_rep", values_to = "coverage")

  #merge
  training_data <- inner_join(coverage, inner_join(true_genotype, pred_genotype))
  
  #add AF
  training_data$AF <- snps$allele_freq[match(training_data$chr_start, chr_start)]

  return(training_data)
}

raw_to_training_data <- function(snps, ref, alt, pred_genotype, true_genotype, coverageThreshold, trainSampling) {
  
  training_data <- merge_raw_data_to_long(snps, ref, alt, pred_genotype, true_genotype)

  #filter based on coverage and non-missing true genotype
  training_data <- training_data %>% 
    filter(coverage >= coverageThreshold) %>% 
    filter(true_genotype != "./." & true_genotype != ".") 

  #indicate whether the predicted genotype is correct or not
  training_data$correct <- as.numeric(training_data$pred_genotype == training_data$true_genotype)

  #compute major allele fraction
  training_data$major_AF <- case_when(training_data$AF <= .5  ~ 1 - training_data$AF,
                                     training_data$AF > .5  ~ training_data$AF)

  #split our training data into whether they have MAF >= .95 or not.
  training_high_majorAF <- training_data %>% filter(major_AF >= .95)
  training_low_majorAF <- training_data %>% filter(major_AF < .95)

  #downsample based on `trainSampling`.
  training_high_majorAF <- training_high_majorAF[sample(1:nrow(training_high_majorAF), trainSampling * nrow(training_high_majorAF)) ,]
  training_low_majorAF <- training_low_majorAF[sample(1:nrow(training_low_majorAF), trainSampling * nrow(training_low_majorAF)) ,]

  #create majorAF_bin
  low_majorAF_quantile <- quantile(training_low_majorAF$major_AF, seq(0, 1, by = 1/5))
  training_low_majorAF$majorAF_bin <- as.character(cut(training_low_majorAF$major_AF, low_majorAF_quantile, include.lowest=TRUE))
  training_high_majorAF$majorAF_bin <- "[0.95, 1]"

  return(list(low_majorAF = training_low_majorAF,
              high_majorAF = training_high_majorAF))
}

train_model <- function(training_data) {
  low_majorAF_accuracy_model <- rms::Glm(correct ~ rcs(coverage, knots = c(4, 6, 10, 16, 40)) * majorAF_bin, 
                                       family = "binomial", 
                                       data = training_data[["low_majorAF"]])

  high_majorAF_accuracy_model <- rms::Glm(correct ~ rcs(coverage, knots = c(4, 6, 10, 16, 40)), 
                                         family = "binomial", 
                                         data = training_data[["high_majorAF"]])

  return(list(low_majorAF = low_majorAF_accuracy_model,
              high_majorAF = high_majorAF_accuracy_model))
}

generate_lattice <- function(training_data, accuracy_model) {
  majorAF_bins <- c(unique(training_data[["low_majorAF"]]$majorAF_bin),
                    unique(training_data[["high_majorAF"]]$majorAF_bin))
  lattice <- expand.grid(majorAF_bin = majorAF_bins, coverage = seq(4, 100))

  lattice_low_majorAF <- lattice %>% filter(majorAF_bin != "[0.95, 1]")
  lattice_low_majorAF$predicted_accuracy <- expit(predict(accuracy_model[["low_majorAF"]], 
                                                          newdata = as.data.frame(lattice_low_majorAF)))
  
  lattice_high_majorAF <- lattice %>% filter(majorAF_bin == "[0.95, 1]")
  lattice_high_majorAF$predicted_accuracy <- expit(predict(accuracy_model[["high_majorAF"]], 
                                                           newdata = as.data.frame(lattice_high_majorAF)))

  #put the two together
  lattice <- rbind(lattice_low_majorAF, lattice_high_majorAF)
  return(lattice)
}


#######################
#Main code starts here.

report_time_and_mem("Loading and transforming data for training.")
training_data <- raw_to_training_data(snps = readRDS(opt$filteredSNPs),
                                      ref = readRDS(opt$ref),
                                      alt = readRDS(opt$alt),
                                      pred_genotype = readRDS(opt$predictedGenotype),
                                      true_genotype = readRDS(opt$trueGenotype),
                                      coverageThreshold = as.numeric(opt$coverageThreshold),
                                      trainSampling = as.numeric(opt$trainSampling))

report_time_and_mem("Training accuracy model.")
accuracy_model <- train_model(training_data)

report_time_and_mem("Generating lattice.")
accuracy_lattice <- generate_lattice(training_data, accuracy_model)

saveRDS(accuracy_model, file = opt$model)
saveRDS(accuracy_lattice, file = opt$lattice)
report_time_and_mem("Finished.")
