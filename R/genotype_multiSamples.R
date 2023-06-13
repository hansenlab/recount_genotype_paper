library(optparse)

option_list <- list(
    make_option(c("--modelLattice"), type = "character",
    help = "Input: R object: Mean and variance model in lattice format."),
    make_option(c("--model"), type = "character",
    help = "Input: R object: Mean and variance model."),
    make_option(c("--pZg"), type = "character",
    help = "Input: R object: Estimation of Pi prior for model."),
    make_option(c("-r", "--ref"), type = "character",
    help = "Input: R object (data.table): SNP x samples data.table of
    ref counts."),
    make_option(c("-a", "--alt"), type = "character",
    help = "Input: R object (data.table): SNP x samples data.table of
    alt counts."),
    make_option(c("-g", "--trueGenotype"), type = "character", help = "Input:
    R object (data.table) SNP x samples of true genotypes."),
    make_option(c("-p", "--predictedGenotype"), type = "character", help = "Output:
    Predicted genotypes."))

opt <- parse_args(OptionParser(option_list = option_list))

if (length(opt$model) == 0 |
    length(opt$modelLattice) == 0 | length(opt$pZg) == 0 |
    length(opt$alt) == 0 | length(opt$ref) == 0 |
    length(opt$predictedGenotype) == 0 |
    length(opt$trueGenotype) == 0) {
    stop("Not all arguments provided. Check --help for description.")
}

library(dplyr)
library(stats4)
library(mvtnorm)
library(caret)
library(parallel)
library(data.table)

ptm_all <- proc.time()

model_MS_lattice <- readRDS(opt$modelLattice)
model_MS <- readRDS(opt$model)

#load in genotype priors. If this is NA, we use the true genotype priors.
p_Zg <- readRDS(opt$pZg)
oraclePrior <- F
trueGenotype <- NA
if(all(is.na(p_Zg)) == T) {
  oraclePrior <- T
  trueGenotype <- readRDS(opt$trueGenotype) 
}

ref_matrix <- readRDS(opt$ref)
alt_matrix <- readRDS(opt$alt)
sample_names <- colnames(ref_matrix)

ref_matrix <- ref_matrix + 1  #to prevent divde by 0 in M, S calculation
alt_matrix <- alt_matrix + 1

M <- log2(ref_matrix / alt_matrix)
S <- log2(sqrt(ref_matrix * alt_matrix))

rm(ref_matrix)
rm(alt_matrix)

lattice_max <- 8.5 #hardcoded max value for S. 
#S[S >= lattice_max] <- lattice_max


predicted_genotype = as.data.table(matrix(NA, nrow = nrow(M), 
                                              ncol = ncol(M)))

colnames(predicted_genotype) <- sample_names

for(sample_i in seq_len(ncol(predicted_genotype))) {

  cat(sample_i, " out of ", ncol(predicted_genotype), "\n")
  ptm <- proc.time()

  if(oraclePrior == T) {
    genotype_i <- unlist(trueGenotype[, ..sample_i])
    if(all(!is.na(genotype_i))) { 
      p_Zg <- c(sum(genotype_i == "0/0") / length(genotype_i), #0 is reference homogyzous
                sum(genotype_i == "0/1") / length(genotype_i), #1 is het
                sum(genotype_i == "1/1") / length(genotype_i)) #2 is alt homozygous
      cat("Oracle prior:\n")
      print(p_Zg)
    }else { #if there is no true genotype, use hardcoded prior.
      cat("Using hardcoded prior.\n")
      p_Zg = c(.94, .04, .02) #hardcoded prior when there isn't a true genotype sample. This is a reasonable one based on what we see in genotype data.
    }
  }
  #calculate posterior probabilites
  #pi_g * N(M, sigma_g)

  withinLattice <- which(unlist(S[, ..sample_i]) <= lattice_max)
  outsideLattice <- which(unlist(S[, ..sample_i]) > lattice_max)

  p_Zg_dens <- sapply(1:3, function(k){

    #use our lattice look-up table to get our prediction
    mu_M <- rep(NA, nrow(S))
    sd_M <- rep(NA, nrow(S))
    if(k == 1) {
      mu_M[withinLattice] <- model_MS_lattice[J(S[withinLattice, ..sample_i])]$mu_M_prediction_1
      sd_M[withinLattice] <- model_MS_lattice[J(S[withinLattice, ..sample_i])]$sd_M_prediction_1
      mu_M[outsideLattice] <- predict(model_MS[[1]][[1]], data.frame(mu_S = unlist(S[outsideLattice, ..sample_i])))
      sd_M[outsideLattice] <- predict(model_MS[[1]][[2]], data.frame(mu_S = unlist(S[outsideLattice, ..sample_i])))
    }else if(k == 2) {
      mu_M[withinLattice] <- model_MS_lattice[J(S[withinLattice, ..sample_i])]$mu_M_prediction_2
      sd_M[withinLattice] <- model_MS_lattice[J(S[withinLattice, ..sample_i])]$sd_M_prediction_2
      mu_M[outsideLattice] <- predict(model_MS[[2]][[1]], data.frame(mu_S = unlist(S[outsideLattice, ..sample_i])))
      sd_M[outsideLattice] <- predict(model_MS[[2]][[2]], data.frame(mu_S = unlist(S[outsideLattice, ..sample_i])))
    }else if(k == 3) {
      mu_M[withinLattice] <- model_MS_lattice[J(S[withinLattice, ..sample_i])]$mu_M_prediction_3
      sd_M[withinLattice] <- model_MS_lattice[J(S[withinLattice, ..sample_i])]$sd_M_prediction_3
      mu_M[outsideLattice] <- predict(model_MS[[3]][[1]], data.frame(mu_S = unlist(S[outsideLattice, ..sample_i])))
      sd_M[outsideLattice] <- predict(model_MS[[3]][[2]], data.frame(mu_S = unlist(S[outsideLattice, ..sample_i])))
    }
    sd_M[sd_M <= 0] = 0 #if we predict (extrapolate) any variance to be <= 0, set it to 0. 

    dens <- p_Zg[k] * dnorm(x = unlist(M[, ..sample_i]), mean = mu_M, sd = sd_M, log = FALSE)
    
    return(dens)

  })


  #p(Z=g|M,S)
  p_Zg_MS <- sapply(1:3, function(k){
    p_Zg_dens[, k] / rowSums(p_Zg_dens)
  })


  inference <- lapply(seq_len(nrow(predicted_genotype)), function(snp_j) {
    high_prob_id <- which(p_Zg_MS[snp_j ,] == max(p_Zg_MS[snp_j ,]))
    if (high_prob_id == 1) geno <- "0/0"
    if (high_prob_id == 2) geno <- "0/1"
    else if (high_prob_id == 3) geno <- "1/1"
    return(list(geno,
                1 - p_Zg_MS[snp_j, high_prob_id]))
  })

  predicted_genotype[, (sample_i)] <- unlist(lapply(inference, function(x) x[1]))
   
  print(proc.time() - ptm)
  print(gc())
}

saveRDS(predicted_genotype, opt$predictedGenotype, compress = TRUE)


print(proc.time() - ptm_all)
