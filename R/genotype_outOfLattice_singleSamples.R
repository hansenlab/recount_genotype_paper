library(optparse)


option_list <- list(
  make_option(c("--genotypedSamples"), type = "character",
              help = "Input: List of genotyped samples."),
  make_option(c("--outOfLatticeSamples"), type = "character",
              help = "Input: List of out of lattice samples."),
  make_option(c("--model"), type = "character",
              help = "Input: R object: Mean and variance model."),
  make_option(c("--prior"), type = "character",
              help = "Input: R object: Estimation of Pi prior for model."),
  make_option(c("--genotypedSamplesComplete"), type = "character", 
              help = "Output: list of genotyped files to write back into with out of lattice SNPs.")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (length(opt$genotypedSamples) == 0 |
    length(opt$outOfLatticeSamples) == 0 |
    length(opt$model) == 0 |
    length(opt$prior) == 0 |
    length(opt$genotypedSamplesComplete) == 0) {
  stop("Not all arguments provided. Check --help for description.")
}

library(data.table)
library(tidyverse)
library(caret)

genotypedSamples_path <- opt$genotypedSamples
outOfLatticeSamples_path <- opt$outOfLatticeSamples



  cat("Loading in ", outOfLatticeSamples_path, "\n")
  all_samples <- read.csv(outOfLatticeSamples_path)
  
if(nrow(all_samples) > 0) {


all_samples <- all_samples %>% filter(chr != "chrX" & chr != "")

if(nrow(all_samples) > 0) {
  cat("Loading in model.\n")
  model_MS <- readRDS(opt$model)
  prior <- readRDS(opt$prior)
  if(all(is.na(prior)) == T) {
    prior = c(.94, .04, .02) 
  }
  
  cat("Performing inference.\n")
  likelihood_times_prior <- sapply(1:3, function(k){
    #Inference on mu_M and sd_M
    if(k == 1) {
      mu_M <- predict(model_MS[[1]][[1]], data.frame(mu_S = unlist(all_samples$S)))
      sd_M <- predict(model_MS[[1]][[2]], data.frame(mu_S = unlist(all_samples$S)))
    }else if(k == 2) {
      mu_M <- predict(model_MS[[2]][[1]], data.frame(mu_S = unlist(all_samples$S)))
      sd_M <- predict(model_MS[[2]][[2]], data.frame(mu_S = unlist(all_samples$S)))
    }else if(k == 3) {
      mu_M <- predict(model_MS[[3]][[1]], data.frame(mu_S = unlist(all_samples$S)))
      sd_M <- predict(model_MS[[3]][[2]], data.frame(mu_S = unlist(all_samples$S)))
    }
    sd_M[sd_M <= 0] = 0 #if we predict (extrapolate) any variance to be <= 0, set it to 0. 
    return(prior[k] * dnorm(x = all_samples$M, mean = mu_M, sd = sd_M, log = FALSE))
  })
  
  #Compute posterior probability of genotype.
  posterior <- sapply(1:3, function(k){
    if(is.null(dim(likelihood_times_prior))) { #if we only have one SNP that needs to be genotyped
      return(likelihood_times_prior[k] / sum(likelihood_times_prior))
    }else { #else, we have many SNPs in a dataframe
      return(likelihood_times_prior[, k] / rowSums(likelihood_times_prior))
    }
  })
  
  #Pick genotype with the max posterior genotype as our prediction.
  if(is.null(dim(posterior))) { #if we only have one SNP that needs to be genotyped
    all_samples$pred_genotype[1] <- which(posterior == max(posterior))
  }else { #else, we have many SNPs in a dataframe
    all_samples$pred_genotype <- apply(posterior, MARGIN = 1, FUN = function(x) which(x == max(x)))
  }
}
}
#Two options of output: write to an output file,
# if(length(opt$result) > 0) {
#   fwrite(all_samples, opt$result)
# }

#and/or write back to a complete genotyped sample file.
if(length(opt$genotypedSamplesComplete) > 0) {
  genotypedSamplesComplete_path <- opt$genotypedSamplesComplete
  

      sample_i_lattice <- all_samples
      #read in original
      sample_i_original <- read.csv(genotypedSamples_path)
      #bind original with lattice
      sample_i <- rbind(sample_i_original, sample_i_lattice)
      #write to complete genotyped file
      cat("Writing to ", genotypedSamplesComplete_path, "\n")
      fwrite(sample_i, genotypedSamplesComplete_path)
  }


