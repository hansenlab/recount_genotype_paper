library(optparse)

option_list <- list(
    make_option(c("--genotypedSamples"), type = "character",
    help = "Input: List of genotyped samples."),
    make_option(c("--outOfLatticeSamples"), type = "character",
    help = "Input: List of out of lattice samples."),
    make_option(c("-s", "--allSampleIDRep"), type = "character",
    help = "Input: List of sample ID reps."),
    make_option(c("--model"), type = "character",
    help = "Input: R object: Mean and variance model."),
    make_option(c("--prior"), type = "character",
    help = "Input: R object: Estimation of Pi prior for model."),
    make_option(c("--result"), type = "character", 
    help = "Output csv: sample_id_rep, chr, start, AF, M, S, pred_genotype"),
    make_option(c("--genotypedSamplesComplete"), type = "character", 
    help = "Output: list of genotyped files to write back into with out of lattice SNPs.")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (length(opt$genotypedSamples) == 0 |
    length(opt$outOfLatticeSamples) == 0 |
    length(opt$allSampleIDRep) == 0 |
	length(opt$model) == 0 |
	length(opt$prior) == 0 |
	(length(opt$result) == 0 & length(opt$genotypedSamplesComplete) == 0)) {
        stop("Not all arguments provided. Check --help for description.")
}

library(data.table)
library(tidyverse)
library(caret)

genotypedSamples_path <- scan(opt$genotypedSamples, what = "character")
outOfLatticeSamples_path <- scan(opt$outOfLatticeSamples, what = "character")
all_sample_id_rep <- scan(opt$allSampleIDRep, what = "character")
stopifnot(length(all_sample_id_rep) == length(outOfLatticeSamples_path))
stopifnot(length(all_sample_id_rep) == length(genotypedSamples_path))

all_samples <- list()
for(i in 1:length(outOfLatticeSamples_path)) {
    cat("Loading in ", outOfLatticeSamples_path[i], "\n")
    sample_i <- fread(outOfLatticeSamples_path[i])
    sample_i$sample_id_rep <- all_sample_id_rep[i]
    all_samples[[i]] <- sample_i
}
all_samples <- do.call(rbind, all_samples)
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

#Two options of output: write to an output file,
if(length(opt$result) > 0) {
    fwrite(all_samples, opt$result)
}

#and/or write back to a complete genotyped sample file.
if(length(opt$genotypedSamplesComplete) > 0) {
    genotypedSamplesComplete_path <- scan(opt$genotypedSamplesComplete, what = "character")
    stopifnot(length(genotypedSamplesComplete_path) == length(outOfLatticeSamples_path))
    for(i in 1:length(genotypedSamplesComplete_path)) {
        sample_i_lattice <- all_samples %>% filter(sample_id_rep == all_sample_id_rep[i])
        if(nrow(sample_i_lattice) > 0) {
            sample_i_lattice <- sample_i_lattice %>% select(!sample_id_rep)
            #read in original
            sample_i_original <- fread(genotypedSamples_path[i])
            #bind original with lattice
            sample_i <- rbind(sample_i_original, sample_i_lattice)
            #write to complete genotyped file
            cat("Writing to ", genotypedSamplesComplete_path[i], "\n")
            fwrite(sample_i, genotypedSamplesComplete_path[i])
        }
    }
}

