library(optparse)

option_list <- list(
    make_option(c("-a", "--allGenotypedSamples"), type = "character",
    help = "Input: List of predicted genotyped files."),
    make_option(c("-s", "--allSampleIDRep"), type = "character",
    help = "Input: List of sample ID reps."),
    make_option(c("-u", "--outOfLatticeGenotypedSamples"), type = "character",
    help = "Input: R object (data.table): out of lattice genotyped file."),
    make_option(c("-c", "--commonSNPs"), type = "character",
    help = "Input: tab delim file: chr, start: aggregrate SNPs that intersected with BigWigs."),
    make_option(c("-g", "--genotypes"), type = "character",
    help = "Input: R object (data.table): location x sample dataframe of
    genotype (filtered version of genotype matrix)."),
    make_option(c("-o", "--conditionalAccuracy"), type = "character",
    help = "Output R object: conditional accuracy, averaged across every SNP."),
    make_option(c("--allGenotypedSamplesAgg"), type = "character",
    help = "Output R object: predicted genotyped files with true genotype information, in one data.table.")
)



opt <- parse_args(OptionParser(option_list = option_list))

if (length(opt$allGenotypedSamples) == 0 | length(opt$allSampleIDRep) == 0 |
    length(opt$commonSNPs) == 0 | length(opt$genotypes) == 0 | length(opt$outOfLatticeGenotypedSamples) == 0 |
    length(opt$conditionalAccuracy) == 0 | length(opt$allGenotypedSamplesAgg) == 0) {
  stop("Not all arguments provided. Check --help for description.")
}

library(tidyverse)
library(data.table)

compute_conditional_accuracy <- function(dt, major_allele_is_reference) {
	if(major_allele_is_reference) {
		dt <- dt %>% filter(AF <= .5)
	}else {
		dt <- dt %>% filter(AF > .5)
	}
	#full description of condtional accuracy is described in Table 1 in Methods of our preprint.
	acc <- dt %>% mutate(
                    accuracy1 = case_when(
                      true_genotype == 1 & pred_genotype == 1 ~ 1,
                      true_genotype == 1 & pred_genotype == 2 ~ .5,
                      true_genotype == 1 & pred_genotype == 3 ~ 0,
                      true_genotype == 2 & pred_genotype == 1 ~ 1,
                      true_genotype == 2 & pred_genotype == 2 ~ 1,
                      true_genotype == 2 & pred_genotype == 3 ~ 0,
                      true_genotype == 3 & pred_genotype == 1 ~ NA,
                      true_genotype == 3 & pred_genotype == 2 ~ NA,
                      true_genotype == 3 & pred_genotype == 3 ~ NA,
                    ),
                    accuracy2 = case_when(
                      true_genotype == 1 & pred_genotype == 1 ~ NA,
                      true_genotype == 1 & pred_genotype == 2 ~ NA,
                      true_genotype == 1 & pred_genotype == 3 ~ NA,
                      true_genotype == 2 & pred_genotype == 1 ~ 0,
                      true_genotype == 2 & pred_genotype == 2 ~ 1,
                      true_genotype == 2 & pred_genotype == 3 ~ 1,
                      true_genotype == 3 & pred_genotype == 1 ~ 0,
                      true_genotype == 3 & pred_genotype == 2 ~ .5,
                      true_genotype == 3 & pred_genotype == 3 ~ 1,
                    ),
                    overall_accuracy = case_when(
                      true_genotype == 1 & pred_genotype == 1 ~ 1,
                      true_genotype == 1 & pred_genotype == 2 ~ .5,
                      true_genotype == 1 & pred_genotype == 3 ~ 0,
                      true_genotype == 2 & pred_genotype == 1 ~ .5,
                      true_genotype == 2 & pred_genotype == 2 ~ 1,
                      true_genotype == 2 & pred_genotype == 3 ~ .5,
                      true_genotype == 3 & pred_genotype == 1 ~ 0,
                      true_genotype == 3 & pred_genotype == 2 ~ .5,
                      true_genotype == 3 & pred_genotype == 3 ~ 1,
                    )
                ) %>%
                summarise(sum_overall_accuracy = sum(overall_accuracy),
                          sum_accuracy_1 = sum(accuracy1, na.rm=T),
                          sum_accuracy_2 = sum(accuracy2, na.rm=T),
                          n = n(),
                          n1 = sum(!is.na(accuracy1)),
                          n2 = sum(!is.na(accuracy2))
                          )
    
  major_allele_is_reference_columns <- c("sum_overall_accuracy", "sum_major_accuracy", "sum_minor_accuracy", "n", "n1", "n2")
	not_major_allele_is_reference_columns <- c("sum_overall_accuracy", "sum_minor_accuracy", "sum_major_accuracy", "n", "n2", "n1")
    
  if(major_allele_is_reference) {
		names(acc) <- major_allele_is_reference_columns
	}else {
		names(acc) <- not_major_allele_is_reference_columns
		acc <- acc[, match(not_major_allele_is_reference_columns, major_allele_is_reference_columns)] #reorder columns
	}
	return(acc)
}


genotypes <- readRDS(opt$genotypes)
genotypes_chr_pos <- paste0(genotypes$CHROM, "_", genotypes$POS)

agg_snps <- fread(opt$commonSNPs, header = F)
agg_snps <- agg_snps[agg_snps$V1 != "" ,] #trailing whitespace issue from empty out of lattice entries - fixed upstream also.

genotypedSamples_path <- scan(opt$allGenotypedSamples, what = "character")
sample_id_rep <- scan(opt$allSampleIDRep, what = "character")

found_idx <- which(sample_id_rep %in% colnames(genotypes))
genotypedSamples_path <- genotypedSamples_path[found_idx]
sample_id_rep <- sample_id_rep[found_idx]

#pull in true genotype for lattice genotyped samples.
all_samples = list()
for(i in 1:length(genotypedSamples_path)) {
	cat(i, " out of ", length(genotypedSamples_path), ": Loading in ", genotypedSamples_path[i], "\n")
	sample_i <- fread(genotypedSamples_path[i])
	sample_i <- sample_i %>% filter(chr != "chrX" & chr != "")
	sample_id_rep_i <- sample_id_rep[i]
	match_idx <- match(paste0(sample_i$chr, "_", sample_i$start), genotypes_chr_pos)
	sample_i <- sample_i[which(!is.na(match_idx)) ,]
	match_idx <- match_idx[which(!is.na(match_idx))]
	sample_i$true_genotype <- unlist(genotypes[match_idx, ..sample_id_rep_i])
	sample_i <- sample_i %>% filter(true_genotype != "./.")
	if(length(which(is.na(sample_i$true_genotype))) > 0) {
		warning("Some of ", genotypedSamples_path[i], " SNPs contain NA for true genotype.\n")
	}
	sample_i$sample_id_rep <- sample_id_rep_i
	all_samples[[i]] <- sample_i
}
all_samples <- do.call(rbind, all_samples)
print(gc())

#pull in true genotype for non-lattice genotyped samples.
outOfLattice_samples <- fread(opt$outOfLatticeGenotypedSamples)
outOfLattice_samples <- outOfLattice_samples[outOfLattice_samples$sample_id_rep %in% colnames(genotypes) ,]
outOfLattice_samples <- outOfLattice_samples %>% filter(chr != "chrX" & chr != "")
outOfLattice_samples_updated <- list()
for(i in 1:length(unique(outOfLattice_samples$sample_id_rep))) {
	sample_id_rep_i <- unique(outOfLattice_samples$sample_id_rep)[i]
	cat(i, " out of ", length(unique(outOfLattice_samples$sample_id_rep)), ": Processing ", sample_id_rep_i, "\n")
	outOfLattice_samples_i <- outOfLattice_samples %>% filter(sample_id_rep == sample_id_rep_i)
	match_idx <- match(paste0(outOfLattice_samples_i$chr, "_", outOfLattice_samples_i$start), genotypes_chr_pos)
	outOfLattice_samples_i <- outOfLattice_samples_i[which(!is.na(match_idx)) ,]
	match_idx <- match_idx[which(!is.na(match_idx))]
	outOfLattice_samples_i$true_genotype <- unlist(genotypes[match_idx, ..sample_id_rep_i])
	outOfLattice_samples_i <- outOfLattice_samples_i %>% filter(true_genotype != "./.")
	if(length(which(is.na(outOfLattice_samples_i$true_genotype))) > 0) {
		warning("Some of ", sample_id_rep_i, " SNPs contain NA for true genotype.\n")
	}
	outOfLattice_samples_updated[[i]] <- outOfLattice_samples_i
}
outOfLattice_samples_updated <- do.call(rbind, outOfLattice_samples_updated)

rm(genotypes) #clear large genotype matrix from working memotry.


#put the two together.
all_samples <- rbind(all_samples, outOfLattice_samples_updated)


#clean up varaiables.
all_samples$pred_genotype <- as.character(all_samples$pred_genotype)
all_samples$true_genotype[all_samples$true_genotype %in% c("0/0", "0|0")] <- "1"
all_samples$true_genotype[all_samples$true_genotype %in% c("0/1", "1/0", "0|1", "1|0")] <- "2"
all_samples$true_genotype[all_samples$true_genotype %in% c("1/1", "1|1")] <- "3"
all_samples <- all_samples %>% filter(true_genotype != "." & true_genotype != "./.")


#compute conditional accuracy
acc1 <- compute_conditional_accuracy(all_samples, major_allele_is_reference = T)
acc2 <- compute_conditional_accuracy(all_samples, major_allele_is_reference = F)
result <- rbind(acc1, acc2)
result <- data.frame(major_accuracy = sum(result$sum_major_accuracy)/sum(result$n1),
                    minor_accuracy = sum(result$sum_minor_accuracy)/sum(result$n2),
                    overall_accuracy = sum(result$sum_overall_accuracy)/sum(result$n))


#save.
saveRDS(result, opt$conditionalAccuracy, compress = TRUE)
saveRDS(all_samples, opt$allGenotypedSamplesAgg, compress = TRUE)
