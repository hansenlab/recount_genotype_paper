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
    make_option(c("-d", "--dbFile"), type = "character",
    help = "Output SQLite .db file: predicted genotypes and true genotypes in 
    long format for SNPs and samples in long format."),
    make_option(c("-o", "--conditionalAccuracy"), type = "character",
    help = "Output R object: conditional accuracy"),
    make_option(c("-m", "--majorAlleleModelConditionalAccuracy"), type = "character",
    help = "Output R object: conditional accuracy based on major allele model")
)

opt <- parse_args(OptionParser(option_list = option_list))


if (length(opt$allGenotypedSamples) == 0 | length(opt$allSampleIDRep) == 0 |
    length(opt$commonSNPs) == 0 | length(opt$genotypes) == 0 |
    length(opt$dbFile) == 0 | length(opt$outOfLatticeGenotypedSamples) == 0 |
    length(opt$conditionalAccuracy) == 0 | length(opt$majorAlleleModelConditionalAccuracy) == 0) {
  stop("Not all arguments provided. Check --help for description.")
}

library(tidyverse)
library(data.table)
library(DBI)
library(RSQLite)

if(!file.exists(opt$dbFile)) {
	cat("Loading genotypes...\n")
	genotypes <- readRDS(opt$genotypes)
	genotypes_chr_pos <- paste0(genotypes$CHROM, "_", genotypes$POS)

	cat("Loading SNPs and sample IDs...\n")
	agg_snps <- fread(opt$commonSNPs, header = F)
	agg_snps <- agg_snps[agg_snps$V1 != "" ,] #trailing whitespace issue from empty out of lattice entries - fixed upstream also.

	genotypedSamples_path <- scan(opt$allGenotypedSamples, what = "character")
	sample_id_rep <- scan(opt$allSampleIDRep, what = "character")

	found_idx <- which(sample_id_rep %in% colnames(genotypes))
	genotypedSamples_path <- genotypedSamples_path[found_idx]
	sample_id_rep <- sample_id_rep[found_idx]

	all_predGenotypes_w_trueGenotypes_col <- c("chr", "start", "AF", "M", "S", "coverage", "pred_genotype", "true_genotype", "sample_id_rep")

	cat("Connecting to SQLite database...\n")
	con <- dbConnect(RSQLite::SQLite(), opt$dbFile)

	#pull in true genotype for lattice genotyped samples.
	for(i in 1:length(genotypedSamples_path)) {
		cat(i, " out of ", length(genotypedSamples_path), ": Loading in ", genotypedSamples_path[i], "\n")
		print(gc())
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

		sample_i$pred_genotype <- as.character(sample_i$pred_genotype)
		sample_i$true_genotype[sample_i$true_genotype %in% c("0/0", "0|0")] <- "1"
		sample_i$true_genotype[sample_i$true_genotype %in% c("0/1", "1/0", "0|1", "1|0")] <- "2"
		sample_i$true_genotype[sample_i$true_genotype %in% c("1/1", "1|1")] <- "3"
		sample_i <- sample_i[sample_i$true_genotype != "./." ,]
		sample_i <- sample_i[sample_i$true_genotype != "." ,]

		names(sample_i) <- all_predGenotypes_w_trueGenotypes_col
		if(i == 1) {
			dbWriteTable(con, "all_predGenotypes_w_trueGenotypes", sample_i, overwrite = TRUE, append = FALSE)
		}else {
			dbWriteTable(con, "all_predGenotypes_w_trueGenotypes", sample_i, append = TRUE)
		}
	}
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

		outOfLattice_samples_i$pred_genotype <- as.character(outOfLattice_samples_i$pred_genotype)
		outOfLattice_samples_i$true_genotype[outOfLattice_samples_i$true_genotype %in% c("0/0", "0|0")] <- "1"
		outOfLattice_samples_i$true_genotype[outOfLattice_samples_i$true_genotype %in% c("0/1", "1/0", "0|1", "1|0")] <- "2"
		outOfLattice_samples_i$true_genotype[outOfLattice_samples_i$true_genotype %in% c("1/1", "1|1")] <- "3"
		outOfLattice_samples_i <- outOfLattice_samples_i[outOfLattice_samples_i$true_genotype != "./." ,]
		outOfLattice_samples_i <- outOfLattice_samples_i[outOfLattice_samples_i$true_genotype != "." ,]
		outOfLattice_samples_updated[[i]] <- outOfLattice_samples_i
	}
	outOfLattice_samples_updated <- do.call(rbind, outOfLattice_samples_updated)
	names(outOfLattice_samples_updated) <- all_predGenotypes_w_trueGenotypes_col
	dbWriteTable(con, "all_predGenotypes_w_trueGenotypes", outOfLattice_samples_updated, append = TRUE)

	cat("Done adding in entires. Creating index. \n")
	dbSendQuery(con, "CREATE INDEX chrposIdx ON all_predGenotypes_w_trueGenotypes (chr, start, sample_id_rep);")

	dbDisconnect(con)
	rm(genotypes) 
}

print(gc())

cat("Starting to compute aconditional accuracy by querying SQLite database.\n")

con <- DBI::dbConnect(RSQLite::SQLite(), opt$dbFile)
geno_db <- tbl(con, "all_predGenotypes_w_trueGenotypes")


cat("Computing major accuracy.\n")

accuracy1 = geno_db %>% 
	filter(AF <= .5) %>%
	mutate(
		major_accuracy = ifelse(true_genotype == "1" & pred_genotype == "1", 1, 0) +
						 ifelse(true_genotype == "1" & pred_genotype == "2", .5, 0) +
						 ifelse(true_genotype == "1" & pred_genotype == "3", 0, 0) +
						 ifelse(true_genotype == "2" & pred_genotype == "1", 1, 0) +
						 ifelse(true_genotype == "2" & pred_genotype == "2", 1, 0) +
						 ifelse(true_genotype == "2" & pred_genotype == "3", 0, 0) +
						 ifelse(true_genotype == "3" & pred_genotype == "1", NA, 0) +
						 ifelse(true_genotype == "3" & pred_genotype == "2", NA, 0) +
						 ifelse(true_genotype == "3" & pred_genotype == "3", NA, 0),
		minor_accuracy = ifelse(true_genotype == "1" & pred_genotype == "1", NA, 0) +
						 ifelse(true_genotype == "1" & pred_genotype == "2", NA, 0) +
						 ifelse(true_genotype == "1" & pred_genotype == "3", NA, 0) +
						 ifelse(true_genotype == "2" & pred_genotype == "1", 0, 0) +
						 ifelse(true_genotype == "2" & pred_genotype == "2", 1, 0) +
						 ifelse(true_genotype == "2" & pred_genotype == "3", 1, 0) +
						 ifelse(true_genotype == "3" & pred_genotype == "1", 0, 0) +
						 ifelse(true_genotype == "3" & pred_genotype == "2", .5, 0) +
						 ifelse(true_genotype == "3" & pred_genotype == "3", 1, 0),
		overall_accuracy = ifelse(true_genotype == "1" & pred_genotype == "1", 1, 0) +
						   ifelse(true_genotype == "1" & pred_genotype == "2", .5, 0) +
						   ifelse(true_genotype == "1" & pred_genotype == "3", 0, 0) +
						   ifelse(true_genotype == "2" & pred_genotype == "1", .5, 0) +
						   ifelse(true_genotype == "2" & pred_genotype == "2", 1, 0) +
						   ifelse(true_genotype == "2" & pred_genotype == "3", .5, 0) +
						   ifelse(true_genotype == "3" & pred_genotype == "1", 0, 0) +
						   ifelse(true_genotype == "3" & pred_genotype == "2", .5, 0) +
						   ifelse(true_genotype == "3" & pred_genotype == "3", 1, 0)) %>%
	summarise(sum_overall_accuracy = sum(overall_accuracy),
              sum_major_accuracy = sum(major_accuracy, na.rm=T),
              sum_minor_accuracy = sum(minor_accuracy, na.rm=T),
              n = n(),
              n_major = sum(!is.na(major_accuracy)),
              n_minor = sum(!is.na(minor_accuracy))
              ) %>%
	collect()


cat("Computing minor accuracy.\n")

accuracy2 = geno_db %>%
	filter(AF > .5) %>%
	mutate(
		minor_accuracy = ifelse(true_genotype == "1" & pred_genotype == "1", 1, 0) +
						 ifelse(true_genotype == "1" & pred_genotype == "2", .5, 0) +
						 ifelse(true_genotype == "1" & pred_genotype == "3", 0, 0) +
						 ifelse(true_genotype == "2" & pred_genotype == "1", 1, 0) +
						 ifelse(true_genotype == "2" & pred_genotype == "2", 1, 0) +
						 ifelse(true_genotype == "2" & pred_genotype == "3", 0, 0) +
						 ifelse(true_genotype == "3" & pred_genotype == "1", NA, 0) +
						 ifelse(true_genotype == "3" & pred_genotype == "2", NA, 0) +
						 ifelse(true_genotype == "3" & pred_genotype == "3", NA, 0),
		major_accuracy = ifelse(true_genotype == "1" & pred_genotype == "1", NA, 0) +
						 ifelse(true_genotype == "1" & pred_genotype == "2", NA, 0) +
						 ifelse(true_genotype == "1" & pred_genotype == "3", NA, 0) +
						 ifelse(true_genotype == "2" & pred_genotype == "1", 0, 0) +
						 ifelse(true_genotype == "2" & pred_genotype == "2", 1, 0) +
						 ifelse(true_genotype == "2" & pred_genotype == "3", 1, 0) +
						 ifelse(true_genotype == "3" & pred_genotype == "1", 0, 0) +
						 ifelse(true_genotype == "3" & pred_genotype == "2", .5, 0) +
						 ifelse(true_genotype == "3" & pred_genotype == "3", 1, 0),
		overall_accuracy = ifelse(true_genotype == "1" & pred_genotype == "1", 1, 0) +
						   ifelse(true_genotype == "1" & pred_genotype == "2", .5, 0) +
						   ifelse(true_genotype == "1" & pred_genotype == "3", 0, 0) +
						   ifelse(true_genotype == "2" & pred_genotype == "1", .5, 0) +
						   ifelse(true_genotype == "2" & pred_genotype == "2", 1, 0) +
						   ifelse(true_genotype == "2" & pred_genotype == "3", .5, 0) +
						   ifelse(true_genotype == "3" & pred_genotype == "1", 0, 0) +
						   ifelse(true_genotype == "3" & pred_genotype == "2", .5, 0) +
						   ifelse(true_genotype == "3" & pred_genotype == "3", 1, 0)) %>%
	summarise(sum_overall_accuracy = sum(overall_accuracy),
              sum_major_accuracy = sum(major_accuracy, na.rm=T),
              sum_minor_accuracy = sum(minor_accuracy, na.rm=T),
              n = n(),
              n_major = sum(!is.na(major_accuracy)),
              n_minor = sum(!is.na(minor_accuracy))
              ) %>%
	collect()


result = rbind(accuracy1, accuracy2)
result = data.frame(major_accuracy = sum(result$sum_major_accuracy)/sum(result$n_major),
                    minor_accuracy = sum(result$sum_minor_accuracy)/sum(result$n_minor),
                    overall_accuracy = sum(result$sum_overall_accuracy)/sum(result$n))

saveRDS(result, opt$conditionalAccuracy)


#####

cat("Computing major alelle model major accuracy.\n")

null_accuracy1 = geno_db %>% 
	filter(AF <= .5) %>%
	mutate(
		major_accuracy = ifelse(true_genotype == "1", 1, 0) +
						 ifelse(true_genotype == "2", 1, 0) +
						 ifelse(true_genotype == "3", NA, 0),
		minor_accuracy = ifelse(true_genotype == "1", NA, 0) +
						 ifelse(true_genotype == "2", 0, 0) +
						 ifelse(true_genotype == "3", 0, 0),
		overall_accuracy = ifelse(true_genotype == "1", 1, 0) +
						   ifelse(true_genotype == "2", .5, 0) +
						   ifelse(true_genotype == "3", 0, 0))  %>%
	summarise(sum_overall_accuracy = sum(overall_accuracy),
              sum_major_accuracy = sum(major_accuracy, na.rm=T),
              sum_minor_accuracy = sum(minor_accuracy, na.rm=T),
              n = n(),
              n_major = sum(!is.na(major_accuracy)),
              n_minor = sum(!is.na(minor_accuracy))
              ) %>%
	collect()


cat("Computing major allele model minor accuracy.\n")
null_accuracy2 = geno_db %>% 
	filter(AF > .5) %>%
	mutate(
		major_accuracy = ifelse(true_genotype == "1", 1, 0) +
						 ifelse(true_genotype == "2", 1, 0) +
						 ifelse(true_genotype == "3", NA, 0),
		minor_accuracy = ifelse(true_genotype == "1", NA, 0) +
						 ifelse(true_genotype == "2", 0, 0) +
						 ifelse(true_genotype == "3", 0, 0),
		overall_accuracy = ifelse(true_genotype == "1", 1, 0) +
						   ifelse(true_genotype == "2", .5, 0) +
						   ifelse(true_genotype == "3", 0, 0))  %>%
	summarise(sum_overall_accuracy = sum(overall_accuracy),
              sum_major_accuracy = sum(major_accuracy, na.rm=T),
              sum_minor_accuracy = sum(minor_accuracy, na.rm=T),
              n = n(),
              n_major = sum(!is.na(major_accuracy)),
              n_minor = sum(!is.na(minor_accuracy))
              ) %>%
	collect()


null_result = rbind(null_accuracy1, null_accuracy2)
null_result = data.frame(major_accuracy = sum(null_result$sum_major_accuracy)/sum(null_result$n_major),
                    minor_accuracy = sum(null_result$sum_minor_accuracy)/sum(null_result$n_minor),
                    overall_accuracy = sum(null_result$sum_overall_accuracy)/sum(null_result$n))

saveRDS(null_result, opt$majorAlleleModelConditionalAccuracy)