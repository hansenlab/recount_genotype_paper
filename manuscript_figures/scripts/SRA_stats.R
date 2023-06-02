library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())

if(!file.exists("../ready_to_plot/all_SRA.annotated.rds")) {
	library(vroom)
	metadata = read.csv("/dcs04/hansen/data/recount_genotype/pipeline/AggregateFiles/all_SRA.csv")
	result = lapply(metadata$genotypedSample, function(x) {
		cat(x, "\n")
		genotypedSample = vroom(x, show_col_types=F)
		n_SNPs_pass_coverage = length(which(genotypedSample$coverage >= 15))
		return(list(nSNPs = nrow(genotypedSample),
					n_SNPs_pass_coverage = n_SNPs_pass_coverage))
	})
	result = as.data.frame(do.call(rbind, result))
	metadata = cbind(metadata, result)
	metadata$nSNPs = unlist(metadata$nSNPs)
	metadata$n_SNPs_pass_coverage = unlist(metadata$n_SNPs_pass_coverage)
	saveRDS(metadata, file = "../ready_to_plot/all_SRA.annotated.rds", compress = T)
}

