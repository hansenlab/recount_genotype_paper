library(tidyverse)
library(cowplot)
library(data.table)
theme_set(theme_cowplot())


recount3_metadata<-fread("Recount3_metadata.tsv", header= T, sep = "\t",quote="") #metadata dowloaded from recount3 website

file_path<-"ready_to_plot/all_SRA.annotated.rds"
if(!file.exists(file_path)) {
	library(vroom)
	metadata = read.csv("all_SRA.csv") #path to SRA genotype files
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
metadata<-readRDS(file_path)
metadata$seq_type<-recount3_metadata$seq_type[match(metadata$sample_id,recount3_metadata$external_id)]


metadata %>% group_by(seq_type) %>% 
  summarize(n(),mean(nSNPs), mean(n_SNPs_pass_coverage/nSNPs))


