library(tidyverse)
library(data.table)
library(cowplot)
theme_set(theme_cowplot())

plot_file = "../ready_to_plot/GTEx_LCL_vs_Geuvadis_M_S.rds"

if(!file.exists(plot_file)) {
	GTEX_metadata = read.csv("/dcs04/hansen/data/recount_genotype/pipeline/AggregateFiles/GTEx_Blended_Tissue_Testing_metadata.csv")
	GTEx_lcl_snp = readRDS(GTEX_metadata$allGenotypesOutput[GTEX_metadata$study == "Cells_EBV-transformed_lymphocytes"])

	geuvadis_metadata = read.csv("/dcs04/hansen/data/recount_genotype/pipeline/AggregateFiles/Geuvadis_testing_metadata.csv")
	geuvadis_snp = readRDS(geuvadis_metadata$allGenotypesOutput[1])

	GTEx_lcl_snp$study = "GTEx (LCL)"
	geuvadis_snp$study = "Geuvadis (LCL)"

	GTEx_lcl_snp = GTEx_lcl_snp %>% select(-coverage_cut)
	geuvadis_snp = geuvadis_snp %>% select(-coverage_cut)

	GTEx_geuvadis_snp = rbind(GTEx_lcl_snp, geuvadis_snp)

	mean_M_S_plot = GTEx_geuvadis_snp[, .(mean_S = mean(S), mean_M = mean(M), nSamples = .N), by = .(chr, start, study)]
	mean_M_S_plot_GTEx_sample = mean_M_S_plot %>% filter(study == "GTEx (LCL)")
	mean_M_S_plot_GTEx_sample = mean_M_S_plot_GTEx_sample[sample(1:nrow(mean_M_S_plot_GTEx_sample), .01 * nrow(mean_M_S_plot_GTEx_sample)) ,]
	mean_M_S_plot_Geuvadis_sample = mean_M_S_plot %>% filter(study == "Geuvadis (LCL)")
	mean_M_S_plot_Geuvadis_sample = mean_M_S_plot_Geuvadis_sample[sample(1:nrow(mean_M_S_plot_Geuvadis_sample), .01 * nrow(mean_M_S_plot_Geuvadis_sample)) ,]
	mean_M_S_plot = rbind(mean_M_S_plot_GTEx_sample, mean_M_S_plot_Geuvadis_sample)

	M_S_plot_GTEx_sample = GTEx_geuvadis_snp %>% filter(study == "GTEx (LCL)")
	M_S_plot_GTEx_sample = M_S_plot_GTEx_sample[sample(1:nrow(M_S_plot_GTEx_sample), .001 * nrow(M_S_plot_GTEx_sample)) ,]
	M_S_plot_Geuvadis_sample = GTEx_geuvadis_snp %>% filter(study == "Geuvadis (LCL)")
	M_S_plot_Geuvadis_sample = M_S_plot_Geuvadis_sample[sample(1:nrow(M_S_plot_Geuvadis_sample), .001 * nrow(M_S_plot_Geuvadis_sample)) ,]
	M_S_plot = rbind(M_S_plot_GTEx_sample, M_S_plot_Geuvadis_sample)
	M_S_plot = M_S_plot %>% filter(true_genotype != "." & true_genotype != "./.")
	names(M_S_plot)[names(M_S_plot) == "true_genotype"] = "genotype"
	M_S_plot$genotype = gsub("1", "AA", M_S_plot$genotype)
	M_S_plot$genotype = gsub("2", "AB", M_S_plot$genotype)
	M_S_plot$genotype = gsub("3", "BB", M_S_plot$genotype)

	save(M_S_plot, mean_M_S_plot, file = plot_file)
}else {
	load(plot_file)
}


pdf("../figures/GTEx_LCL_vs_Geuvadis_M_S.pdf", width = 6, height = 3)

ggplot(mean_M_S_plot, aes(x = mean_S, y = mean_M)) + 
	geom_point(alpha = .05) + 
	facet_wrap(~study) +
	labs(x=expression(bar("S")), y=expression(bar("M"))) +
	scale_x_continuous(limits = c(1.1, 8))

ggplot(M_S_plot, aes(x = S, y = M, color = genotype)) + 
	geom_point(alpha = .05) + 
	facet_wrap(~study) +
	scale_x_continuous(limits = c(1.1, 8))

dev.off()

