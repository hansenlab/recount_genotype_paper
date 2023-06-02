library(tidyverse)
library(data.table)
library(cowplot)
theme_set(theme_cowplot())
plot_file = "../ready_to_plot/GTEx_LCL_vs_Geuvadis_accuracy.rds"

if(!file.exists(plot_file)) {
	GTEx_metadata = read.csv("/dcs04/hansen/data/recount_genotype/pipeline/AggregateFiles/GTEx_Blended_Tissue_Testing_metadata.csv")
	GTEx_LCL = readRDS(GTEx_metadata$conditionalAccuracy[GTEx_metadata$study == "Cells_EBV-transformed_lymphocytes"])
	GTEx_LCL = GTEx_LCL %>% pivot_longer(everything(), names_to = "type", values_to = "accuracy")
	GTEx_LCL$study = "GTEx LCL"

	Geuvadis_metadata = read.csv("/dcs04/hansen/data/recount_genotype/pipeline/AggregateFiles/Geuvadis_testing_metadata.csv")
	Geuvadis_LCL = readRDS(Geuvadis_metadata$conditionalAccuracy[1])
	Geuvadis_LCL = Geuvadis_LCL %>% pivot_longer(everything(), names_to = "type", values_to = "accuracy")
	Geuvadis_LCL$study = "Geuvadis LCL"

	GTEx_vs_Geuvadis = rbind(GTEx_LCL, Geuvadis_LCL)
	GTEx_vs_Geuvadis$type = gsub("major_accuracy", "major", GTEx_vs_Geuvadis$type)
	GTEx_vs_Geuvadis$type = gsub("minor_accuracy", "minor", GTEx_vs_Geuvadis$type)
	GTEx_vs_Geuvadis$type = gsub("overall_accuracy", "overall", GTEx_vs_Geuvadis$type)
	GTEx_vs_Geuvadis$type <- factor(GTEx_vs_Geuvadis$type, ordered = TRUE, levels = c("overall", "major", "minor"))
	GTEx_vs_Geuvadis$study <- factor(GTEx_vs_Geuvadis$study, ordered = TRUE, levels = c("GTEx LCL", "Geuvadis LCL"))
	saveRDS(GTEx_vs_Geuvadis, file = plot_file)
}else {
	GTEx_vs_Geuvadis = readRDS(plot_file)
}

my_color_manual = c("#000000", "#F8766D", "#00BFC4")

pdf("../figures/GTEx_LCL_vs_Geuvadis_accuracy.pdf", width = 6, height = 3)
ggplot(GTEx_vs_Geuvadis, aes(x = type, y = accuracy, fill = type)) + 
	geom_bar(stat = "identity") +
	geom_text(aes(label = round(accuracy, 3), y = round(accuracy, 3) + 0.05, color = type)) +
	labs(x = "", y = "Accuracy") +
	scale_fill_manual(values = my_color_manual) +
	scale_color_manual(values = my_color_manual) +
	theme(legend.position = "none") +
	theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = NA, fill = NA)) +
    facet_wrap(~study)
dev.off()

# #GTEx vs Geuvadis

# GTEx_lcl_snp_ca = readRDS(GTEx_metadata$SNP_conditionalAccuracy[GTEx_metadata$study == "Cells_EBV-transformed_lymphocytes"])
# GTEx_lcl_snp_ca$study = "GTEx LCL"
# GTEx_lcl_snp_ca$pred_major_alleles_correct = GTEx_lcl_snp_ca$major_accuracy * GTEx_lcl_snp_ca$true_major_alleles
# GTEx_lcl_snp_ca$pred_major_alleles_correct[is.na(GTEx_lcl_snp_ca$pred_major_alleles_correct)] = 0
# GTEx_lcl_snp_ca$pred_minor_alleles_correct = GTEx_lcl_snp_ca$minor_accuracy * GTEx_lcl_snp_ca$true_minor_alleles
# GTEx_lcl_snp_ca$pred_minor_alleles_correct[is.na(GTEx_lcl_snp_ca$pred_minor_alleles_correct)] = 0
# GTEx_lcl_snp_ca$overall_accuracy = (GTEx_lcl_snp_ca$pred_major_alleles_correct + GTEx_lcl_snp_ca$pred_minor_alleles_correct) / (GTEx_lcl_snp_ca$true_major_alleles + GTEx_lcl_snp_ca$true_minor_alleles)

# geuvadis_metadata  = read.csv("/dcs04/hansen/data/recount_genotype/pipeline/AggregateFiles/Geuvadis_testing_metadata.csv")
# geuvadis_snp_ca = readRDS(geuvadis_metadata$SNP_conditionalAccuracy[1])
# geuvadis_snp_ca$study = "Geuvadis LCL"
# geuvadis_snp_ca$pred_major_alleles_correct = geuvadis_snp_ca$major_accuracy * geuvadis_snp_ca$true_major_alleles
# geuvadis_snp_ca$pred_major_alleles_correct[is.na(geuvadis_snp_ca$pred_major_alleles_correct)] = 0
# geuvadis_snp_ca$pred_minor_alleles_correct = geuvadis_snp_ca$minor_accuracy * geuvadis_snp_ca$true_minor_alleles
# geuvadis_snp_ca$pred_minor_alleles_correct[is.na(geuvadis_snp_ca$pred_minor_alleles_correct)] = 0
# geuvadis_snp_ca$overall_accuracy = (geuvadis_snp_ca$pred_major_alleles_correct + geuvadis_snp_ca$pred_minor_alleles_correct) / (geuvadis_snp_ca$true_major_alleles + geuvadis_snp_ca$true_minor_alleles)


# GTEx_geuvadis_compare = rbind(GTEx_lcl_snp_ca, geuvadis_snp_ca)
# GTEx_geuvadis_compare = GTEx_geuvadis_compare %>% pivot_longer(cols = c("minor_accuracy", "major_accuracy", "overall_accuracy"), names_to = "type", values_to = "Accuracy")
# GTEx_geuvadis_compare$type[GTEx_geuvadis_compare$type == "major_accuracy"] = "major"
# GTEx_geuvadis_compare$type[GTEx_geuvadis_compare$type == "minor_accuracy"] = "minor"
# GTEx_geuvadis_compare$type[GTEx_geuvadis_compare$type == "overall_accuracy"] = "overall"
# GTEx_geuvadis_compare$type <- factor(GTEx_geuvadis_compare$type, ordered = TRUE, levels = c("overall", "major", "minor"))
# GTEx_geuvadis_compare$study <- factor(GTEx_geuvadis_compare$study, ordered = TRUE, levels = c("GTEx LCL", "Geuvadis LCL"))
# GTEx_geuvadis_compare_means <- aggregate(Accuracy ~  type * study, GTEx_geuvadis_compare, mean) 
# GTEx_geuvadis_compare_means$Accuracy = round(GTEx_geuvadis_compare_means$Accuracy, 2)

# my_color_manual = c("#000000", "#F8766D", "#00BFC4")

# pdf("../figures/GTEx_LCL_vs_Geuvadis_accuracy.pdf", width = 5, height = 3)
# ggplot(GTEx_geuvadis_compare, aes(type, Accuracy)) + 
# 	geom_boxplot(aes(color = type), outlier.shape = NA) +
# 	facet_wrap(~study) +
# 	stat_summary(fun=mean, geom="point", shape=20, size=5, aes(color = type, fill = type)) +
# 	geom_text(data = GTEx_geuvadis_compare_means, aes(label = Accuracy, y = Accuracy + 0.07, color = type)) +
# 	theme(legend.position = "none", axis.title.x=element_blank()) +
# 	scale_y_continuous(limits = c(0, 1)) + 
# 	scale_color_manual(values = my_color_manual)
# dev.off()
