library(tidyverse)
library(data.table)
library(cowplot)
library(gridExtra)
theme_set(theme_cowplot())

plot_file = "../ready_to_plot/GTEx_accuracy_overall.rds"

if(!file.exists(plot_file)) {
	metadata = read.csv("/dcs04/hansen/data/recount_genotype/pipeline/AggregateFiles/GTEx_Blended_Tissue_Testing_metadata.all_tissues.csv")
	model_result = readRDS(metadata$conditionalAccuracy[1])
	model_result$model = "GTEx Test"
	null_result = readRDS(metadata$majorAlleleModelConditionalAccuracy[1])
	null_result$model = "Major allele model"

	Geuvadis_metadata = read.csv("/dcs04/hansen/data/recount_genotype/pipeline/AggregateFiles/Geuvadis_testing_metadata.csv")
	Geuvadis_LCL = readRDS(Geuvadis_metadata$conditionalAccuracy[1])
	Geuvadis_LCL$model = "Geuvadis Out-of-study"

	result = rbind(model_result, rbind(null_result, Geuvadis_LCL))
	result = result %>% pivot_longer(!model, names_to = "type", values_to = "accuracy")

	result$type = gsub("major_accuracy", "major", result$type)
	result$type = gsub("minor_accuracy", "minor", result$type)
	result$type = gsub("overall_accuracy", "overall", result$type)
	result$type <- factor(result$type, ordered = TRUE, levels = c("overall", "major", "minor"))
	result$model <- factor(result$model, ordered = TRUE, levels = c("Major allele model", "GTEx Test", "Geuvadis Out-of-study"))
	saveRDS(result, file = plot_file)
} else {
	result = readRDS(plot_file)
}


plot_file = "../ready_to_plot/GTEX_accuracy_by_tissue.rds"
if(!file.exists(plot_file)) {

	metadata = read.csv("/dcs04/hansen/data/recount_genotype/pipeline/AggregateFiles/GTEx_Blended_Tissue_Testing_metadata.csv")

	major_minor_acc = list()
	for(i in 1:nrow(metadata)) {
		cat(metadata$study[i], "\n")
		acc = readRDS(metadata$conditionalAccuracy[i])
		acc$study = metadata$study[i]
		major_minor_acc[[i]] = as.data.frame(acc)
	}
	major_minor_acc = do.call(rbind, major_minor_acc)

	major_minor_acc_long  = major_minor_acc %>% 
		select(major_accuracy, minor_accuracy, overall_accuracy, study) %>%
		pivot_longer(cols = c("overall_accuracy", "major_accuracy", "minor_accuracy"), names_to = "type", values_to = "Accuracy")

	major_minor_acc_long$type = gsub("major_accuracy", "major", major_minor_acc_long$type)
	major_minor_acc_long$type = gsub("minor_accuracy", "minor", major_minor_acc_long$type)
	major_minor_acc_long$type = gsub("overall_accuracy", "overall", major_minor_acc_long$type)
	major_minor_acc_long$type <- factor(major_minor_acc_long$type, ordered = TRUE, levels = c("overall", "major", "minor"))
	major_minor_acc_long$study_short = substr(major_minor_acc_long$study, 1, 14)

	saveRDS(major_minor_acc_long, file = plot_file)

}else {
	major_minor_acc_long = readRDS(plot_file)
}



my_color_manual = c("#000000", "#F8766D", "#00BFC4")

pdf("../figures/GTEx_accuracy_overall.pdf", width = 8, height = 4.65)

p1_data = result %>% filter(model == "Major allele model")
p1_data$facet = "GTEx Test"
p1 = ggplot(p1_data, aes(x = type, y = accuracy, fill = type)) + 
	geom_bar(stat = "identity") +
	geom_text(aes(label = round(accuracy, 3), y = round(accuracy, 3) + 0.05, color = type)) +
	facet_wrap(~facet) +
	labs(x = "", y = "Accuracy") +
	scale_fill_manual(values = my_color_manual) +
	scale_color_manual(values = my_color_manual) +
	ggtitle("Major Allele Model") +
	#theme(plot.title = element_text(margin=margin(0,0,28,0)))+
	theme(legend.position = "none") +
	theme(plot.title = element_text(hjust = 0.5)) +
	theme(plot.title = element_text(size=14)) +
	theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = NA, fill = NA)) 

p2 = ggplot(result %>% filter(model != "Major allele model"), aes(x = type, y = accuracy, fill = type)) + 
	geom_bar(stat = "identity") +
	geom_text(aes(label = round(accuracy, 3), y = round(accuracy, 3) + 0.05, color = type)) +
	facet_wrap(~model) +
	labs(x = "", y = "Accuracy") +
	scale_fill_manual(values = my_color_manual) +
	scale_color_manual(values = my_color_manual) +
	theme(legend.position = "none") +
	theme(plot.title = element_text(hjust = 0.5)) +
	theme(plot.title = element_text(size=14)) +
	ggtitle("Our Model") +
		theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = NA, fill = NA)) 

major_minor_acc_long$facet = "GTEx Test by tissue"
p3 = ggplot(major_minor_acc_long, aes(x = type, y = Accuracy, color = type)) + 
	geom_boxplot(outlier.shape=NA) +
	geom_jitter(position=position_jitter(height=0)) +
	facet_wrap(~facet) +
	scale_y_continuous(limits = c(.85, 1)) +
	scale_color_manual(values = my_color_manual) +
	theme(legend.position = "none") +
	theme(axis.title.x=element_blank()) +
	theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = NA, fill = NA)) 

lay <- rbind(c(1,2,2),
             c(NA,3,3))

grid.arrange(p1, p2, p3, layout_matrix = lay, heights = 3:2)

dev.off()
