#Code to make the figure 5 of the paper

library(tidyverse)
library(data.table)
library(grid)
library(gridExtra)
library(GenomicRanges)
library(cowplot)
theme_set(theme_cowplot())
theme1= theme(axis.text=element_text(size = 12), axis.title=element_text(size = 20), axis.text.x = element_text(angle = 30, vjust = 0.5, hjust=1))

metadata = read.csv("/dcs04/hansen/data/recount_genotype/pipeline/AggregateFiles/GTEx_Blended_Tissue_Testing_metadata.csv")

#Find the correctly predicted allele:
GTEx_lcl = readRDS(metadata$SNP_conditionalAccuracy[metadata$study == "Cells_EBV-transformed_lymphocytes"])
GTEx_lcl$study = "GTEx (LCL)"
GTEx_lcl$pred_major = GTEx_lcl$major_accuracy * GTEx_lcl$true_major_alleles
GTEx_lcl$pred_major[is.na(GTEx_lcl$pred_major)] = 0
GTEx_lcl$pred_minor = GTEx_lcl$minor_accuracy * GTEx_lcl$true_minor_alleles
GTEx_lcl$pred_minor[is.na(GTEx_lcl$pred_minor)] = 0
GTEx_lcl$overall_accuracy = (GTEx_lcl$pred_major + GTEx_lcl$pred_minor) / (GTEx_lcl$true_major_alleles + GTEx_lcl$true_minor_alleles)

geuvadis_metadata  = read.csv("/dcs04/hansen/data/recount_genotype/pipeline/AggregateFiles/Geuvadis_testing_metadata.csv")
geuvadis_lcl = readRDS(geuvadis_metadata$SNP_conditionalAccuracy[1])
geuvadis_lcl$study = "Geuvadis (LCL)"
geuvadis_lcl$pred_major = geuvadis_lcl$major_accuracy * geuvadis_lcl$true_major_alleles
geuvadis_lcl$pred_major[is.na(geuvadis_lcl$pred_major)] = 0
geuvadis_lcl$pred_minor = geuvadis_lcl$minor_accuracy * geuvadis_lcl$true_minor_alleles
geuvadis_lcl$pred_minor[is.na(geuvadis_lcl$pred_minor)] = 0
geuvadis_lcl$overall_accuracy = (geuvadis_lcl$pred_major + geuvadis_lcl$pred_minor) / (geuvadis_lcl$true_major_alleles + geuvadis_lcl$true_minor_alleles)



#3 GTEx minor accuracy as a function of coverage. 
i=20
GTEx_cut = quantile(GTEx_lcl$coverage_mean, probs = seq(0, 1, by = 1/i)) 


geuvadis_lcl$coverage_cut = cut(geuvadis_lcl$coverage_mean, breaks=unique(GTEx_cut), include.lowest = T)
GTEx_lcl$coverage_cut = cut(GTEx_lcl$coverage_mean, breaks=unique(GTEx_cut), include.lowest = T)

#geuvadis_lcl_plot = geuvadis_lcl[sample(1:nrow(geuvadis_lcl), .05 * nrow(geuvadis_lcl)) ,]
#GTEx_snp_ca_plot = GTEx_lcl[sample(1:nrow(GTEx_lcl), .05 * nrow(GTEx_lcl)) ,]

GTEx_geuvadis_compare2 = rbind(geuvadis_lcl, GTEx_lcl)
GTEx_geuvadis_compare2_means <- aggregate(minor_accuracy ~  coverage_cut * study, GTEx_geuvadis_compare2, mean) 
GTEx_geuvadis_compare2_means$minor_accuracy = round(GTEx_geuvadis_compare2_means$minor_accuracy, 2)


GTEx_geuvadis_compare2_means_join <- merge(aggregate(major_accuracy ~  coverage_cut * study, GTEx_geuvadis_compare2, mean),GTEx_geuvadis_compare2_means) 
GTEx_geuvadis_compare2_means_join$major_accuracy = round(GTEx_geuvadis_compare2_means_join$major_accuracy, 2)


compare_plot  = GTEx_geuvadis_compare2_means_join %>% pivot_longer(cols = c("major_accuracy", "minor_accuracy"), names_to = "type", values_to = "Accuracy")


pdf("~/recount_genotype/manuscript_figures/GTEx.3.2.pdf", width = 12, height = 5)
ggplot(compare_plot, aes(coverage_cut, Accuracy, group=type)) + geom_line(aes(color=type)) + geom_point()+ 
  #scale_x_discrete(labels=c("[4,4.33]" = "1", "(4.33,5.25]" = "2","(5.25,6.42]"="3", "(6.42,8.11]"="4", "(8.11,14.5]"= "5",  "(14.5,7.26e+04]" = "6"))+
  geom_text(data = compare_plot, color = "black", aes(label = Accuracy, y = Accuracy + 0.03)) +
  facet_wrap(~study)+ theme1+ labs(x= "Coverage bin")
dev.off()

SD<-aggregate(minor_accuracy ~  coverage_cut * study, GTEx_geuvadis_compare2, sd)
SD$minor_accuracy = round(SD$minor_accuracy, 2)
SD <- merge(aggregate(major_accuracy ~  coverage_cut * study, GTEx_geuvadis_compare2, sd),SD) 
SD$major_accuracy = round(SD$major_accuracy, 2)

SD_plot  = SD %>% pivot_longer(cols = c("major_accuracy", "minor_accuracy"), names_to = "type", values_to = "SD")


merged_plot<-merge(SD_plot, compare_plot, by=c("coverage_cut","study","type"))

head(merged_plot)
pdf("~/recount_genotype/manuscript_figures/GTEx.3.2.SD.pdf", width = 12, height = 5)
ggplot(SD_plot, aes(coverage_cut, SD, group=type)) + geom_line(aes(color=type)) + geom_point()+ 
  #scale_x_discrete(labels=c("[4,4.33]" = "1", "(4.33,5.25]" = "2","(5.25,6.42]"="3", "(6.42,8.11]"="4", "(8.11,14.5]"= "5",  "(14.5,7.26e+04]" = "6"))+
  geom_text(data = SD_plot, color = "black", aes(label = SD, y = SD + 0.03)) +
  facet_wrap(~study)+ theme1+ labs(x= "Coverage bin", y="SD", title="SD of major and minor accuracy for each coverage bin")
dev.off()


pdf("~/recount_genotype/manuscript_figures/GTEx.3.2.accSD.pdf", width = 12, height = 8)
ggplot(merged_plot, aes(coverage_cut, Accuracy, group=type)) + geom_line(aes(color=type)) + geom_point()+ 
  #scale_x_discrete(labels=c("[4,4.33]" = "1", "(4.33,5.25]" = "2","(5.25,6.42]"="3", "(6.42,8.11]"="4", "(8.11,14.5]"= "5",  "(14.5,7.26e+04]" = "6"))+
  geom_text(data = merged_plot, color = "black", aes(label = Accuracy, y = Accuracy + 0.03)) +
  facet_wrap(~study)+ theme1+ labs(x= "Coverage bin", y="Accuracy", title="Accuracy with SD for each coverage bin")+
  geom_errorbar(aes(x=coverage_cut, ymin=Accuracy-SD, ymax=1), width=0.25)
dev.off()

