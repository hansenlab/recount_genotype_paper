library(tidyverse)
library(data.table)
library(grid)
library(gridExtra)

theme_set(theme_cowplot())


metadata  = read.csv("/dcs04/hansen/data/recount_genotype/pipeline/AggregateFiles/GTEx_Blended_Tissue_Testing_metadata.csv")
metadata<-metadata[metadata$study=="Cells_EBV-transformed_lymphocytes",]
snp_ca = readRDS(metadata$SNP_conditionalAccuracy)
snp_ca = snp_ca %>% filter(chr != "chrX")

metadata_geu  = read.csv("/dcs04/hansen/data/recount_genotype/pipeline/AggregateFiles/Geuvadis_testing_metadata.csv")
snp_ca_geu = readRDS(metadata_geu$SNP_conditionalAccuracy[1])
snp_ca_geu = snp_ca_geu %>% filter(chr != "chrX")

cov_cut<-quantile(snp_ca$coverage_mean, probs = seq(0, 1, by = 1/6))
nSamp_cut<-quantile(snp_ca$nSamples, probs = seq(0, 1, by = 1/5))


snp_ca$coverage_cut = cut(snp_ca$coverage_mean, breaks=cov_cut, include.lowest = T)
snp_ca$nSamples_cut = cut(snp_ca$nSamples, breaks=unique(nSamp_cut), include.lowest = T)


snp_ca_geu$coverage_cut = cut(snp_ca_geu$coverage_mean, breaks=cov_cut, include.lowest = T)
snp_ca_geu$nSamples_cut = cut(snp_ca_geu$nSamples, breaks=quantile(snp_ca_geu$nSamples, probs = seq(0, 1, by = 1/5)), include.lowest = T) 

snp_ca_plot = snp_ca[!is.na(snp_ca$minor_accuracy)]
plot_gtex<-snp_ca_plot %>% group_by(coverage_cut, nSamples_cut) %>% mutate(mean_minor= mean(minor_accuracy)) %>% ungroup()


snp_ca_geu_plot = snp_ca_geu[!is.na(snp_ca_geu$minor_accuracy)]
plot_geu<-snp_ca_geu_plot %>% group_by(coverage_cut, nSamples_cut) %>% mutate(mean_minor= mean(minor_accuracy)) %>% ungroup()

theme1= theme(axis.text=element_text(size = 12), axis.title=element_text(size = 20))
pdf("~/recount_genotype/manuscript_figures/GTEx.5.pdf", width=8, height=5)
ggplot(plot_gtex, aes(nSamples_cut, mean_minor, group= coverage_cut)) + geom_line(aes(color=coverage_cut), linewidth=2)+
  theme1+ labs(x="Number of samples", y= "Mean minor allele accuracy")+
  ylim(0,1)
dev.off()

pdf("~/recount_genotype/manuscript_figures/GTEx.5.1.pdf", width=8, height=5)
ggplot(plot_geu, aes(nSamples_cut, mean_minor, group= coverage_cut)) + geom_line(aes(color=coverage_cut), linewidth=2)+
  theme1+ labs(x="Number of samples", y= "Mean minor allele accuracy") +
  ylim(0,1)
dev.off()
