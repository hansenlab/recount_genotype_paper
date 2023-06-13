library(matrixStats)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(data.table)
theme_set(theme_cowplot())
theme1= theme(axis.text=element_text(size = 10), axis.title=element_text(size = 15))
####In this R session, I looked at SRA samples.
##I projected the SRA data onto 1000 Genome PCA that I have previously made in make_1kg_pca.R

pca<-readRDS("/dcs04/hansen/data/recount_genotype/PCA/1kGenome_pIII/pca_plot/1000GenomePCA.rds")
p<-readRDS("/dcs04/hansen/data/recount_genotype/PCA/SRA/pca_plot.rds")
recount3_metadata<-fread("/dcs04/hansen/data/recount_genotype/PCA/SRA/Recount3_metadata.tsv", header= T, sep = "\t",quote="")
location<-readRDS("/dcs04/hansen/data/recount_genotype/PCA/1kGenome_pIII/pca_plot/1000GenomePCA_location.rds")
location<-30875
varian<-as.numeric(summary(pca)$importance[2,1:10])

p$percent_zero<-(p$zero_geno/length(location))*100
p$github_label<- recount3_metadata$seq_type[match(p$sub_pop,recount3_metadata$external_id)]
p$github_label[which(p$github_label== "smartseq")]<-"scRNA-seq"
p$github_label[which(p$github_label== "bulk")]<-"Bulk"

dim(p[p$percent_zero<10,])

p$missingness<-"10-30"
p$missingness[p$percent_zero>10 & p$percent_zero <=30]<-"10-30"
p$missingness[p$percent_zero>30 & p$percent_zero <=70]<-"30-70"
p$missingness[p$percent_zero>70]<-"70-100"

#Make the plots
pdf(file="~/recount_genotype/manuscript_figures/figures/PCA_sra.pdf", width = 8, height = 3.5)
ggplot(p[which(p$percent_zero <10),], aes(pc1,pc2))+ geom_point(alpha=0.4, color="dark grey", size= 0.7)+ geom_point(data=p[1:2504,], aes(pc1,pc2, color=pop),alpha=0.6, size=1)+
  coord_fixed(ratio = varian[2]/varian[1])+labs(x=paste0("PC1(", varian[1]*100, "%)"),y=paste0("PC2(", varian[2]*100, "%)"), color="Population")+theme1
dev.off()

q=p[-c(1:2504),]
q$pop<-NA
q<-rbind(q,p[1:2504,])

pdf(file="~/recount_genotype/manuscript_figures/figures/PCA_sra2.pdf", width = 9, height = 5)
ggplot(q[!is.na(q$missingness),])+ geom_point(aes(pc1,pc2,color=pop ),alpha=0.7, size= 0.6)+
  coord_fixed(ratio = varian[2]/varian[1])+ facet_grid(vars(missingness))+
  labs(x=paste0("PC1(", varian[1]*100, "%)"),y=paste0("PC2(", varian[2]*100, "%)"), color="Population")+theme1
dev.off()

