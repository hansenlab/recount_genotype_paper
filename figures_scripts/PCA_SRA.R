library(matrixStats)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(data.table)
theme_set(theme_cowplot())
theme1= theme(axis.text=element_text(size = 10), axis.title=element_text(size = 15))




####In this R session, I looked at SRA samples.
##I projected the SRA data onto 1000 Genome PCA that I have previously made in make_1kg_pca.R
pca<-readRDS("1000GenomePCA.rds")
sra_pca<-fread("sra_pca.csv.gz")
location<-readRDS("1000GenomePCA_location.rds")
recount3_metadata<-fread("Recount3_metadata.tsv", header= T, sep = "\t",quote="")
#location<-30875
varian<-round(as.numeric(summary(pca)$importance[2,1:10]),3)

sra_pca$percent_zero<-(sra_pca$zero_geno/length(location))*100
sra_pca$github_label<- recount3_metadata$seq_type[match(sra_pca$sub_pop,recount3_metadata$external_id)]
sra_pca$github_label[which(sra_pca$github_label== "smartseq")]<-"scRNA-seq"
sra_pca$github_label[which(sra_pca$github_label== "bulk")]<-"Bulk"

dim(sra_pca[sra_pca$percent_zero<10,]) #[1] 65596     7

sra_pca$missingness<-"10-30"
sra_pca$missingness[sra_pca$percent_zero>10 & sra_pca$percent_zero <=30]<-"10-30"
sra_pca$missingness[sra_pca$percent_zero>30 & sra_pca$percent_zero <=70]<-"30-70"
sra_pca$missingness[sra_pca$percent_zero>70]<-"70-100"

#Make the plots
pdf(file="figure/PCA_sra.pdf", width = 8, height = 3.5)
ggplot(sra_pca[which(sra_pca$percent_zero <10),], aes(pc1,pc2))+ geom_point(alpha=0.4, color="dark grey", size= 0.7)+ geom_point(data=sra_pca[1:2504,], aes(pc1,pc2, color=pop),alpha=0.6, size=1)+
  coord_fixed(ratio = varian[2]/varian[1])+labs(x=paste0("PC1(", varian[1]*100, "%)"),y=paste0("PC2(", varian[2]*100, "%)"), color="Population")+theme1
dev.off()



