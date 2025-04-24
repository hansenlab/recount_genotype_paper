library(matrixStats)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(data.table)
library(vroom)
theme_set(theme_cowplot())

####In this R session, I looked at SRA samples.
##I projected the SRA data onto 1000 Genome PCA that I have previously made in make_1kg_pca.R

pca<-readRDS("1000GenomePCA.rds")
pca_plot<-readRDS("pca_plot_1KGenome.rds")
pca_plot$zero_geno <- 0
pca_plot$sub_pop <- "1KG"
location<-readRDS("1000GenomePCA_location.rds")

#Read in the metadata:
metadata = read.csv("all_SRA.csv")
positions2<-location

for (i in 1:nrow(metadata)){
  print(i)
  cat(metadata$sample_id[i], "\n")
  geno<-as_tibble(vroom(metadata$genotypedSamples[i])) %>% select("chr","start","pred_genotype")
  if(nrow(geno)>1){
    geno$chr <- paste0(geno$chr,"_" ,geno$start)
    geno$start<- NULL
    #Make sure that we have the same number of dimentions in SRA as the original 1KG data:
    geno <- geno[geno$chr %in% positions2,]
    if(length(geno$chr)< length(positions2)){ 
      df<-data.frame(chr= positions2[!positions2 %in% geno$chr])
      geno<-full_join(geno, df)
    }
    #Make suer the order of the rows are the same as the original 1KG matrix:
    geno<-geno[match(positions2,geno$chr),]
    #replace the missing genotype with 0. Our genotype information is homo_ref=1,het=2,homo_alt=3.
    stopifnot(all.equal(geno$chr, positions2))
    #replace NA with 0 after making a matrix
    m<-as.matrix(geno[,-1])
    colnames(m)<-metadata$sample_id[i]
    m[is.na(m)]<- -1
    ##count number of 0s for each sample
    zero<-colSums(m== -1)
  
    # project new data onto the PCA space
    pca_projec <- scale(t(m), pca$center, pca$scale) %*% pca$rotation
    #Add to the plot
    pca_plot2<-data.frame(pc1=as.numeric(pca_projec[,1]), pc2=as.numeric(pca_projec[,2]), pop=as.factor(metadata$study[i]), sub_pop=as.factor(metadata$sample_id[i]), zero_geno=zero)
    pca_plot<-rbind(pca_plot,pca_plot2)
  }else{
  pca_plot2<-data.frame(pc1=NA, pc2=NA, pop=as.factor(metadata$study[i]), sub_pop=as.factor(metadata$sample_id[i]), zero_geno=NA)
  pca_plot<-rbind(pca_plot,pca_plot2)
  }}

fwrite(pca_plot,"sra_pca.csv.gz")
