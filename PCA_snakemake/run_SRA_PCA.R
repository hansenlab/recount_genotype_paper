library(matrixStats)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(data.table)
theme_set(theme_cowplot())

####In this R session, I looked at SRA samples.
##I projected the SRA data onto 1000 Genome PCA that I have previously made in make_1kg_pca.R

pca<-readRDS("/dcs04/hansen/data/recount_genotype/PCA/1kGenome_pIII/pca_plot/1000GenomePCA.rds")
pca_plot<-readRDS("/dcs04/hansen/data/recount_genotype/PCA/1kGenome_pIII/pca_plot/pca_plot_1KGenome.rds")
pca_plot$zero_geno = pca_plot$sub_pop = NA
location<-readRDS("/dcs04/hansen/data/recount_genotype/PCA/1kGenome_pIII/pca_plot/1000GenomePCA_location.rds")

#Read in the metadata:
metadata = read.csv("/dcs04/hansen/data/recount_genotype/pipeline/AggregateFiles/all_SRA.csv")
positions2<-location

for (i in 1:dim(metadata)[1]){
  print(i)
  cat(metadata$sample_id[i], "\n")
  geno<-as_tibble(read.csv(metadata$genotypedSamples[i])) %>% unique() %>% select("chr","start","pred_genotype")
  if(dim(geno)[1]>1){
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
  }}

saveRDS(pca_plot, file = "/dcs04/hansen/data/recount_genotype/PCA/SRA/pca_plot.rds")
