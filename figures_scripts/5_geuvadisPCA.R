library(matrixStats)
library(tidyverse)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
theme1= theme(axis.text=element_text(size = 10), axis.title=element_text(size = 15))

pca<-readRDS("/dcs04/hansen/data/recount_genotype/PCA/1kGenome_pIII/pca_plot/1000GenomePCA.rds")

setwd("~/recount_genotype/redo_manuscript_figures/")
plot_file = "~/recount_genotype/redo_manuscript_figures/ready_to_plot/gauvadis_PCA.rds"

if(!file.exists(plot_file)) {
  
geu_met<-read.delim("~/recount_genotype/study_metadata/geuvadis/input/EMBL_geuvadis_metadata.txt")

sra<-read.csv("/dcs07/hansen/data/recount_genotype/new_count_pipeline/new_count_pipeline/AggregateFiles/all_SRA.csv")
sra<-sra %>% filter(study=="ERP001942", sample_id %in% geu_met$Comment.ENA_RUN.)

#-------------------------------------------
#Do PCA
#load in the PCs and the locations
pca_plot<-readRDS("/dcs04/hansen/data/recount_genotype/PCA/1kGenome_pIII/pca_plot/pca_plot_1KGenome.rds")
location<-readRDS("/dcs04/hansen/data/recount_genotype/PCA/1kGenome_pIII/pca_plot/1000GenomePCA_location.rds")
#Define a function to make the pca plot:

pca_df<-function(df){
  positions<-location
  
  geno_wide<-df %>% select(chr,start,pred_genotype, sample_id_rep) %>%
    pivot_wider(names_from = sample_id_rep, values_from = pred_genotype, values_fill =0)
  
  
  geno_wide$chr <- paste0(geno_wide$chr,"_" ,geno_wide$start)
  geno_wide$start<- NULL
  #Make sure that all the positions from 1KG are present in the Geuvadis geno data, if not just add them and make sure that they are in the same order as the 1KG data
  geno <- geno_wide[geno_wide$chr %in% positions,]
  if(length(geno$chr)< length(positions)){ 
    df<-data.frame(chr= positions[!positions %in% geno$chr])
    geno<-full_join(geno, df)
  }
  
  geno<-geno[match(positions,geno$chr),]
  stopifnot(all.equal(geno$chr, positions))
  geno[is.na(geno)]<-0
  m<-as.matrix(geno[,-1])
  
  # project new data onto the PCA space
  pca_projec <- scale(t(m), pca$center, pca$scale) %*% pca$rotation
  population2<- as.factor(sapply(1:ncol(geno[,-1]), function(XX) { geu_met$Characteristics.ancestry.category.[which(geu_met$Comment.ENA_RUN. == colnames(geno[,-1])[XX])[1]] }))
  
  pca_plot2<-data.frame(pc1=as.numeric(pca_projec[,1]), pc2=as.numeric(pca_projec[,2]), pop=population2)#, sample_id=sam)
  return(pca_plot2)
}





#----------
#fixed geno
geu_plot<-data.frame(pc1=numeric(),pc2=numeric(), pop=character())
for(i in 1:nrow(sra)){
  print(i)
  run<-sra$sample_id[i]
  path<-paste0("/dcs07/hansen/data/recount_genotype/new_count_pipeline/new_count_pipeline/output/ERP001942/predict_genotype_accuracy/",
               run, "_predGenotypes_final.csv.gz")
  one_samp<-read.csv(path)
  one_samp$sample_id_rep<-run
  geu_plot<-rbind(geu_plot,pca_df(one_samp))
}

geuvadis_PCA<-rbind(pca_plot,geu_plot)

geuvadis_PCA$Study<-"GEUVADIS"
geuvadis_PCA$Study[1:2504]<-"1000Genome"

saveRDS(geuvadis_PCA, "~/recount_genotype/redo_manuscript_figures/ready_to_plot/gauvadis_PCA.rds")
}

geuvadis_PCA<-readRDS(plot_file)
varian<-as.numeric(summary(pca)$importance[2,1:10])
pca_plot<-readRDS("/dcs04/hansen/data/recount_genotype/PCA/1kGenome_pIII/pca_plot/pca_plot_1KGenome.rds")

my_color<-c("#DD77E4","#962638","#FFA81A","#585CD4","#5ba965","gray57")
#my_color<-c("tomato","yellow3","springgreen3","steelblue1","orchid2","gray57")
p1=ggplot(pca_plot)+
  geom_point(aes(pc1,pc2, color=pop),alpha=0.6)+
  scale_color_manual(values = my_color)+
  coord_fixed(ratio = varian[2]/varian[1])+
  labs(x=paste0("PC1(", round(varian[1],3)*100, "%)"),y=paste0("PC2(", round(varian[2],3)*100, "%)"), color="1K Genome DNA")+theme1



my_color<-c("#1290e3","#3246a8","#ba11b5","steelblue1","#3279a8","#01548a")
p2=ggplot(geuvadis_PCA[1:2504,])+
  geom_point(aes(pc1,pc2),alpha=0.4, color="grey")+
  scale_color_manual(values = my_color)+
  geom_point(data= geuvadis_PCA[-c(1:2504),], aes(pc1,pc2, color=pop),alpha=0.6)+coord_fixed(ratio = varian[2]/varian[1])+
  labs(x="PC1",y="PC2", color="Geuvadis RNA")+theme1





figure5<-plot_grid(p1, p2,labels = "AUTO", ncol = 1)
pdf(file="~/recount_genotype/redo_manuscript_figures/figure/figure5.pdf", width = 7, height = 4.5)
print(figure5)
dev.off()


