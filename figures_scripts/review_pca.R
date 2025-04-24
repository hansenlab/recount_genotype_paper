library(matrixStats)
library(tidyverse)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
theme1= theme(axis.text=element_text(size = 10), axis.title=element_text(size = 15))
setwd("~/recount_genotype/redo_manuscript_figures")

plot_file = "ready_to_plot/review_pca.rds"

if(!file.exists(plot_file)) {
  
  geu_met<-read.delim("~/recount_genotype/study_metadata/geuvadis/input/EMBL_geuvadis_metadata.txt")
  
  sra<-read.csv("/dcs07/hansen/data/recount_genotype/new_count_pipeline/new_count_pipeline/AggregateFiles/all_SRA.csv")
  sra<-sra %>% filter(study=="ERP001942", sample_id %in% geu_met$Comment.ENA_RUN.)
  
  #-------------------------------------------
  #Do PCA
  #load in the PCs and the locations
  pca<-readRDS("/dcs04/hansen/data/recount_genotype/PCA/1kGenome_pIII/pca_plot/1000GenomePCA.rds")
  pca_plot<-readRDS("/dcs04/hansen/data/recount_genotype/PCA/1kGenome_pIII/pca_plot/pca_plot_1KGenome.rds")
  location<-readRDS("/dcs04/hansen/data/recount_genotype/PCA/1kGenome_pIII/pca_plot/1000GenomePCA_location.rds")
  
  #---------------------------------------------
  predict_genotype <- function(model_MS, model_MS_lattice, prior, M, S) {
    
    if(all(is.na(prior)) == T) {
      prior = c(.94, .04, .02) 
    }
    #Compute likelihood * prior for each genotype.
    likelihood_times_prior <- sapply(1:3, function(k){
      #Inference on mu_M and sd_M
      mu_M <- rep(NA, length(M))
      sd_M <- rep(NA, length(M))
      if(k == 1) {
        #model_MS_lattice is a data.table with key set to its S values, 
        #so we use the key to get lattice values.
        mu_M <- model_MS_lattice[J(S)]$mu_M_prediction_1
        sd_M <- model_MS_lattice[J(S)]$sd_M_prediction_1
        indx<-which(is.na(mu_M))
        mu_M[indx] <- predict(model_MS[[1]][[1]], data.frame(mu_S = S[indx]))
        sd_M[indx] <- predict(model_MS[[1]][[2]], data.frame(mu_S = S[indx]))
      }else if(k == 2) {
        mu_M <- model_MS_lattice[J(S)]$mu_M_prediction_2
        sd_M <- model_MS_lattice[J(S)]$sd_M_prediction_2
        indx<-which(is.na(mu_M))
        mu_M[indx] <- predict(model_MS[[2]][[1]], data.frame(mu_S = S[indx]))
        sd_M[indx] <- predict(model_MS[[2]][[2]], data.frame(mu_S = S[indx]))
      }else if(k == 3) {
        mu_M <- model_MS_lattice[J(S)]$mu_M_prediction_3
        sd_M <- model_MS_lattice[J(S)]$sd_M_prediction_3
        indx<-which(is.na(mu_M))
        mu_M[indx] <- predict(model_MS[[3]][[1]], data.frame(mu_S = S[indx]))
        sd_M[indx] <- predict(model_MS[[3]][[2]], data.frame(mu_S = S[indx]))
      }
      sd_M[sd_M <= 0] = 0 #if we predict (extrapolate) any variance to be <= 0, set it to 0. 
      return(prior[k] * dnorm(x = M, mean = mu_M, sd = sd_M, log = FALSE))
    })
    
    
    
    posterior <- sapply(1:3, function(k){
      if(is.null(dim(likelihood_times_prior))) { #if we only have one SNP that needs to be genotyped
        return(likelihood_times_prior[k] / sum(likelihood_times_prior))
      }else { #else, we have many SNPs in a dataframe
        return(likelihood_times_prior[, k] / rowSums(likelihood_times_prior))
      }
    })
    
    #Pick genotype with the max posterior genotype as our prediction.
    if(is.null(dim(posterior))) { #if we only have one SNP that needs to be genotyped
      predicted_genotype <- which(posterior == max(posterior))
    }else { #else, we have many SNPs in a dataframe
      predicted_genotype <- apply(posterior, MARGIN = 1, FUN = function(x) which(x == max(x)))
    }
    
    
    return(predicted_genotype)
  }
  modelLattice<-readRDS("/dcs07/hansen/data/recount_genotype/new_count_pipeline/new_count_pipeline/model/GTEx_Blended_Tissue_Training_genotyping_lattice.rds")
  prior<-readRDS("/dcs07/hansen/data/recount_genotype/new_count_pipeline/new_count_pipeline/model/GTEx_Blended_Tissue_Training_genotyping_prior.rds")
  model<-readRDS("/dcs07/hansen/data/recount_genotype/new_count_pipeline/new_count_pipeline/model/GTEx_Blended_Tissue_Training_genotyping_model.rds")
  
  
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
  
  
  
  
  i=1
  #----------
  #fixed geno
  geu_plot<-data.frame(pc1=numeric(),pc2=numeric(), pop=character())
  for(i in 62:nrow(sra)){
    print(i)
    run<-sra$sample_id[i]
    path<-paste0("/dcs07/hansen/data/recount_genotype/new_count_pipeline/new_count_pipeline/output/ERP001942/predict_genotype_accuracy/",
                 run, "_predGenotypes_final.csv.gz")
    one_samp<-read.csv(path)
    one_samp$sample_id_rep<-run
    
    total<-round(one_samp$coverage/2,0)
    id<-which(total<4)
    ref_count=round(one_samp$ref_count[id]/2,0)
    alt_count=total[id]-ref_count
    
    M <- log2((ref_count + 1) / (alt_count + 1))
    S <- log2(sqrt((ref_count + 1) * (alt_count + 1)))
    one_samp_try1<-one_samp[id,]
   
    one_samp_try1$pred_genotype <- predict_genotype(model_MS = model, 
                                          model_MS_lattice = modelLattice, 
                                          prior = prior,
                                          M = M, 
                                          S = S)
   
    
    geu_plot<-rbind(geu_plot,pca_df(one_samp_try1))
  }
  
  geuvadis_PCA<-rbind(pca_plot,geu_plot)
  
  geuvadis_PCA$Study<-"GEUVADIS"
  geuvadis_PCA$Study[1:2504]<-"1000Genome"
  
  saveRDS(geuvadis_PCA, "~/recount_genotype/redo_manuscript_figures/ready_to_plot/review_pca.rds")
}




varian<-as.numeric(summary(pca)$importance[2,1:10])
my_color<-c("#1290e3","#3246a8","tomato","steelblue1","#3279a8","#01548a")

pdf(file="~/recount_genotype/redo_manuscript_figures/figure/geuvadis_halfCov.pdf", width = 8, height = 4)
p2=ggplot(geuvadis_PCA[1:2504,])+
  geom_point(aes(pc1,pc2),alpha=0.4, color="grey")+
  scale_color_manual(values = my_color)+
  geom_point(data= geuvadis_PCA[-c(1:2504),], aes(pc1,pc2, color=pop),alpha=0.6)+coord_fixed(ratio = varian[2]/varian[1])+
  labs(x="PC1",y="PC2", color="Geuvadis RNA")+theme1
print(p2)
dev.off()

