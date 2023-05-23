#Find the locations with the genotyoe information in all the GTEx tissues
#We only looked at tissues with more than 80 samples bc tissues with small sample sizes might have had bad library prep

metadata = read.csv("/dcs04/hansen/data/recount_genotype/pipeline/AggregateFiles/GTEx_Blended_Tissue_Testing_metadata.csv")
length(metadata$study)
for (i in 1:length(metadata$study)){
  cat(metadata$study[i], "\n")
  if(metadata$study[i] != "Whole_Blood" & metadata$study[i] != "Testis"){
  geno<-readRDS(metadata$allGenotypedSamplesAgg[i])
  
  geno$pred_genotype <- as.numeric(geno$pred_genotype)
  geno_wide<-geno %>% select(chr,start,pred_genotype, sample_id_rep) %>%
    pivot_wider(names_from = sample_id_rep, values_from = pred_genotype)#, values_fill =0)
  cat(ncol(geno_wide[-c(1:2)]), "\n")
  if(ncol(geno_wide[-c(1:2)]) > 80){
  geno_wide$chr <- paste0(geno_wide$chr,"_" ,geno_wide$start)
  geno_wide$start<-NULL
  
  geno_m<- as.matrix(geno_wide[,-1])
  id= which(is.na(rowSums(geno_m)))
  geno_m<-geno_m[-id,]
  loc<-geno_wide$chr[-id]
  
  #find intercept with geuvadis
  stopifnot(length(positions1)>0)
  positions1<-intersect(loc,positions1)
  print(length(positions1))
  n<-n+1
  }}}

saveRDS(positions1, file="/dcs04/hansen/data/recount_genotype/PCA/1kGenome_pIII/pca_plot/pos_GTEx.rds")



