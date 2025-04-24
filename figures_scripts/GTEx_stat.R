library(tidyverse)
library(data.table)



plot_file = "ready_to_plot/GTEx_stat.rds"
if(!file.exists(plot_file)) {
  GTEx_metadata = read.csv("GTEx_testing.csv") #Path to GTex testing metadata
  all_GTEx = list()
  for(i in 1:nrow(GTEx_metadata)) {
    print(i)
    GTEx_tissue = fread(GTEx_metadata$genotypedSamples[i])
    GTEx_tissue<-GTEx_tissue %>% group_by(sample_id_rep) %>%summarize(nSNP= n())
    GTEx_tissue$tissue = GTEx_metadata$tissue[i]
    all_GTEx[[i]] = GTEx_tissue
  }
  all_GTEx = do.call(rbind, all_GTEx)
  
  saveRDS(all_GTEx, plot_file)
}else {
  all_GTEx = readRDS(plot_file)
}



all_GTEx %>% summarize(mean=mean(nSNP),min=min(nSNP),max(nSNP))
all_GTEx %>% group_by(tissue) %>% summarize(nsample=n()) %>% ungroup() %>% summarize(mean=mean(nsample),min(nsample),max(nsample))
