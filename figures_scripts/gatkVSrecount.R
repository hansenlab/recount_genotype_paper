library(vroom)
library(data.table)
library(tidyverse)
library(MetBrewer)

cc_fill <-scale_fill_manual(values=met.brewer("Monet"))
cc_color <-scale_color_manual(values=met.brewer("Monet"))
source("~/recount_genotype/redo_manuscript_figures/scripts/source.R")
setwd("~/recount_genotype/redo_manuscript_figures")

gatk_metadata<-read.csv("/dcs07/hansen/data/recount_genotype/new_count_pipeline/new_count_pipeline/AggregateFiles/gatk_geuvadis_gt.csv")

#geu_met<-read.delim("recount_genotype/study_metadata/geuvadis/input/EMBL_geuvadis_metadata.txt")
geu_met<-read.csv("/users/arazi/recount_genotype/GATK_snakemake/test.csv")
opt<-c()
opt$allGenotypedSamplesAgg<-"/dcs07/hansen/data/recount_genotype/new_count_pipeline/new_count_pipeline/geuvadis/aggregated/agg_geuvadis_w_true_gt.csv.gz"

geu<-fread(opt$allGenotypedSamplesAgg)
geu<-geu %>% 
  mutate(true_genotype=case_when(
    true_genotype == "0/0" | true_genotype == "0|0" ~ 1,
    true_genotype == "0/1" | true_genotype == "1/0" | true_genotype == "1|0" | true_genotype == "0|1"~ 2,
    true_genotype == "1/1" | true_genotype == "1|1" ~ 3))

plotFile<-"~/recount_genotype/redo_manuscript_figures/ready_to_plot/gatkVSrecount.csv"
plotFile2<-"~/recount_genotype/redo_manuscript_figures/ready_to_plot/gatkVSrecount_stat.csv"

if(!file.exists(c(plotFile,plotFile2))){
gatkVSrecount<-c()
stat_df<-c()
samp_unique<-unique(geu$sample_id_rep)
for(i in 1:nrow(geu_met)){
  print(i)
  samp_id<-geu_met$sample_id[i]
  
  if(sum(samp_unique %in% samp_id)>0){
  sample_wide<-geu %>% 
    filter(sample_id_rep %in% samp_id) %>% 
    select(chr,start,AF,pred_genotype,true_genotype,sample_id_rep) %>% 
    pivot_wider(names_from = sample_id_rep, values_from = c(pred_genotype,true_genotype))
  colnames(sample_wide)<-c("chr", "start","AF", "pred_genotype", "true_genotype")
  
  gatk_geno<-fread(gatk_metadata$gatk_gt_format[gatk_metadata$run_id==samp_id])
  colnames(gatk_geno)<-c("chr", "start", "gatk_geno")
  
  gatk_geno<-gatk_geno%>% mutate(gatk_geno=case_when(
    gatk_geno == "0/0" | gatk_geno == "0|0" ~ 1,
    gatk_geno == "0/1" | gatk_geno == "1/0" | gatk_geno == "1|0" | gatk_geno == "0|1"~ 2,
    gatk_geno == "1/1" | gatk_geno == "1|1" ~ 3)) 
  
  joined<-left_join(sample_wide,gatk_geno)
  
  #------------------
  #Get some stats
  gatk_filter<-gatk_geno %>% filter(!is.na(gatk_geno), chr %in% paste0("chr", 1:22))
  
  stat = data.frame(sample_id=samp_id,
                    recount_nrow=nrow(joined),
                    gatk_nrow= nrow(gatk_filter),
                    have_true= sum(!is.na(joined$true_genotype)),
                    comparison_nrow= sum(!is.na(joined$true_genotype) &  !is.na(joined$gatk_geno)))
  stat_df<-rbind(stat_df,stat)

  #-----------------------------


  joined<-joined %>%  filter(!is.na(true_genotype), !is.na(gatk_geno))

  joined_gatk<- joined %>% select(-pred_genotype)
  colnames(joined_gatk)[which(colnames(joined_gatk)=="gatk_geno")]<-"pred_genotype"


  #----------------------
  #calculate conditional accuracies for Recount3:
  #----------------------
  accuracy1<-conditional1(joined)
  accuracy2<-conditional2(joined)
  result = rbind(accuracy1, accuracy2)
  result_recount3 = data.frame(pipeline="Recount3",
                               sample_id=samp_id,
                      major_accuracy = sum(result$sum_major_accuracy)/sum(result$n_major),
                      minor_accuracy = sum(result$sum_minor_accuracy)/sum(result$n_minor),
                      overall_accuracy = sum(result$sum_overall_accuracy)/sum(result$n))


  #----------------------
  #calculate conditional accuracies for GATK:
  #----------------------
  df<-joined_gatk
  accuracy1<-conditional1(df)
  accuracy2<-conditional2(df)
  result = rbind(accuracy1, accuracy2)
  result_gatk = data.frame(pipeline="GATK",
                           sample_id=samp_id,
                      major_accuracy = sum(result$sum_major_accuracy)/sum(result$n_major),
                      minor_accuracy = sum(result$sum_minor_accuracy)/sum(result$n_minor),
                      overall_accuracy = sum(result$sum_overall_accuracy)/sum(result$n))
  result<-rbind(result_recount3,result_gatk)
  gatkVSrecount<-rbind(gatkVSrecount,result)


  fwrite(gatkVSrecount,plotFile)
  fwrite(stat_df, plotFile2)
  
}
}
}else{
  gatkVSrecount<-fread(plotFile)
  stat_df<-fread(plotFile2)
}

plot_df<-gatkVSrecount %>% group_by(pipeline) %>% summarize(overall= round(mean(overall_accuracy),3),
                                                major=round(mean(major_accuracy),3),
                                                minor=round(mean(minor_accuracy),3))

plot_df<-plot_df %>% pivot_longer(!pipeline,names_to = "type", values_to = "Accuracy")
plot_df$type <- factor(plot_df$type , levels=c("overall","major","minor"))

plot_df$pipeline <- factor(plot_df$pipeline , levels=c("Recount3","GATK"))


pdf(file="figure/gatkVSrecount.pdf", width = 10, height = 6)
ggplot(plot_df, aes(y=Accuracy, x=type, color=type, fill=type))+
  geom_bar(stat="identity")+
  geom_text(aes(label=Accuracy), vjust = -0.3, color="black", size=3.5)+
  cc_fill+
  cc_color+
  theme_classic()+
  labs(title="GATK vs Recount3")+
  facet_wrap(vars(pipeline))
dev.off()
