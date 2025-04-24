library(tidyverse)
library(data.table)

setwd("~/recount_genotype/redo_manuscript_figures")
source("scripts/source.R")
my_color_manual = c("#000",
                    "#6777cf",
                    "#D85E5E")

#Calculate GTEX conditional accuracies:
#GTEx_metadata = read.csv("/dcs07/hansen/data/recount_genotype/new_count_pipeline/new_count_pipeline/AggregateFiles/GTEx_testing.csv")
conditional_accuracies<-c()
figure<-"ready_to_plot/geuvadis_conditional_acc.rds"
if(!file.exists(figure)){
  predict_gt_final<-fread("/dcs07/hansen/data/recount_genotype/new_count_pipeline/new_count_pipeline/geuvadis/aggregated/agg_geuvadis_w_true_gt.csv.gz")
  
  predict_gt_final<-predict_gt_final %>% mutate(true_genotype=case_when(
    true_genotype == "0/0" | true_genotype == "0|0" ~ 1,
    true_genotype == "0/1" | true_genotype == "1/0" | true_genotype == "1|0" | true_genotype == "0|1"~ 2,
    true_genotype == "1/1" | true_genotype == "1|1" ~ 3,
  ))
  run_id<-unique(predict_gt_final$sample_id_rep)
  for(i in 1:length(run_id)){
    print(i)
    sun_name=run_id[i]
    df<-predict_gt_final %>% filter(sample_id_rep== sun_name)
    
    accuracy1<-conditional1(df)
    accuracy2<-conditional2(df)
    
    
    result = rbind(accuracy1, accuracy2)
    result = data.frame(sample_id=sun_name,
                        major_accuracy = sum(result$sum_major_accuracy)/sum(result$n_major),
                        minor_accuracy = sum(result$sum_minor_accuracy)/sum(result$n_minor),
                        overall_accuracy = sum(result$sum_overall_accuracy)/sum(result$n))
    
    
    
    conditional_accuracies<-rbind(conditional_accuracies,result)
  }
  
  saveRDS(conditional_accuracies, figure)
}

conditional_accuracies<-readRDS(figure)

plot<-data.frame(overall= round(mean(conditional_accuracies$overall_accuracy),3),
                 major=round(mean(conditional_accuracies$major_accuracy),3),
                 minor=round(mean(conditional_accuracies$minor_accuracy),3))

plot<-plot %>% pivot_longer(everything(),names_to = "type", values_to = "Accuracy")
plot$type <- factor(plot$type , levels=c("overall","major","minor"))


pdf(file="figure/Geuvadis_accuracy.pdf", width = 6, height = 4)
ggplot(plot, aes(x=type,y=Accuracy, color= type, fill=type) )+ 
  geom_bar(stat="identity")+
  geom_text(aes(label=Accuracy), vjust = -0.3, color="black", size=3.5)+
  scale_fill_manual(values = my_color_manual) +
  scale_color_manual(values = my_color_manual)+
  theme_classic()

dev.off()

