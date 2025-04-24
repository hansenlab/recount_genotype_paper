library(tidyverse)
library(data.table)

setwd("~/recount_genotype/redo_manuscript_figures")
source("scripts/source.R")


#Calculate GTEX conditional accuracies:
GTEx_metadata = read.csv("/dcs07/hansen/data/recount_genotype/new_count_pipeline/new_count_pipeline/AggregateFiles/GTEx_testing.csv")
conditional_accuracies<-c()
figure<-"ready_to_plot/GTEx_testing_conditional_acc.rds"
if(!file.exists(figure)){
for(i in 1:nrow(GTEx_metadata)){
  print(i)
GTEx_tissue = read.csv(GTEx_metadata$genotypedSamples[i])
tissue<-GTEx_metadata$tissue[i]
print(tissue)

accuracy1<-conditional1(GTEx_tissue)
accuracy2<-conditional2(GTEx_tissue)


result = rbind(accuracy1, accuracy2)
result = data.frame(tissue=tissue,
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

plot2<-conditional_accuracies %>% pivot_longer(!tissue,names_to = "type", values_to = "Accuracy")
plot2$type[plot2$type=="major_accuracy"]<-"major"
plot2$type[plot2$type=="minor_accuracy"]<-"minor"
plot2$type[plot2$type=="overall_accuracy"]<-"overall"
plot2$type <- factor(plot2$type , levels=c("overall","major","minor"))


my_color_manual = c("#000",
                    "#6777cf",
                    "#D85E5E")
pdf(file="figure/GTEx_accuracy_all.pdf", width = 6, height = 3)
ggplot(plot, aes(x=type,y=Accuracy, color= type, fill=type) )+ 
  geom_bar(stat="identity")+
  geom_text(aes(label=Accuracy), vjust = -0.3, color="black", size=3.5)+
  scale_fill_manual(values = my_color_manual) +
  scale_color_manual(values = my_color_manual)+
  theme_classic()

dev.off()


pdf(file="figure/GTEx_accuracy.pdf", width = 6, height = 3)

ggplot(plot2, aes(x=type,y=Accuracy,color=type))+
geom_boxplot(outlier.shape=NA) +
  geom_jitter(position=position_jitter(height=0))+
  scale_fill_manual(values = my_color_manual) +
  scale_color_manual(values = my_color_manual)+
  theme_classic()
dev.off()

