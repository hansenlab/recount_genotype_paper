library(tidyverse)
library(data.table)

source("scripts/source.R")



GTEx_metadata = read.csv("GTEx_testing.csv") #Path to GTex testing metadata
major<-c()

figure<-"ready_to_plot/major_model.rds"
if(!file.exists(figure)){
  for(i in 1:nrow(GTEx_metadata)){
    print(i)
    GTEx_tissue = read.csv(GTEx_metadata$genotypedSamples[i])
    tissue<-GTEx_metadata$tissue[i]
    print(tissue)
    
    #Always predict the major allele: AF is the alternative allele frequency:
    GTEx_tissue<-GTEx_tissue %>% 
      mutate(pred_genotype=case_when(AF <= .5 ~ "1",
                                     AF > .5 ~ "3"))
    
    
    accuracy1<-conditional1(GTEx_tissue)
    accuracy2<-conditional2(GTEx_tissue)
    
    
    result = rbind(accuracy1, accuracy2)
    result = data.frame(tissue=tissue,
                        major_accuracy = sum(result$sum_major_accuracy)/sum(result$n_major),
                        minor_accuracy = sum(result$sum_minor_accuracy)/sum(result$n_minor),
                        overall_accuracy = sum(result$sum_overall_accuracy)/sum(result$n))
    
    print(result)
    
    major<-rbind(major,result)
  }
  
  saveRDS(major, figure)
}

major<-readRDS(figure)
plot<-data.frame(overall= round(mean(major$overall_accuracy),3),
                 major=round(mean(major$major_accuracy),3),
                 minor=round(mean(major$minor_accuracy),3))

plot<-plot %>% pivot_longer(everything(),names_to = "Type", values_to = "Accuracy")
plot$Type <- factor(plot$Type , levels=c("overall","major","minor"))


my_color_manual = c("#000",
                    "#6777cf",
                    "#D85E5E")
pdf(file="figure/major_model.pdf", width = 6, height = 3)
ggplot(plot, aes(x=Type,y=Accuracy, color= Type, fill=Type) )+ 
  geom_bar(stat="identity")+
  geom_text(aes(label=Accuracy), vjust = -0.3, color="black", size=3.5)+
  scale_fill_manual(values = my_color_manual) +
  scale_color_manual(values = my_color_manual)+
  theme_classic()
dev.off()


