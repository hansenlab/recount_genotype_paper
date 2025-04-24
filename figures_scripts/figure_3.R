library(tidyverse)
library(data.table)
library(cowplot)
library(MetBrewer)


source("scripts/source.R")


cc_fill<- scale_fill_manual(values=c(met.brewer("Cassatt2",8)[4],met.brewer("Cassatt2",9)[2],met.brewer("Cassatt2",9)[7]))
cc_color <-scale_color_manual(values=c(met.brewer("Cassatt2",8)[4],met.brewer("Cassatt2",9)[2],met.brewer("Cassatt2",9)[7]))

make_plot<-function(conditional_accuracies){
  plot<-data.frame(overall= round(mean(conditional_accuracies$overall_accuracy),3),
                   major=round(mean(conditional_accuracies$major_accuracy),3),
                   minor=round(mean(conditional_accuracies$minor_accuracy),3))
  
  plot<-plot %>% pivot_longer(everything(),names_to = "Type", values_to = "Accuracy")
  plot$Type <- factor(plot$Type , levels=c("overall","major","minor"))
  return(plot)
}

major<-readRDS("ready_to_plot/major_model.rds")
gtex<-readRDS("ready_to_plot/GTEx_testing_conditional_acc.rds")
geuvadis<-readRDS("ready_to_plot/geuvadis_conditional_acc.rds")
gatk<-fread("ready_to_plot/gatkVSrecount.csv")


geu_plot<-make_plot(geuvadis)
geu_plot$study<-"Geuvadis Out-of-study"
gtex_plot<-make_plot(gtex)
gtex_plot$study<-"GTEx Test"

final_plot<-rbind(gtex_plot,geu_plot)
final_plot$study <- factor(final_plot$study , levels=c("GTEx Test","Geuvadis Out-of-study"))


p1<-ggplot(final_plot, aes(x=Type,y=Accuracy, color= Type, fill=Type) )+ 
  geom_bar(stat="identity")+
  geom_text(aes(label=Accuracy), vjust = -0.3, color="black", size=3.5)+
  cc_fill+
  cc_color+
  theme_classic()+
  theme(legend.position="none",strip.background = element_blank())+
  facet_wrap(vars(study))


major_plot<-make_plot(major)
p2<-ggplot(major_plot, aes(x=Type,y=Accuracy, color= Type, fill=Type) )+ 
  geom_bar(stat="identity")+
  geom_text(aes(label=Accuracy), vjust = -0.3, color="black", size=3.5)+
  cc_fill+
  cc_color+
  theme_classic()+
  theme(legend.position="none")



plot2<-gtex %>% pivot_longer(!tissue,names_to = "Type", values_to = "Accuracy")
plot2$Type[plot2$Type=="major_accuracy"]<-"major"
plot2$Type[plot2$Type=="minor_accuracy"]<-"minor"
plot2$Type[plot2$Type=="overall_accuracy"]<-"overall"
plot2$Type <- factor(plot2$Type , levels=c("overall","major","minor"))



p3<-ggplot(plot2, aes(x=Type,y=Accuracy,color=Type))+
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(position=position_jitter(height=0))+
  cc_fill+
  cc_color+
  theme_classic()+
  theme(legend.position="none")




gatk_plot<-gatk %>%
  group_by(pipeline) %>% 
  summarize(overall= round(mean(overall_accuracy),3),
            major=round(mean(major_accuracy),3),
            minor=round(mean(minor_accuracy),3)) %>%
  pivot_longer(!pipeline,names_to = "Type", values_to = "Accuracy")

gatk_plot$Type <- factor(gatk_plot$Type , levels=c("overall","major","minor"))
gatk_plot$pipeline <- factor(gatk_plot$pipeline , levels=c("Recount3","GATK"))

p4<-ggplot(gatk_plot, aes(y=Accuracy, x=Type, color=Type, fill=Type))+
  geom_bar(stat="identity")+
  geom_text(aes(label=Accuracy), vjust = -0.3, color="black", size=3.5)+
  cc_fill+
  cc_color+
  theme_classic()+
  theme(legend.position="none",strip.background = element_blank())+
  facet_wrap(vars(pipeline))

pdf(file="figure/figure_3.pdf", width = 7, height = 4.5)
plot_grid(p1,p2,p3,p4, labels = c('A', 'B', 'C', 'D'), label_size = 12, ncol = 2,rel_widths = c(1.25, 1))

dev.off()

#-------------------------------------------------------------------------------------

figure<-"ready_to_plot/GTEx_testing_conditional_acc.rds"
conditional_accuracies<-readRDS(figure)
gtex<-conditional_accuracies[conditional_accuracies$tissue=="Cells_EBV-transformed_lymphocytes",]



figure<-"ready_to_plot/geuvadis_conditional_acc.rds"
conditional_accuracies<-readRDS(figure)

plot<-data.frame(overall= round(mean(conditional_accuracies$overall_accuracy),3),
                 major=round(mean(conditional_accuracies$major_accuracy),3),
                 minor=round(mean(conditional_accuracies$minor_accuracy),3))
plot$tissue<-"Geuvadis LCL"
gtex$tissue<- "GTEx LCL"
colnames(gtex)<-c("tissue","major", "minor","overall")
plot<-rbind(plot,gtex)

plot<-plot %>% pivot_longer(!tissue,names_to = "Type", values_to = "Accuracy")
plot$Type <- factor(plot$Type , levels=c("overall","major","minor"))
plot$Accuracy<-round(plot$Accuracy,3)

pdf(file="figure/GTEx_vs_geuvadis.pdf", width = 5, height = 2.5)
ggplot(plot, aes(x=Type,y=Accuracy, color= Type, fill=Type) )+ 
  geom_bar(stat="identity")+
  geom_text(aes(label=Accuracy), vjust = -0.3, color="black", size=3.5)+
  cc_fill+
  cc_color+
  theme_classic()+
  theme(strip.background = element_blank())+
  facet_wrap(vars(tissue))

dev.off()




