library(ggplot2)
library(cowplot)
library(scattermore)
library(dplyr)
theme_set(theme_cowplot())
theme1= theme(axis.text=element_text(size = 10), axis.title=element_text(size = 20))
setwd("~/recount_genotype/redo_manuscript_figures")

my_color<-c("#DD77E4","#962638","#FFA81A","#585CD4","#5ba965","gray57")
sra_population<-readRDS("pca_plot_predPop.rds") #path to predicted population for SRA
p_no_pop<-fread("sra_pca.csv.gz") #path to PC loadings for SRA


#------------------------------------------------------------------------------------------------------------
#PCA plot for all sra samples
#------------------------------------------------------------------------------------------------------------


#Add missingness categories:
sra_population$missingness<-factor(sra_population$missingness, levels = c("0-10%", "10-30%", "30-70%",">70%"))
sra_population$missingness<-"0-10%"
sra_population$missingness[sra_population$percent_zero>10 & sra_population$percent_zero <=30]<-"10-30%"
sra_population$missingness[sra_population$percent_zero>30 & sra_population$percent_zero <=70]<-"30-70%"
sra_population$missingness[sra_population$percent_zero>70]<-">70%"


 ref<-p_no_pop[1:2505,]
 levels(ref$pop) <- c(levels(ref$pop), "SRA")
 ref$pop[2505]<-"SRA"
 
 #Make PCA plot for all SRA samples
 p0<- ggplot(sra_population[!is.na(sra_population$assay),], aes(x = pc1, y = pc2))+theme1+labs(x="PC1",y="PC2", color="Population")+scale_color_manual(values = my_color)+coord_fixed(ratio = 0.443609)+theme(plot.margin = margin(0, 0.5,0, 0.8))
 p1<- p0+geom_scattermore(data=sra_population[sra_population$missingness=="0-10%",], pointsize = 2.5, alpha =0.6, color="light grey",pixels=c(700,700))
 p1.1<-p1+geom_scattermore(data=ref,aes(x = pc1, y = pc2, color=pop), pointsize = 2.5,pixels=c(700,700))

 pdf(file="figure/PCA_allSRA.pdf", width = 6, height = 2.5)
 print(p1.1)
 dev.off()
 
 #------------------------------------------------------------------------------------------------------------
 #PCA plot2 
 #------------------------------------------------------------------------------------------------------------
 
 pca_plot<-readRDS("1KG_1percentSimulation.rds")
 
 pca_plot$percent_zero<-(pca_plot$zero_geno/30875)*100
 pca_plot$missingness<-"0%"
 pca_plot$missingness[pca_plot$percent_zero>16 & pca_plot$percent_zero <=20]<-"20%"
 pca_plot$missingness[pca_plot$percent_zero>45 & pca_plot$percent_zero <=50]<-"50%"
 pca_plot$missingness[pca_plot$percent_zero>75]<-"75%"
 pca_plot$missingness<-factor(pca_plot$missingness, levels = c("0%", "20%", "50%","75%"))
 
 dim(sra_population[sra_population$missingness=="30-70%"& sra_population$assay=="scrna-seq",])
 p0<- ggplot(sra_population[!is.na(sra_population$assay),], aes(x = pc1, y = pc2))+theme1+theme(plot.margin = margin(15, 0.5, 0, 0.8))+labs(x="PC1",y="PC2")+lims(x=c(min(sra_population$pc1),max(sra_population$pc1)),y=c(min(sra_population$pc2),max(sra_population$pc2)))
 p1<- p0+geom_scattermore(data=sra_population[sra_population$missingness=="0-10%"& sra_population$assay=="rna-seq",], pointsize = 3, alpha =0.2)+ theme(axis.title.x=element_blank(),axis.title.y =element_blank())
 p2<- p0+geom_scattermore(data=sra_population[sra_population$missingness=="10-30%"&sra_population$assay=="rna-seq",], pointsize = 3, alpha =0.2)+ theme(axis.title.x=element_blank(),axis.title.y =element_blank())
 p3<- p0+geom_scattermore(data=sra_population[sra_population$missingness=="30-70%"& sra_population$assay=="rna-seq",], pointsize = 3, alpha =0.2)+ theme(axis.title.x=element_blank(),axis.title.y =element_blank())
 p4<- p0+geom_scattermore(data=sra_population[sra_population$missingness==">70%"&sra_population$assay=="rna-seq",], pointsize = 3, alpha =0.2)+theme(axis.title.y =element_blank())
 p1.1<-p0+geom_scattermore(data=sra_population[sra_population$missingness=="0-10%"& sra_population$assay=="scrna-seq",], pointsize = 3, alpha =0.2)+ theme(axis.title.x=element_blank(), axis.title.y =element_blank())+scale_y_continuous( sec.axis = sec_axis(~ . + 10, name = "0-10%"),limits = c(min(sra_population$pc2),max(sra_population$pc2)))+ theme(axis.line.y.right=element_blank(), axis.text.y.right = element_blank(),axis.ticks.y.right =element_blank(),
                                                                                                                                                                                                                                                                                               axis.title.y.right = element_text( face="bold", size=20,angle =90))+theme(plot.margin = margin(15, 0.5, 5, 0.8))
 p2.1<- p0+geom_scattermore(data=sra_population[sra_population$missingness=="10-30%"&sra_population$assay=="scrna-seq",], pointsize = 3, alpha =0.2)+ theme(axis.title.x=element_blank(), axis.title.y =element_blank())+scale_y_continuous( sec.axis = sec_axis(~ . + 10, name = "10-30%"),limits = c(min(sra_population$pc2),max(sra_population$pc2)))+ theme(axis.line.y.right=element_blank(), axis.text.y.right = element_blank(),axis.ticks.y.right =element_blank(),
                                                                                                                                                                                                                                                                                                 axis.title.y.right = element_text( face="bold", size=20,angle =90))
 p3.1<- p0+geom_scattermore(data=sra_population[sra_population$missingness=="30-70%"& sra_population$assay=="scrna-seq",], pointsize = 3, alpha =0.2)+ theme(axis.title.x=element_blank(), axis.title.y =element_blank())+scale_y_continuous( sec.axis = sec_axis(~ . + 10, name = "30-70%"),limits = c(min(sra_population$pc2),max(sra_population$pc2)))+ theme(axis.line.y.right=element_blank(), axis.text.y.right = element_blank(),axis.ticks.y.right =element_blank(),
                                                                                                                                                                                                                                                                                                 axis.title.y.right = element_text( face="bold", size=20,angle =90))
 p4.1<- p0+geom_scattermore(data=sra_population[sra_population$missingness==">70%"&sra_population$assay=="scrna-seq",], pointsize = 3, alpha =0.2)+ theme(axis.title.y =element_blank())+scale_y_continuous( sec.axis = sec_axis(~ . +10, name = ">70%"),limits = c(min(sra_population$pc2),max(sra_population$pc2)))+ theme(axis.line.y.right=element_blank(), axis.text.y.right = element_blank(),axis.ticks.y.right =element_blank(),axis.title.y.right = element_text( face="bold", size=15,angle =90))
 
 p0<- ggplot(pca_plot, aes(x = pc1, y = pc2,color=pop))+theme1+theme(plot.margin = margin(15, 0.5, 0,0.8),legend.position = "none")+labs(x="PC1",y="PC2")+lims(x=c(min(sra_population$pc1),max(sra_population$pc1)),y=c(min(sra_population$pc2),max(sra_population$pc2)))+scale_color_manual(values = my_color)
 p1.2<- p0+geom_scattermore(data=pca_plot[pca_plot$percent_zero==0,], pointsize = 3, alpha =0.9)+ theme(axis.title.x=element_blank())#+scale_y_continuous( sec.axis = sec_axis(~ . + 10, name = "0%"),limits = c(min(p$pc2),max(p$pc2)))+ theme(axis.line.y.right=element_blank(), axis.text.y.right = element_blank(),axis.ticks.y.right =element_blank(),
                                                                                                                                                                                                                                               #  axis.title.y.right = element_text( face="bold", size=10,angle =90))
 p2.2<- p0+geom_scattermore(data=pca_plot[round(pca_plot$percent_zero,0) ==20,], pointsize = 3, alpha =0.9)+ theme(axis.title.x=element_blank())#+scale_y_continuous( sec.axis = sec_axis(~ . + 10, name = "20%"),limits = c(min(p$pc2),max(p$pc2)))+ theme(axis.line.y.right=element_blank(), axis.text.y.right = element_blank(),axis.ticks.y.right =element_blank(),
                                                                                                                                                                                                                                                  # axis.title.y.right = element_text( face="bold", size=10,angle =90))
 p3.2<- p0+geom_scattermore(data=pca_plot[round(pca_plot$percent_zero,0) ==50,], pointsize = 3, alpha =0.9)+ theme(axis.title.x=element_blank())#+scale_y_continuous( sec.axis = sec_axis(~ . + 10, name = "50%"),limits = c(min(p$pc2),max(p$pc2)))+ theme(axis.line.y.right=element_blank(), axis.text.y.right = element_blank(),axis.ticks.y.right =element_blank(),
                                                                                                                                                                                                                                                 #  axis.title.y.right = element_text( face="bold", size=10,angle =90))
 p4.2<- p0+geom_scattermore(data=pca_plot[round(pca_plot$percent_zero,0) ==75,], pointsize = 3, alpha =0.9)#+scale_y_continuous( sec.axis = sec_axis(~ . + 10, name = "75%"),limits = c(min(p$pc2),max(p$pc2)))+ theme(axis.line.y.right=element_blank(), axis.text.y.right = element_blank(),axis.ticks.y.right =element_blank(),
                                                                                                                                                                                                              #axis.title.y.right = element_text( face="bold", size=10,angle =90))
 m<-plot_grid(p1.2,p1,p1.1,p2.2, p2, p2.1,p3.2, p3,p3.1,p4.2, p4,p4.1,ncol=3, nrow=4, labels = c("1KG DNA","Bulk RNA-seq", "scRNA-seq"), label_size = 20)
 pdf(file="figure/PCA_final.pdf", width = 12, height = 9)
 print(m)
 dev.off()

 #------------------------------------------------------------------------------------------------------------
 #PCA plot3
 #------------------------------------------------------------------------------------------------------------
 
 p_prob<-readRDS("SRA_pop_acc.rds")
 sra_population$acc<-NA
 sra_population$acc<-p_prob$max_prob[match(sra_population$sub_pop, p_prob$sample_id)]
 pp2<-sra_population %>% filter(acc>0.7, !is.na(assay)) %>% group_by(assay,pred_pop) %>% summarise(n=n())%>% ungroup() %>% group_by(assay) %>%  mutate(percent=round(n/sum(n)*100,3))
 pp2$pred_pop<-factor(pp2$pred_pop, levels=c("SAS","EAS","AFR", "EUR","AMR"))
 pp2$assay[pp2$assay=="rna-seq"]<-"Bulk"
 pp2$assay[pp2$assay=="scrna-seq"]<-"scRNA-seq"
 sra_population%>% filter(assay=="scrna-seq") %>% summarise(n())

 
 my_color<-c("#5ba965","#FFA81A","#DD77E4","#585CD4","#962638","gray57")
 
 #my_color<-c("orchid2","springgreen3","tomato","steelblue1","yellow3","gray57")
 s<-ggplot(pp2, aes(x=assay,y=percent,fill=pred_pop))+ geom_bar(position="stack", stat="identity", alpha=0.9)+
   geom_text(aes(label =paste0(round(percent,1),"%")), position=position_stack(vjust = .5), fontface="bold")+
   scale_fill_manual(values = my_color)+
   labs(fill="Predicted\nPopulation",x="Assay type", y="Percentage")
 
pdf(file="figure/PCA_stackBar.pdf", width = 6, height = 6.5)
 print(s)
 dev.off()
 
 
 #------------------------------------------------------------------------------------------------------------
 #PCA plot4
 #------------------------------------------------------------------------------------------------------------
 p<-readRDS("pca_plot_predPop.rds")
 
 my_color<-c("#DD77E4","#962638","#FFA81A","#585CD4","#5ba965","gray57")
 plot_data_column = function (data, miss, type) {
   ggplot(data[data$missingness==miss,], aes(x = pc1, y = pc2))+
     geom_scattermore(data=data[which(data$missingness==miss & data$assay==type),],aes(x = pc1, y = pc2,color=pred_pop), pointsize = 2, alpha =0.4)+ 
     theme(plot.margin = margin(0, 0.5, 0, 0.8),legend.position = "none")+
     scale_color_manual(values = my_color)+
     labs(x="PC1",y="PC2")+lims(x=c(min(data$pc1),max(data$pc1)),y=c(min(data$pc2),max(data$pc2)))+
     annotate(geom = "text", label=paste0(miss,"%"), x=-38,y=10)
 }
 
 myplots <- lapply(levels(p$missingness), plot_data_column, data=p, type="rna-seq")
 for(i in 1:8){myplots[[i]]<- myplots[[i]]+theme(axis.title.x=element_blank())}
 for(i in c(2,4,6,8,10)){myplots[[i]]<- myplots[[i]]+theme(axis.title.y=element_blank())}
 mp <- plot_grid(plotlist = myplots,nrow = 5, ncol = 2, label_size = 13, align = "v", axis = "lr")
 
 
 pdf(file="figure/PCA_sup_bulk.pdf", width = 7, height = 9)
 print(mp)
 dev.off()
 
 myplots <- lapply(levels(p$missingness), plot_data_column, data=p, type="scrna-seq")
 for(i in 1:8){myplots[[i]]<- myplots[[i]]+theme(axis.title.x=element_blank())}
 for(i in c(2,4,6,8,10)){myplots[[i]]<- myplots[[i]]+theme(axis.title.y=element_blank())}
 mp <- plot_grid(plotlist = myplots,nrow = 5, ncol = 2, label_size = 13, align = "v", axis = "lr")
 
 
 pdf(file="figure/PCA_sup_scRNA.pdf", width = 7, height = 9)
 print(mp)
 dev.off()
 
 