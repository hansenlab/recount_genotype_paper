library(ggplot2)
library(cowplot)
library(scattermore)
library(dplyr)
theme_set(theme_cowplot())
theme1= theme(axis.text=element_text(size = 20), axis.title=element_text(size = 30))


p<-readRDS("/Users/afroozrazi/Desktop/jhpce/plot/SRA_pop.rds")
p_no_pop<-readRDS("/Users/afroozrazi/Desktop/jhpce/plot/SRA.rds")
my_color<-c("tomato","yellow3","springgreen3","steelblue1","orchid2","gray57")

p$missingness<-"0-10%"
p$missingness[p$percent_zero>10 & p$percent_zero <=30]<-"10-30%"
p$missingness[p$percent_zero>30 & p$percent_zero <=70]<-"30-70%"
p$missingness[p$percent_zero>70]<-">70%"
p$missingness<-factor(p$missingness, levels = c("0-10%", "10-30%", "30-70%",">70%"))


 ref<-p_no_pop[1:2505,]
 levels(ref$pop) <- c(levels(ref$pop), "SRA")
 ref$pop[2505]<-"SRA"
 
 p0<- ggplot(p[!is.na(p$assay),], aes(x = pc1, y = pc2))+theme1+labs(x="PC1",y="PC2", color="Population")+scale_color_manual(values = my_color)+coord_fixed(ratio = 0.443609)+theme(plot.margin = margin(0, 0.5,0, 0.8))
 p1<- p0+geom_scattermore(data=p[p$missingness=="0-10%",], pointsize = 2.5, alpha =0.6, color="light grey",pixels=c(700,700))
 p1.1<-p1+geom_scattermore(data=ref,aes(x = pc1, y = pc2, color=pop), pointsize = 2.5, alpha =0.6,pixels=c(700,700))
 save_plot(filename = "/Users/afroozrazi/Desktop/plots/PCA_allSRA.pdf",p1.1,base_asp=2.254)
 
 pca_plot<-readRDS("/Users/afroozrazi/Desktop/jhpce/plot/PCA_SRA/1KG_1percentSimulation.rds")
 
 
 
 p0<- ggplot(p[!is.na(p$assay),], aes(x = pc1, y = pc2))+theme1+theme(plot.margin = margin(15, 0.5, 0, 0.8))+labs(x="PC1",y="PC2")+lims(x=c(min(p$pc1),max(p$pc1)),y=c(min(p$pc2),max(p$pc2)))
 p1<- p0+geom_scattermore(data=p[p$missingness=="0-10%"& p$assay=="rna-seq",], pointsize = 3, alpha =0.2)+ theme(axis.title.x=element_blank(),axis.title.y =element_blank())
 p2<- p0+geom_scattermore(data=p[p$missingness=="10-30%"&p$assay=="rna-seq",], pointsize = 3, alpha =0.2)+ theme(axis.title.x=element_blank(),axis.title.y =element_blank())
 p3<- p0+geom_scattermore(data=p[p$missingness=="30-70%"& p$assay=="rna-seq",], pointsize = 3, alpha =0.2)+ theme(axis.title.x=element_blank(),axis.title.y =element_blank())
 p4<- p0+geom_scattermore(data=p[p$missingness==">70%"&p$assay=="rna-seq",], pointsize = 3, alpha =0.2)+theme(axis.title.y =element_blank())
 p1.1<-p0+geom_scattermore(data=p[p$missingness=="0-10%"& p$assay=="scrna-seq",], pointsize = 3, alpha =0.2)+ theme(axis.title.x=element_blank(), axis.title.y =element_blank())+scale_y_continuous( sec.axis = sec_axis(~ . + 10, name = "0-10%"),limits = c(min(p$pc2),max(p$pc2)))+ theme(axis.line.y.right=element_blank(), axis.text.y.right = element_blank(),axis.ticks.y.right =element_blank(),
                                                                                                                                                                                                                                                                                               axis.title.y.right = element_text( face="bold", size=20,angle =90))+theme(plot.margin = margin(15, 0.5, 5, 0.8))
 p2.1<- p0+geom_scattermore(data=p[p$missingness=="10-30%"&p$assay=="scrna-seq",], pointsize = 3, alpha =0.2)+ theme(axis.title.x=element_blank(), axis.title.y =element_blank())+scale_y_continuous( sec.axis = sec_axis(~ . + 10, name = "10-30%"),limits = c(min(p$pc2),max(p$pc2)))+ theme(axis.line.y.right=element_blank(), axis.text.y.right = element_blank(),axis.ticks.y.right =element_blank(),
                                                                                                                                                                                                                                                                                                 axis.title.y.right = element_text( face="bold", size=20,angle =90))
 p3.1<- p0+geom_scattermore(data=p[p$missingness=="30-70%"& p$assay=="scrna-seq",], pointsize = 3, alpha =0.2)+ theme(axis.title.x=element_blank(), axis.title.y =element_blank())+scale_y_continuous( sec.axis = sec_axis(~ . + 10, name = "30-70%"),limits = c(min(p$pc2),max(p$pc2)))+ theme(axis.line.y.right=element_blank(), axis.text.y.right = element_blank(),axis.ticks.y.right =element_blank(),
                                                                                                                                                                                                                                                                                                 axis.title.y.right = element_text( face="bold", size=20,angle =90))
 p4.1<- p0+geom_scattermore(data=p[p$missingness==">70%"&p$assay=="scrna-seq",], pointsize = 3, alpha =0.2)+ theme(axis.title.y =element_blank())+scale_y_continuous( sec.axis = sec_axis(~ . +10, name = ">70%"),limits = c(min(p$pc2),max(p$pc2)))+ theme(axis.line.y.right=element_blank(), axis.text.y.right = element_blank(),axis.ticks.y.right =element_blank(),axis.title.y.right = element_text( face="bold", size=15,angle =90))
 
 p0<- ggplot(pca_plot, aes(x = pc1, y = pc2,color=pop))+theme1+theme(plot.margin = margin(15, 0.5, 0,0.8),legend.position = "none")+labs(x="PC1",y="PC2")+lims(x=c(min(p$pc1),max(p$pc1)),y=c(min(p$pc2),max(p$pc2)))
 p1.2<- p0+geom_scattermore(data=pca_plot[pca_plot$missingness=="0%",], pointsize = 3, alpha =0.9)+ theme(axis.title.x=element_blank())#+scale_y_continuous( sec.axis = sec_axis(~ . + 10, name = "0%"),limits = c(min(p$pc2),max(p$pc2)))+ theme(axis.line.y.right=element_blank(), axis.text.y.right = element_blank(),axis.ticks.y.right =element_blank(),
                                                                                                                                                                                                                                               #  axis.title.y.right = element_text( face="bold", size=10,angle =90))
 p2.2<- p0+geom_scattermore(data=pca_plot[pca_plot$missingness=="20%",], pointsize = 3, alpha =0.9)+ theme(axis.title.x=element_blank())#+scale_y_continuous( sec.axis = sec_axis(~ . + 10, name = "20%"),limits = c(min(p$pc2),max(p$pc2)))+ theme(axis.line.y.right=element_blank(), axis.text.y.right = element_blank(),axis.ticks.y.right =element_blank(),
                                                                                                                                                                                                                                                  # axis.title.y.right = element_text( face="bold", size=10,angle =90))
 p3.2<- p0+geom_scattermore(data=pca_plot[pca_plot$missingness=="50%",], pointsize = 3, alpha =0.9)+ theme(axis.title.x=element_blank())#+scale_y_continuous( sec.axis = sec_axis(~ . + 10, name = "50%"),limits = c(min(p$pc2),max(p$pc2)))+ theme(axis.line.y.right=element_blank(), axis.text.y.right = element_blank(),axis.ticks.y.right =element_blank(),
                                                                                                                                                                                                                                                 #  axis.title.y.right = element_text( face="bold", size=10,angle =90))
 p4.2<- p0+geom_scattermore(data=pca_plot[pca_plot$missingness=="75%",], pointsize = 3, alpha =0.9)#+scale_y_continuous( sec.axis = sec_axis(~ . + 10, name = "75%"),limits = c(min(p$pc2),max(p$pc2)))+ theme(axis.line.y.right=element_blank(), axis.text.y.right = element_blank(),axis.ticks.y.right =element_blank(),
                                                                                                                                                                                                              #axis.title.y.right = element_text( face="bold", size=10,angle =90))
 m<-plot_grid(p1.2,p1,p1.1,p2.2, p2, p2.1,p3.2, p3,p3.1,p4.2, p4,p4.1,ncol=3, nrow=4, labels = c("1KG DNA","Bulk RNA-seq", "scRNA-seq"), label_size = 20)
 save_plot(filename = "/Users/afroozrazi/Desktop/plots/PCA_final.pdf",m , nrow = 4, ncol = 3,base_asp=2.254)
 
 ggplot(p[!is.na(p$assay),], aes(x = pc1, y = pc2,color=pop))+geom_scattermore(pointsize = 2.5, alpha =0.2)+theme1+ facet_grid(missingness~assay)
 ##########
 
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
 save_plot(filename = "/Users/afroozrazi/Desktop/plots/PCA_s.pdf",mp , nrow = 5, ncol = 2,base_height = 2, base_width =4,base_asp=2.254)
 
 #####
 ###Pie plot:
 p_prob<-readRDS("/Users/afroozrazi/Desktop/jhpce/plot/SRA_pop_acc.rds")
 p$acc<-p_prob$max_prob
 #pp<-p[which(p$assay=="rna-seq" & p$acc>0.7),]
 #df_lab<-as.data.frame(table(pp$pred_pop))
 #colnames(df_lab)<-c("Population","N")
 #df_lab$prop<-paste0(round((df_lab$N/sum(df_lab$N))*100,1),"%")
 # my_color<-c("tomato","yellow3","springgreen3","steelblue1","orchid2","gray57")
 # 
 # 
 # s<- ggplot(df_lab, aes(x="", y=reorder(N, +N), fill=Population)) +
 #   geom_bar(stat="identity", width=1, alpha=0.7) +
 #   coord_polar("y", start=0) +
 #   geom_text(aes(label = paste0(Population, "\n",prop)), position = position_stack(vjust=0.5, reverse = FALSE)) +
 #   labs(x = NULL, y = NULL) +
 #   theme_classic() +
 #   theme(axis.line = element_blank(),
 #         axis.text = element_blank(),
 #         axis.ticks = element_blank(),
 #         legend.position = "none") +
 #   scale_fill_manual(values = my_color)
 # save_plot(filename = "/Users/afroozrazi/Desktop/plots/PCA_pie.pdf",s, base_height = 5, base_width =5)
 
 p$acc_cut<-cut(round(p$acc,1),10)
 pp2<-p %>% filter(acc>0.7, !is.na(assay)) %>% group_by(assay,pred_pop) %>% summarise(n=n())%>% ungroup() %>% group_by(assay) %>%  mutate(percent=round(n/sum(n)*100,3))
 pp2$pred_pop<-factor(pp2$pred_pop, levels=c("SAS","EAS","AFR", "EUR","AMR"))
 pp2$assay[pp2$assay=="rna-seq"]<-"Bulk"
 pp2$assay[pp2$assay=="scrna-seq"]<-"scRNA-seq"
 p%>% filter(assay=="scrna-seq") %>% summarise(n())

 
 my_color<-c("orchid2","springgreen3","tomato","steelblue1","yellow3","gray57")
 s<-ggplot(pp2, aes(x=assay,y=percent,fill=pred_pop))+ geom_bar(position="stack", stat="identity", alpha=0.9)+
   geom_text(aes(label =paste0(round(percent,1),"%")), position=position_stack(vjust = .5), fontface="bold")+
   scale_fill_manual(values = my_color)+
   labs(fill="Predicted\nPopulation",x="Assay type", y="Percentage")
 save_plot(filename = "/Users/afroozrazi/Desktop/plots/PCA_stackBar.pdf",s,base_height = 7, base_width =6.5)
 