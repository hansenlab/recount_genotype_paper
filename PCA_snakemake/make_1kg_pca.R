library(VariantAnnotation)


biallelic<-readRDS("pos_GTEx.rds")
st<-strsplit(biallelic, "_")
df=data.frame(chr=rep(NA,length(st) ), start=rep(0,length(st) ))
for(XX in 1:length(st)) {
  df$chr[XX]<-st[[XX]][1]
df$start[XX]<-as.numeric(st[[XX]][2])}

biallelic<-GRanges(seqnames= df$chr, IRanges(start=as.numeric(df$start), end=as.numeric(df$start)+1))

#load in parts of the VCF from 1K genome that we are interested in:
fl<-"1k_phase3.hg38_geno_proteinCoding.vcf.gz"
tab<-TabixFile(fl)
svp<-ScanVcfParam(geno="GT", which=biallelic)
vcf_rng<-readVcf(tab, "hg38", param= svp)
#get the genotype info:
geno<-geno(vcf_rng)$GT

#reformat the location:
loc=NA
chr<-strsplit(rownames(geno), "_")
for(XX in 1:length(chr)) {loc[XX]<-chr[[XX]][1]}
loc<-gsub(":","_",loc)
rownames(geno)<-loc
geno<-geno[!duplicated(loc),]
geno<-geno[order(rownames(geno)),]

#Make a matrix to be able to run the prcomp()
class(geno) <- "numeric"
id= which(is.na(rowSums(geno)))
geno_num<-geno[-id,]
location<-rownames(geno_num)
pca<-prcomp(t(geno_num))

pop<-read.csv("integrated_call_samples_v3.20130502.ALL.panel", sep = "\t")
population<- as.factor(sapply(1:ncol(geno), function(XX) { pop$super_pop[which(pop$sample == colnames(geno)[XX])[1]] }))
pca_plot$sub_pop<-as.factor(sapply(1:ncol(geno), function(XX) { pop$pop[which(pop$sample == colnames(geno)[XX])[1]] }))

pca_plot<-data.frame(pc1=as.numeric(pca$x[,1]), pc2=as.numeric(pca$x[,2]), pop=population)

#---------------------------------------------------
#Make plots
varian<-as.numeric(summary(pca)$importance[2,1:10])

pdf(file="~/plot/PCA1K.pdf", width = 12, height = 8)
p1=ggplot(pca_plot)+ geom_point(aes(pc1,pc2, color=pop))+
  labs(x=paste0("pc1(", varian[1]*100, "%)"),y=paste0("pc2(", varian[2]*100, "%)"), title=paste0("1000 Genome" ),subtitle= paste0("# of SNPs =",dim(geno_num)[1]))
print(p1)
dev.off()


pdf(file="~/plot/PCA1K_PC1_5.pdf", width = 12, height = 8)
p1=ggplot(pca_plot[1:2504,])+ geom_point(aes(pc1,pc2, color=sub_pop))+
  labs(x=paste0("pc1(", varian[1]*100, "%)"),y=paste0("pc2(", varian[2]*100, "%)"),
       title=paste0("1000 Genome" ),subtitle= paste0("# of SNPs =",dim(geno_num)[1]))
print(p1)

p1=ggplot(pca_plot[1:2504,])+ geom_point(aes(pca$x[,1],pca$x[,3], color=sub_pop))+
  labs(x=paste0("pc1(", varian[1]*100, "%)"),y=paste0("pc3(", varian[3]*100, "%)"), title=paste0("1000 Genome" ),subtitle= paste0("# of SNPs =",dim(geno_num)[1]))
print(p1)

p1=ggplot(pca_plot[1:2504,])+ geom_point(aes(pca$x[,1],pca$x[,4], color=sub_pop))+
  labs(x=paste0("pc1(", varian[1]*100, "%)"),y=paste0("pc4(", varian[4]*100, "%)"), title=paste0("1000 Genome" ),subtitle= paste0("# of SNPs =",dim(geno_num)[1]))
print(p1)

p1=ggplot(pca_plot[1:2504,])+ geom_point(aes(pca$x[,1],pca$x[,5], color=sub_pop))+
  labs(x=paste0("pc1(", varian[1]*100, "%)"),y=paste0("pc5(", varian[5]*100, "%)"), title=paste0("1000 Genome" ),subtitle= paste0("# of SNPs =",dim(geno_num)[1]))
print(p1)

p1=ggplot(pca_plot[1:2504,])+ geom_point(aes(pca$x[,2],pca$x[,3], color=sub_pop))+
  labs(x=paste0("pc2(", varian[2]*100, "%)"),y=paste0("pc3(", varian[3]*100, "%)"), title=paste0("1000 Genome" ),subtitle= paste0("# of SNPs =",dim(geno_num)[1]))
print(p1)


p1=ggplot(pca_plot[1:2504,])+ geom_point(aes(pca$x[,2],pca$x[,4], color=sub_pop))+
  labs(x=paste0("pc2(", varian[2]*100, "%)"),y=paste0("pc4(", varian[4]*100, "%)"), title=paste0("1000 Genome" ),subtitle= paste0("# of SNPs =",dim(geno_num)[1]))
print(p1)

p1=ggplot(pca_plot[1:2504,])+ geom_point(aes(pca$x[,2],pca$x[,5], color=sub_pop))+
  labs(x=paste0("pc2(", varian[2]*100, "%)"),y=paste0("pc5(", varian[5]*100, "%)"), title=paste0("1000 Genome" ),subtitle= paste0("# of SNPs =",dim(geno_num)[1]))
print(p1)

p1=ggplot(pca_plot[1:2504,])+ geom_point(aes(pca$x[,3],pca$x[,4], color=sub_pop))+
  labs(x=paste0("pc3(", varian[3]*100, "%)"),y=paste0("pc4(", varian[4]*100, "%)"), title=paste0("1000 Genome" ),subtitle= paste0("# of SNPs =",dim(geno_num)[1]))
print(p1)

p1=ggplot(pca_plot[1:2504,])+ geom_point(aes(pca$x[,3],pca$x[,5], color=sub_pop))+
  labs(x=paste0("pc3(", varian[3]*100, "%)"),y=paste0("pc5(", varian[5]*100, "%)"), title=paste0("1000 Genome" ),subtitle= paste0("# of SNPs =",dim(geno_num)[1]))
print(p1)

p1=ggplot(pca_plot[1:2504,])+ geom_point(aes(pca$x[,4],pca$x[,5], color=sub_pop))+
  labs(x=paste0("pc4(", varian[4]*100, "%)"),y=paste0("pc5(", varian[5]*100, "%)"), title=paste0("1000 Genome" ),subtitle= paste0("# of SNPs =",dim(geno_num)[1]))
print(p1)

dev.off()
