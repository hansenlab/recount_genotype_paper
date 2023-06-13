library(AnnotationHub)
library(data.table)
library(GenomicRanges)
library(glue)
library(Biostrings)
library(stringr)
library(rtracklayer)
Sys.setenv(ANNOTATION_HUB_CACHE=path.expand(rappdirs::user_cache_dir(appname="AnnotationHub")))

#get the gene locations that are protein coding only
ah <- AnnotationHub()

info <- query(ah, c("Homo.sapiens","Ensembl","GRCh38","gtf"))
grch38 <- info[['AH51014']]

#save(grch38, file = "/users/swang1/geuvadis_count/gene_anno_hg38.rda")

# subset for genes and chr1
chr_info <- c(1:22,"X","Y")
gene_info_hg38 <- grch38[grch38$type == 'gene']
gene_info_hg38_chr <- gene_info_hg38[seqnames(gene_info_hg38) %in% chr_info]
genetype_protein_coding <- gene_info_hg38_chr[gene_info_hg38_chr$gene_biotype %in% c("protein_coding","lincRNA","snoRNA")]
seqlevels(genetype_protein_coding)<-paste0("chr",seqlevels(genetype_protein_coding))
sum(width(reduce(genetype_protein_coding)))

#save(genetype_protein_coding,file="/users/arazi/extdata/genotype/gene_info_hg38_chr_gr.rda")
##ln -s /dcl01/hansen/data/arazi/gtex_count/ gtex_count
#all files will be in the shared folder on dcl01
#bcftools norm -m +any /dcl01/hansen/data/gtex_private/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.vcf.gz | bcftools view -m2 -M2 -v snps | bcftools view -s GTEX-1117F | bcftools query -f '%CHROM %POS %REF %ALT AF:[ %AF]\n' > /dcl01/hansen/data/arazi/gtex_count/genotype/biallelic_SNP2.txt
biallelic_SNP<-fread("/users/arazi/extdata/genotype/biallelic_SNP.txt")
biallelic_SNP_gr <- GRanges(seqnames=biallelic_SNP$V1,ranges = IRanges(biallelic_SNP$V2,biallelic_SNP$V2),ref_seq= biallelic_SNP$V3, alt_seq= biallelic_SNP$V4, allele_freq=biallelic_SNP$V6)

overlap<-findOverlaps(biallelic_SNP_gr, genetype_protein_coding)
biallelic_SNP_gr<-biallelic_SNP_gr[unique(queryHits(overlap))]

# find ref seq:
chrlist<-as.vector(seqnames(biallelic_SNP_gr))
poslist<-start(biallelic_SNP_gr)

## 1-based -- compared with UCSC Genome Browser
sequence_hg38 <- readDNAStringSet("/dcl02/lieber/ajaffe/recount-pump/refs/hg38.recount_pump.fa")


ref_chr_list <- lapply(c(1:22,"X"), function(x){
  as.character(sequence_hg38[[glue("chr{x}")]])
})

names(ref_chr_list) <- glue("chr{c(1:22,'X')}")

post_seq <- sapply(c(1:length(chrlist)), function(x){
  str_sub(ref_chr_list[[chrlist[x]]], start = poslist[x], end=poslist[x])
})
biallelic_SNP_gr$ref_seq<-as.character(post_seq)

alt_id <- which(sapply(biallelic_SNP_gr$alt_seq,nchar)>1)
ref_id <- which(sapply(biallelic_SNP_gr$ref_seq,nchar)>1)
rep_id <- unique(c(alt_id,ref_id))
if (length(rep_id) > 0) biallelic_SNP_gr <- biallelic_SNP_gr[-rep_id]
#saveRDS(post_seq,"/users/arazi/recount/recount3_ref.rds",compress=TRUE)
saveRDS(biallelic_SNP_gr, "/users/arazi/extdata/genotype/biallelic_SNP_gr.rda")

biallelic_SNP_bed <- data.frame(chromosome=seqnames(biallelic_SNP_gr),start=start(biallelic_SNP_gr),end=end(biallelic_SNP_gr)+1)

write.table(biallelic_SNP_bed,file="/users/arazi/extdata/genotype/biallelic_SNP_gr.bed",row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE)
#write.table(biallelic_SNP_bed,file="/users/arazi/extdata/genotype/biallelic_SNP_gr1.txt",row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE)
#bcftools view -R /users/arazi/extdata/genotype/biallelic_SNP_gr.bed /dcl01/hansen/data/gtex_private/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.vcf.gz | bcftools query -f '%CHROM %POS [ %GT]\n' > /users/arazi/gtex_count/genotype/real_total_geno.vcf
