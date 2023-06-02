library(optparse)

option_list <- list(
  make_option(c("-i", "--input"), type = "character",
              help = "Input: Biallelic Text file to be overlaped with the protein 
              coding region and outputed as a Grange"),
  make_option(c("-o", "--output"), type = "character",
              help = "Output: Granges of the biallelic SNPs that overlap the 
              protein coding region"))

opt <- parse_args(OptionParser(option_list = option_list))
if (is.na(opt$input) | is.na(opt$output)) {
  stop("Not all arguments provided. Check --help for description.")
}

library(data.table)
library(GenomicRanges)
library(glue)
library(Biostrings)
library(stringr)
library(rtracklayer)

#Read in the protein coding region that was saved from annotationHub "Ensembl","GRCh38". Appropriate code can be found in get_biallelic_snp.R
genetype_protein_coding <- readRDS("/users/arazi/extdata/genotype/gene_info_hg38_chr_gr.rda")

#Read in the biallelic SNP text file:
biallelic_SNP<-fread(opt$input)
#create Granges
biallelic_SNP_gr <- GRanges(seqnames=biallelic_SNP$V1,ranges = IRanges(biallelic_SNP$V2,biallelic_SNP$V2),ref_seq= biallelic_SNP$V3, alt_seq= biallelic_SNP$V4, allele_freq=biallelic_SNP$V6)

#find overlapwith the protein coding region
overlap<-findOverlaps(biallelic_SNP_gr, genetype_protein_coding)
biallelic_SNP_gr<-biallelic_SNP_gr[unique(queryHits(overlap))]

#find ref seq:
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

#double check if you have biallelic SNPS:
alt_id <- which(sapply(biallelic_SNP_gr$alt_seq,nchar)>1)
ref_id <- which(sapply(biallelic_SNP_gr$ref_seq,nchar)>1)
rep_id <- unique(c(alt_id,ref_id))
if (length(rep_id) > 0) biallelic_SNP_gr <- biallelic_SNP_gr[-rep_id]

saveRDS(biallelic_SNP_gr, opt$output)