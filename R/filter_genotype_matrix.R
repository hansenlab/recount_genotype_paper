library(optparse)

option_list <- list(
    make_option(c("-m", "--metadata"), type = "character",
    help = "Input: CSV file of Metadata. Contains columns: sample_id, total,
    alt, study. Each row is a sample."),
    make_option(c("-t", "--study"), type = "character",
    help = "Input: study name to run analysis on. Must be a study name
    found in Metadata specified."),
    make_option(c("-f", "--filteredSNPs"), type = "character",
    help = "Input: R object (GRanges): SNPs that intersected with BigWigs."),
    make_option(c("-g", "--genotypes"), type = "character",
    help = "Input: CSV with rows: genotypes that intersected with
    filteredSNPs; columns: all samples. "),
    make_option(c("--filteredGenotypes"), type = "character",
    help = "Output: R object dataframe: genotype matrix that matches dimension
    of alt and ref matricies: rows: filtered SNPs. columns: POS, REF, 
    followed by sample names."))

opt <- parse_args(OptionParser(option_list = option_list))
if (is.na(opt$filteredSNPs) | is.na(opt$genotypes) |
    is.na(opt$filteredGenotypes) | is.na(opt$metadata) |
    is.na(opt$study)) {
    stop("Not all arguments provided. Check --help for description.")
}

library(data.table)
library(rtracklayer)
library(GenomicRanges)
options(scipen=999)

#load everything in
metadata <- read.csv(opt$metadata, stringsAsFactors = FALSE)
metadata <- metadata[metadata$study == opt$study ,]
if(endsWith(opt$filteredSNPs, ".rds")) {
    snps_gr <- readRDS(opt$filteredSNPs)
}else if(endsWith(opt$filteredSNPs, ".txt")) {
     snps_gr <- fread(opt$filteredSNPs, header = F)
     snps_gr <- snps_gr[snps_gr$V1 != "" ,] #weird trailing whitespace issue from python pandas.
}else {
    stop("Bad Filter SNPs provided: needs to be txt or .rds (GRanges) file.")
}
genotypes <- fread(opt$genotypes)

cat("Loaded in files. Processing. \n")
print(gc())

#Genotype rows: 
#We *shouldn't* have to subset our genomic coordinates because we got the genotypes
#from the coordinates of the filtered SNPs in `get_genotype_matrix` step. However, 
#in that step, if there are deletions that overlap the site of interest, they will be 
#pulled out of the VCF in the genotypes. 
#So we actually get a small amount of deletions that are not of our interest that 
#need to be filtered out.
#For GTEx analysis, we would expect that all SNPs are found in genotype file as 
#our list of SNPs came out of filtering from GTEx genotype file (VCF).
#For Geuvadis analysis, only some of the SNPs are found in the genotype file as
#our list of SNPs came out of filtering from GTEx genotype file, which does not
#have to match the Geuvadis VCF file!
genotypes_chr_pos <- paste0(genotypes$CHROM, "_", genotypes$POS)
if(endsWith(opt$filteredSNPs, ".rds")) {
    snps_gr_pos <- paste0(as.character(seqnames(snps_gr)), "_", 
        start(snps_gr))
}else if(endsWith(opt$filteredSNPs, ".txt")) {
    snps_gr_pos <- paste0(snps_gr$V1, "_", snps_gr$V2)
}
match_id <- match(snps_gr_pos, genotypes_chr_pos)
if(sum(is.na(match_id)) > 0) { #segfault?
    warning("Some sites from list of SNPs not found in genotype file:\n", 
        paste(head(snps_gr_pos[is.na(match_id)]), collapse = "\n"))
}
#todo: this step is memory consuming: how to subset rows without doubling memory usage?
snps_gr_idx <- which(!is.na(match_id))
genotypes_idx <- match_id[which(!is.na(match_id))]
snps_gr <- snps_gr[snps_gr_idx ,]
genotypes <- genotypes[genotypes_idx ,]
print(gc())

#Genotype columns:
#Subset genotypes to make sure our samples match genotyped individuals:
#- we only need a subset of individuals for downstream analysis
#- some individuals may need to be repeated, 
#as it is possible to have multiple samples from same individual. 

match_id <- match(metadata$individual_id, colnames(genotypes))   
if(sum(is.na(match_id)) > 0) {
    warning("Some individuals from metadata not found in genotype file.\n
             Will proceed but have NA for those columns:\n", 
        paste(metadata$individual_id[is.na(match_id)], collapse = "\n"))
}
metadata_individals_in_genotype_idx <- which(!is.na(match_id))
metadata <- metadata[metadata_individals_in_genotype_idx ,]
match_id <- match_id[which(!is.na(match_id))]
genotypes_final <- genotypes[, ..match_id]
colnames(genotypes_final) <- metadata$sample_id_rep
#add back CHROM and POS. 
genotypes_final$CHROM <- genotypes$CHROM
genotypes_final$POS <- genotypes$POS
saveRDS(genotypes_final, file = opt$filteredGenotypes, compress = T)
