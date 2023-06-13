library(optparse)

option_list <- list(
    make_option(c("-m", "--metadata"), type = "character",
    help = "Input: CSV file of Metadata. Contains columns: sample_id, total,
    alt, study. Each row is a sample."),
    make_option(c("-t", "--study"), type = "character",
    help = "Input: study name to run analysis on. Must be a study name
    found in Metadata specified."),
    make_option(c("-d", "--tempFolder"), type = "character",
    help = "Input: Temp directory to output bigWig alt files. They
    will be deleted at the end of the script."),
    make_option(c("-s", "--coverage"), type = "character",
    help = "Input: R object (dataframe): Concatenated BigWig dataframe filtered
    by SNP database. For a set samples from a study type."),
    make_option(c("-f", "--filteredSNPs"), type = "character",
    help = "Input: R object (GRanges): SNPs that intersected with BigWigs."),
    make_option(c("-r", "--ref"), type = "character",
    help = "Output: R object (data.table): SNP x samples data.table of
    ref counts."),
    make_option(c("-a", "--alt"), type = "character",
    help = "Output: R object (data.table): SNP x samples data.table of
    alt counts."),
    make_option(c("--major"), type = "character",
    help = "Output: R object (data.table): SNP x samples data.table of
    major allele counts."),
    make_option(c("--minor"), type = "character",
    help = "Output: R object (data.table): SNP x samples data.table of
    minor allele counts."))

opt <- parse_args(OptionParser(option_list = option_list))
if (length(opt$metadata) == 0 | length(opt$study) == 0 |
    length(opt$coverage) == 0 | length(opt$filteredSNPs) == 0 |
    length(opt$ref) == 0 | length(opt$alt) == 0 | length(opt$tempFolder) == 0 |
    length(opt$major) == 0 | length(opt$minor) == 0 ) {
    stop("Not all arguments provided. Check --help for description.")
}

library(data.table)
library(vroom)
library(dplyr)
library(rtracklayer)
library(GenomicRanges)

cat("study name: ", opt$study, "\n")
temp_folder <- paste0(opt$tempFolder, "/")
ptm_all <- proc.time()

#Load in Metadata, Chromosome contig, filtered SNPs.
metadata <- read.csv(opt$metadata, stringsAsFactors = FALSE)
metadata <- metadata[metadata$study == opt$study ,]
if (nrow(metadata) < 2)
    stop("Error: Metadata found < 2 samples given specified study name.")
chr_mapping <- fread("/dcl02/lieber/ajaffe/recount-pump/recount3.alts.chromosome_mappings.tsv", header = FALSE)
snps_gr <- readRDS(opt$filteredSNPs)
coverage_mtx <- as.data.table(readRDS(opt$coverage))

#Create `alt_mtx` data.table. 
alt_mtx <- as.data.table(matrix(0, nrow = length(snps_gr), 
                                   ncol = nrow(metadata)))
colnames(alt_mtx) <- metadata$sample_id_rep


#Set up identifying keys for `coverage_mtx` and `alt` data.tables. 
snps_key <- paste(as.character(seqnames(snps_gr)), start(snps_gr), 
    snps_gr$ref_seq, snps_gr$alt_seq, sep = "_")

print(gc())

for (i in seq_len(nrow(metadata))) {
    ptm <- proc.time()
    cat(i, " out of ", nrow(metadata), "\n")
    
    #Load in alt file to construct `alt_count`:
    #1. Decompress .zst via zstdcat in command line, save output as *.alt.temp.csv
    #2. Read in *.alt.temp.csv via vroom (fast), convert to data.table in order 
    #to collapase alt counts. 
    #3. Collapse to `alt_count`, and create GRanges object for intersection. 
    cat("Reading in: ", metadata$alt[i], "\n")
    temp_altFile <- paste0(temp_folder, metadata$sample_id_rep[i], ".alt.temp.csv")
    system(paste0("zstdcat ", metadata$alt[i], " > ", temp_altFile))
    alt <- as.data.table(vroom(temp_altFile, col_names = F, 
        col_select = c(1, 2, 4), col_types = c("i", "i", "c"), num_threads = 1))
    system(paste0("rm ", temp_altFile))
    colnames(alt) <- c("chr", "pos", "alt")
    #collapse alt counts. 
    alt_count <- alt[, .N, by = .(chr, pos, alt)] 
    rm(alt)
    alt_count = alt_count[!is.na(alt_count$alt) ,]
    #`alt_count` is 0-based, so we add 1 to positions. 
    alt_count$pos <- alt_count$pos + 1
    #use chromosome names that are used in recount3. 
    alt_count$chr <- chr_mapping$V2[match(alt_count$chr, chr_mapping$V1)] 
    alt_count_gr <- GRanges(seqnames = alt_count$chr,
                            ranges = IRanges(alt_count$pos,
                                             alt_count$pos))

    #Get `alt_count` into `alt_mtx`:
    ov <- findOverlaps(alt_count_gr, snps_gr)
    #Subset `alt_count_gr`, `alt_count` to the positions from SNPs, and compute 
    #its own unique key.
    alt_count_gr <- alt_count_gr[queryHits(ov)]
    alt_count <- alt_count[queryHits(ov)]
    alt_key <- paste(as.character(seqnames(alt_count_gr)), start(alt_count_gr), 
        snps_gr$ref_seq[subjectHits(ov)], alt_count$alt, sep = "_")
    #Even though `alt_key` is in the same positions as `snps_key`, many of the 
    #`alt_key` entries refer to alternate alleles that we are not tracking. 
    #We need to match a second time, this time using the entire key.
    idx <- match(snps_key, alt_key)
    snps_key_idx <- which(!is.na(idx))
    alt_key_idx <- idx[!is.na(idx)]
    alt_mtx[snps_key_idx, (i) := alt_count[alt_key_idx]$N]
    
    print(gc())
    print(proc.time() - ptm)
}

ref_mtx <- coverage_mtx - alt_mtx
colnames(ref_mtx) <- metadata$sample_id_rep
#We sometimes have cases where the alt counts > coverage counts (bigWig): 
#the alt reads were processed to keep overlapping pair-end reads
#whereas the coverage counts (bigWig) did not keep overlapping pair-end reads. 
#this is an ad hoc way to deal with negative counts in `ref_mtx`.
ref_mtx[ref_mtx < 0] <- 0 


#Finally, we also generate major and minor allele count matricies. 
#Major: count of major allele, which is determined from population allele fration of GTEx. 
#Minor: count of minor allele, which is determined from population allele fration of GTEx. 
#We first dupliicate `ref_mtx` and `alt_mtx' as `major_mtx` and `minor_mtx` respectively.
#Then, we look at the AF column in GTEx VCF and extracted SNPs in which AF > .5 
#(reference allele is *not* major allele). In `major_mtx` and `minor_mtx`,  
#we swap counts for these SNPs.

majorSNPs_gr <- snps_gr[snps_gr$allele_freq > .5]
ov <- findOverlaps(snps_gr, majorSNPs_gr)
major_mtx <- data.table::copy(ref_mtx)
minor_mtx <- data.table::copy(alt_mtx)
major_mtx[queryHits(ov) ,] <- alt_mtx[queryHits(ov) ,]
minor_mtx[queryHits(ov) ,] <- ref_mtx[queryHits(ov) ,]

saveRDS(ref_mtx, file = opt$ref, compress = T)
saveRDS(alt_mtx, file = opt$alt, compress = T)
saveRDS(major_mtx, file = opt$major, compress = T)
saveRDS(minor_mtx, file = opt$minor, compress = T)

cat("\nTotal time for ", nrow(metadata), "files: \n")
print(proc.time() - ptm_all)
