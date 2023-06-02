
library(optparse)

option_list <- list(
    make_option(c("-m", "--metadata"), type = "character",
    help = "Input: CSV file of Metadata. Contains columns: sample_id, total,
    alt, study Each row is a sample."),
    make_option(c("-t", "--study"), type = "character",
    help = "Input: Study name to run analysis on. Must be a study name
    found in Metadata specified."),
    make_option(c("-d", "--tempFolder"), type = "character",
    help = "Input: Temp directory to output bigWig coverage files. They
    will be deleted at the end of the script."),
    make_option(c("-s", "--SNPs"), type = "character", help = "Input:
    SNP locations in .rds (GRanges object)."),
    make_option(c("--coverage"), type = "character",
    help = "Output: R object (dataframe): Concatenated BigWig dataframe filtered
    by SNP database. For a set samples from a study type."),
    make_option(c("-f", "--filteredSNPs"), type = "character",
    help = "Output: R object (GRanges): SNPs that intersected with BigWigs."),
    make_option(c("-b", "--filteredSNPsList"), type = "character",
    help = "Output: text file: SNPs that intersected with BigWigs.
            First column: chr, Second column: position."))

opt <- parse_args(OptionParser(option_list = option_list))
if (length(opt$metadata) == 0 | length(opt$study) == 0 |
    length(opt$tempFolder) == 0 | length(opt$SNPs) == 0 |
    length(opt$coverage) == 0 | length(opt$filteredSNPs) == 0 |
    length(opt$filteredSNPsList) == 0) {
    stop("Not all arguments provided. Check --help for description.")
}


library(rtracklayer)
library(GenomicRanges)

coverage_threshold = 5

cat("Study name: ", opt$study, "\n")
ptm_all <- proc.time()
temp_folder <- paste0(opt$tempFolder, "/")

#load in Metadata and SNP GRanges.
metadata <- read.csv(opt$metadata, stringsAsFactors = FALSE)
metadata <- metadata[metadata$study == opt$study ,]
if (nrow(metadata) < 2)
    stop("Error: Metadata found < 2 samples given specified study name.")

biallelic_SNP_gr <- readRDS(opt$SNPs)

#online vector for computing mean coverage per site
mean_coverage <- numeric(length(biallelic_SNP_gr))

#First pass through files: read in bigWig, intersect with SNPs, and write
#intersected coverage counts in a temp file. Keep track of total coverage.
for (i in seq_len(nrow(metadata))) {
    ptm <- proc.time()
    cat(i, " out of ", nrow(metadata), "\n")
    cat("Reading in: ", metadata$total[i], "\n")

    bigwig_i <- import(metadata$total[i], format = "bigwig")
    overlap_loci <- findOverlaps(biallelic_SNP_gr, bigwig_i)
    bigwig_counts_i <- bigwig_i$score[subjectHits(overlap_loci)]

    mean_coverage <- mean_coverage + bigwig_counts_i

    temp_file <- file(paste0(temp_folder, metadata$sample_id_rep[i], ".tmp"), "wb")
    writeBin(bigwig_counts_i, temp_file)
    close(temp_file)

    print(gc())
    print(proc.time() - ptm)
}

#finish computing mean coverage
mean_coverage <- mean_coverage / nrow(metadata)
filter_index <- which(mean_coverage > coverage_threshold)

cat("Number of sites kept after filter: ", length(filter_index), "\n")

#Second pass through files: read in temp file, subset to sites that pass filter,
#delete temp file, and keep site counts to construct coverage matrix.
temp_bw <- lapply(seq_len(nrow(metadata)), function(i) {
    ptm <- proc.time()
    file_path_i <- paste0(temp_folder, metadata$sample_id_rep[i], ".tmp")
    cat("Reading in: ", file_path_i, "\n")

    bigwig_counts_i <- readBin(file_path_i,
                               what = "numeric",
                               n = length(biallelic_SNP_gr))
    file.remove(file_path_i) #delete temp file
    print(gc())
    print(proc.time() - ptm)
    return(bigwig_counts_i[filter_index])
})

#construct coverage matrix.
coverage_mtx <- do.call(cbind, temp_bw)
colnames(coverage_mtx) <- metadata$sample_id_rep

#filter original list of SNPs GRanges
biallelic_SNP_gr_filtered <- biallelic_SNP_gr[filter_index]
biallelic_SNP_gr_filtered_list <- data.frame(chromosome = seqnames(biallelic_SNP_gr_filtered),
                                            start = start(biallelic_SNP_gr_filtered))

#save results
saveRDS(coverage_mtx, 
        file = opt$coverage, 
        compress = T)
saveRDS(biallelic_SNP_gr_filtered,
        file = opt$filteredSNPs,
        compress = T)
write.table(biallelic_SNP_gr_filtered_list,
            file = opt$filteredSNPsList,
            row.names = F,
            col.names = F,
            sep="\t",
            quote=F)

cat("\nTotal time for ", nrow(metadata), "files: \n")
print(proc.time() - ptm_all)
