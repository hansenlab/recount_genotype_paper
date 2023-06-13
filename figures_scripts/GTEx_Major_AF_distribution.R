library(tidyverse)
library(GenomicRanges)
library(cowplot)
theme_set(theme_cowplot())

plot_file = "../ready_to_plot/GTEx_Major_AF_distribution.rds"

if(!file.exists(plot_file)) {
    snps = readRDS("/dcl01/hansen/data/arazi/gtex_count/genotype/biallelic_SNP_gr.rda")

    snps$major_AF = case_when(snps$allele_freq <= .5  ~ 1 - snps$allele_freq,
                              snps$allele_freq > .5  ~ snps$allele_freq)

    test_snps = read.delim("/dcs04/hansen/data/recount_genotype/pipeline/GTEx_Blended_Tissue_Testing/all_tissues/GTEX_blended_model/aggregate_commonSNPs/commonSNPs.txt", header=F, sep = "\t")
    test_snps_gr = GRanges(seqnames = test_snps$V1,
                                ranges = IRanges(test_snps$V2,
                                                 test_snps$V2))
    ov = findOverlaps(test_snps_gr, snps)

    test_snps$major_AF = snps$major_AF[subjectHits(ov)]
    saveRDS(test_snps, file = plot_file)
}else {
    test_snps = readRDS(plot_file)
}


pdf("../figures/GTEx_Major_AF_distribution.pdf", width = 5, height = 4)

ggplot(test_snps, aes(x = major_AF)) + geom_histogram(bins = 25) + labs(x = "Major Allele Fraction")

dev.off()

print(summary(test_snps$major_AF))