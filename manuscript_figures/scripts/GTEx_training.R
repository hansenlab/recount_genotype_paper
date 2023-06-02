library(tidyverse)
library(data.table)
library(cowplot)
library(gridExtra)
theme_set(theme_cowplot())

if(!file.exists("../ready_to_plot/GTEx_training.rds")) {
	#Load in trained model
	prior = readRDS("/dcs04/hansen/data/recount_genotype/training/GTEx_Blended_Tissue_Training/p_Zg_GTEx_Blended_Tissue_Training.rds")
	model = readRDS("/dcs04/hansen/data/recount_genotype/training/GTEx_Blended_Tissue_Training/model_MS_GTEx_Blended_Tissue_Training.rds")
	lattice = readRDS("/dcs04/hansen/data/recount_genotype/training/GTEx_Blended_Tissue_Training/model_MS_GTEx_Blended_Tissue_Training_lattice.rds")
	metadata = read.csv("../../training_snakemake/GTEx_Blended_Tissue_Training_metadata.csv")
	#Load in ref, alt, genotype information
	ref_matrix <- readRDS("/dcs04/hansen/data/recount_genotype/pipeline/GTEx_Blended_Tissue_Training/get_alt_ref_matrix/GTEx_Blended_Tissue_Training_ref.rds")
	alt_matrix <- readRDS("/dcs04/hansen/data/recount_genotype/pipeline/GTEx_Blended_Tissue_Training/get_alt_ref_matrix/GTEx_Blended_Tissue_Training_alt.rds")
	snp_genotype <- readRDS("/dcs04/hansen/data/recount_genotype/pipeline/GTEx_Blended_Tissue_Training/get_genotype_matrix/GTEx_Blended_Tissue_Training_genotype.rda")
	#Load in unique SNP ids
	snps_loci = read.delim("/dcs04/hansen/data/recount_genotype/pipeline/GTEx_Blended_Tissue_Training/get_coverage_matrix/GTEx_Blended_Tissue_Training_filteredSNPs.txt", sep = "\t", stringsAsFactors = F, header = F)
	snps_loci = paste0(snps_loci$V1, "_", snps_loci$V2)

	#M, S transformation
	M = log2((ref_matrix + 1) / (alt_matrix + 1))
	S = log2(sqrt((ref_matrix + 1) * (alt_matrix + 1)))

	#Map snp_genotype into numerical values.
	snp_genotype[snp_genotype == "0/0"] <- 0 #Homozygous ref
	snp_genotype[snp_genotype== "0/1" | snp_genotype == "1/0" ] <-1 #Heterozygous
	snp_genotype[snp_genotype == "1/1"] <- 2 #Homozygous alt
	snp_genotype[snp_genotype == "./."] <- 3 # missing

	#Sample random SNPs for plotting.
	set.seed(123)
	SNP_sample_idx = sample(1:nrow(M), 1000)
	ref_samp = ref_matrix[SNP_sample_idx ,]
	ref_samp$chr_pos = snps_loci[SNP_sample_idx]
	alt_samp = alt_matrix[SNP_sample_idx ,]
	alt_samp$chr_pos = snps_loci[SNP_sample_idx]
	M_samp = M[SNP_sample_idx ,]
	M_samp$chr_pos = snps_loci[SNP_sample_idx]
	S_samp = S[SNP_sample_idx ,]
	S_samp$chr_pos = snps_loci[SNP_sample_idx]
	snp_samp = snp_genotype[SNP_sample_idx ,]
	snp_samp$chr_pos = snps_loci[SNP_sample_idx]

	#Merge all datasets together and rename labels. 
	ref_samp_long = ref_samp %>% pivot_longer(!chr_pos, names_to = "sample", values_to = "ref")
	alt_samp_long = alt_samp %>% pivot_longer(!chr_pos, names_to = "sample", values_to = "alt")
	M_samp_long = M_samp %>% pivot_longer(!chr_pos, names_to = "sample", values_to = "M")
	S_samp_long = S_samp %>% pivot_longer(!chr_pos, names_to = "sample", values_to = "S")
	snp_samp_long = snp_samp %>% pivot_longer(!chr_pos, names_to = "sample", values_to = "genotype")
	to_merge = list(ref_samp_long, alt_samp_long, M_samp_long, S_samp_long, snp_samp_long)
	M_S_snp_samp = Reduce(inner_join, to_merge)
	M_S_snp_samp = M_S_snp_samp[M_S_snp_samp$genotype != 3 ,]
	M_S_snp_samp$tissue = metadata$tissue[match(M_S_snp_samp$sample, metadata$sample_id_rep)]
	M_S_snp_samp$genotype = gsub(0, "Ref Homozygous", M_S_snp_samp$genotype)
	M_S_snp_samp$genotype = gsub(1, "Heterozygous", M_S_snp_samp$genotype)
	M_S_snp_samp$genotype = gsub(2, "Alt Homozygous", M_S_snp_samp$genotype)
	M_S_snp_samp$genotype <- factor(M_S_snp_samp$genotype, ordered = TRUE, levels = c("Ref Homozygous", "Heterozygous", "Alt Homozygous"))

	#Construct mu_M, mu_S, sd_M
	M_S_snp_samp_mu = M_S_snp_samp %>% group_by(chr_pos, genotype, tissue) %>% 
										summarize(mu_M = mean(M), mu_S = mean(S), sd_M = sd(M), nSamples = n())


	#Construct decision boundary grid.
	lattice = expand.grid(S = seq(0, 10, .1), M = seq(-10, 10, .1))
	lattice$mu_homRef <- predict(model[[1]][[1]], data.frame(mu_S = lattice$S))
	lattice$sd_homRef <- predict(model[[1]][[2]], data.frame(mu_S = lattice$S))
	lattice$mu_het <- predict(model[[2]][[1]], data.frame(mu_S = lattice$S))
	lattice$sd_het <- predict(model[[2]][[2]], data.frame(mu_S = lattice$S))
	lattice$mu_homAlt <- predict(model[[3]][[1]], data.frame(mu_S = lattice$S))
	lattice$sd_homAlt <- predict(model[[3]][[2]], data.frame(mu_S = lattice$S))

	lattice$dens_1 <- log(prior[1]) + dnorm(x = lattice$M, mean = lattice$mu_homRef, sd = lattice$sd_homRef, log = TRUE)
	lattice$dens_2 <- log(prior[2]) + dnorm(x = lattice$M, mean = lattice$mu_het, sd = lattice$sd_het, log = TRUE)
	lattice$dens_3 <- log(prior[3]) + dnorm(x = lattice$M, mean = lattice$mu_homAlt, sd = lattice$sd_homAlt, log = TRUE)

	lattice$genotype = mapply(function(x, y, z) {
			a = exp(c(x, y, z))
			return(which(a == max(a)))}, 
		lattice$dens_1, lattice$dens_2, lattice$dens_3)
	lattice$genotype = gsub(1, "Ref Homozygous", lattice$genotype)
	lattice$genotype = gsub(2, "Heterozygous", lattice$genotype)
	lattice$genotype = gsub(3, "Alt Homozygous", lattice$genotype)
	lattice$genotype <- factor(lattice$genotype, ordered = TRUE, levels = c("Ref Homozygous", "Heterozygous", "Alt Homozygous"))

	save(M_S_snp_samp, M_S_snp_samp_mu, lattice, file = "../ready_to_plot/M_S_decisionBoundary.rds")
}else {
	load("../ready_to_plot/GTEx_training.rds")
}

levels(M_S_snp_samp$genotype)<- c("AA", "AB", "BB")
levels(lattice$genotype)<- c("AA", "AB", "BB")
levels(M_S_snp_samp_mu$genotype)<- c("AA", "AB", "BB")

M_S_snp_samp = M_S_snp_samp[sample(1:nrow(M_S_snp_samp), .2 * nrow(M_S_snp_samp)) ,]

#Figure S1
pdf("../figures/GTEx_training_ref_alt.pdf", width = 7, height = 2.75)
ggplot(M_S_snp_samp, aes(y = ref, x = alt)) + geom_point(alpha = .1) + facet_wrap(~genotype) + scale_x_continuous(limits = c(0, 500)) + scale_y_continuous(limits = c(0, 500))+ labs(y="Reference", x="Alternative")
dev.off()

pdf("../figures/GTEx_training_S_M.pdf", width = 7, height = 2.75)
ggplot(M_S_snp_samp, aes(x = S, y = M)) + geom_point(alpha = .1) + facet_wrap(~genotype) + scale_x_continuous(limits = c(0, 8)) + scale_y_continuous(limits = c(-10, 10))
dev.off()


#Figure 2

pdf("../figures/GTEx_training.pdf", width = 7, height = 2.75*3)
p1 = ggplot(M_S_snp_samp_mu, aes(x = mu_S, y = mu_M)) + 
	geom_point(alpha = .1) + 
	facet_wrap(~genotype) + 
	geom_smooth(method = "lm", linewidth = 1) + 
	scale_x_continuous(limits = c(1.161, 8)) +
	scale_y_continuous(limits = c(-10, 10)) + 
	labs(x="mean(S)", y="mean(M)") +
	theme(axis.title.y = element_text(margin = margin(t = 0, r = -7, b = 0, l = 0)))

p2 = ggplot(M_S_snp_samp_mu, aes(x = mu_S, y = sd_M)) + geom_point(alpha = .1) + 
	facet_wrap(~genotype) + 
	geom_smooth(method = "gam", linewidth = 1) + 
	scale_x_continuous(limits = c(1.161, 8)) +
	labs(x="mean(S)", y= "sd(M)") +
	theme(axis.title.y = element_text(margin = margin(t = 0, r = 0, b = 0, l = 5)))


p3 = ggplot(lattice) + geom_point(aes(x = S, y = M), color = "dark blue", alpha = .02) + 
	facet_wrap(~genotype) + 
	scale_x_continuous(limits = c(1.161, 8)) +
 	scale_y_continuous(limits = c(-10, 10)) +
	geom_point(data = M_S_snp_samp, aes(x = S, y = M), alpha = .1) +
	theme(legend.position = "none") +
	theme(axis.title.y = element_text(margin = margin(t = 0, r = -7, b = 0, l = 0)))


grid.arrange(p1, p2, p3, nrow = 3)
dev.off()