#This script takes crawls through recount3 outputs of GTEx data
#and pulls out bigwig file path, alt count file path, tissue, 
#individual_id, sample_id, replicate information.

#It filters out samples with missing tissues, and reformats
#for Snakemake pipeline. The metadata should have ~18k samples.

library(tidyverse)

#define constant varaibles of where bigwigs and alts are located. 
bigwig_path <- "/dcl02/lieber/ajaffe/recount-pump/gtex_monorail_output/bigwigs/"
alt_path <- "/dcl02/lieber/ajaffe/recount-pump/gtex_monorail_output/alts/"
tissue_info_path <- "/users/swang1/gtex_count/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"

#get bigwig files
bigwig_files <- system(paste0("/bin/ls ", bigwig_path, " | grep .all.bw"), intern = T)
path_total_bw <- paste0(bigwig_path, bigwig_files)
#get alt files
sample_rep <- sapply(strsplit(bigwig_files, "[.]"), function(x) x[2])
sample <- sub(".all.bw", "", bigwig_files)
sample_id <- sapply(strsplit(bigwig_files, "[.]"), function(x) x[1])
path_alts_bw <- paste0(alt_path, sample, ".bamcount_nonref.csv.zst")

#make a data table with all the sample names and file locations 
metadata <- data.frame(sample_id = sample_id, rep_id=sample_rep, 
					   total= path_total_bw, alt = path_alts_bw)

#get tissue names
sample_attribute <- read.delim(tissue_info_path, sep = "\t", 
						       stringsAsFactors=F)
metadata$tissue <- as.character(sapply(1:nrow(metadata), 
	function(x) sample_attribute$SMTSD[sample_attribute$SAMPID == metadata$sample_id[x]]))

#get individual ID
metadata$individual_id <- unlist(lapply(strsplit(metadata$sample_id, "-"), 
								function(x) paste0(x[1], "-", x[2])))

#add column for combination of sample_id and replicates:
metadata$sample_id_rep <- paste0(metadata$sample_id, ".", 
								 metadata$rep_id)

#filter out missing tissue samples
metadata <- metadata[metadata$tissue != "character(0)" ,]

#change tissue name to have no "-" or " " or "(" or ")" characters.
#those characters will crash Snakemake. 
metadata$tissue <- gsub("\\s*\\([^\\)]+\\)", "", 
						metadata$tissue) #removes ()
metadata$tissue <- trimws(metadata$tissue)
metadata$tissue <- gsub(" - ", "_", metadata$tissue)
metadata$tissue <- gsub(" ", "_", metadata$tissue)

#create study column: the column Snakemake will group each 
#study/tissue in its data processing.
metadata$study <- metadata$tissue

#now, generate metadatas for snakemake
write.csv(metadata, "GTEx_metadata.csv", quote = F, row.names = F)


#generate training set with blended tissues:
# We make sure that all samples in the training set are genotyped.
genotyped_individuals <- readLines("/dcl01/hansen/data/gtex_private/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.header.txt")
genotyped_individuals <- unlist(strsplit(genotyped_individuals[length(genotyped_individuals)], split = "\t"))
genotyped_individuals <- genotyped_individuals[startsWith(genotyped_individuals, "GTEX")]
genotyped_metadata <- metadata[metadata$individual_id %in% genotyped_individuals ,]

#consider tissues that have many individuals > 200
tissue_by_individuals <- genotyped_metadata %>% group_by(tissue) %>% 
			summarise(n_individuals = length(unique(individual_id)))
blended_tissues <- tissue_by_individuals$tissue[tissue_by_individuals$n_individuals >= 200]
blended_tissues_metadata <- genotyped_metadata %>% filter(tissue %in% blended_tissues)

#each individual can only contribue one tissue type to the training:
set.seed(1999)
training_metadata <- data.frame()
for(individual in unique(blended_tissues_metadata$individual_id)) {
	individual_metadata <- blended_tissues_metadata %>% filter(individual_id == individual)
	individual_metadata <- individual_metadata[sample(1:nrow(individual_metadata), 1) ,]
	training_metadata <- rbind(training_metadata, individual_metadata)
}

#leave out 200 individuals in the training set for testing.
leave_out_training <- sample(1:nrow(training_metadata), 200)
testing_individuals <- training_metadata$individual_id[leave_out_training]
training_metadata <- training_metadata[-leave_out_training ,]
training_metadata$study <- "GTEx_Blended_Tissue_Training"
training_individuals <- training_metadata$individual_id
write.csv(training_metadata, "GTEx_Blended_Tissue_Training_metadata.csv", quote = F, row.names = F)
saveRDS(training_individuals, file = "GTEx_Blended_Tissue_Training_individuals.rds")
#write out Test metadata
testing_metadata <- genotyped_metadata[genotyped_metadata$individual_id %in% testing_individuals ,]
write.csv(testing_metadata, "GTEx_Blended_Tissue_Testing_metadata.csv", quote = F, row.names = F)
