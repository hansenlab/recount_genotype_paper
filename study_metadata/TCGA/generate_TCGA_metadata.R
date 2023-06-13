library(tidyverse)
library(recount3)

#pull out TCGA paths from cluster
metadata_source = "/dcl02/lieber/ajaffe/recount-pump/tcga_monorail_output/recount_pump_output/all_bws2.all_bw.paths"
bigwigs = readLines(metadata_source)
bigwigs = sapply(bigwigs, function(x) substr(x, 2, nchar(x)))
names(bigwigs) = NULL
gdc_file_id = sapply(bigwigs, function(x) strsplit(x, "/")[[1]][5])
names(gdc_file_id) = NULL
metadata = data.frame(gdc_file_id = gdc_file_id,
					  total = paste0("/dcl02/lieber/ajaffe/recount-pump/tcga_monorail_output/recount_pump_output", bigwigs))
metadata$alt = gsub("all.bw", "bamcount_nonref.csv.zst", metadata$bigwig)


#match gdc_file_id with TCGA gdc_file_id via recount3 package metadata.
metadata$sample_id = NA
metadata$study = NA

human_projects <- available_projects()
tcga_projects = human_projects %>% filter(file_source == "tcga")

for(i in 1:nrow(tcga_projects)) {
	cat(i, " out of ", nrow(tcga_projects), "\n")
	cat(tcga_projects$project[i], "\n")
	tcga <- create_rse(tcga_projects[i ,], type ='gene')
	tcga_metadata = as.data.frame(colData(tcga)) %>% select(tcga.gdc_file_id, tcga.tcga_barcode)
	match_id = match(tcga_metadata$tcga.gdc_file_id, metadata$gdc_file_id)
	stopifnot(all(!is.na(match_id)))
	metadata$sample_id[match_id] = tcga_metadata$tcga.tcga_barcode
	metadata$study[match_id] = tcga_projects$project[i]
}

#check whether it is a tumor or normal sample
tumor_or_normal = function(x) {
	#based on https://docs.gdc.cancer.gov/Encyclopedia/pages/TCGA_Barcode/
	if(is.na(x)) {
		return(NA)
	}
	sample_vial = strsplit(x, "-")[[1]][4]
	sample = as.numeric(substr(sample_vial, 1, 2))
	if(sample <= 9) {
		return("tumor")
	}else if (sample >= 10 & sample <= 19) {
		return("normal")
	}else if (sample >= 20 & sample <= 29) {
		return("control")
	}else {
		return(NA)
	}
}
metadata$tumor_normal = sapply(metadata$sample_id, tumor_or_normal)

#save
write.csv(metadata, file = "tcga.tumor_and_normal.metadata.csv", quote=F, row.names=F)
write.csv(metadata %>% filter(tumor_normal == "normal"), file = "tcga.normal.metadata.csv", quote=F, row.names=F)