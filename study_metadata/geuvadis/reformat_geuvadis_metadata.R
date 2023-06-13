
#https://www.ebi.ac.uk/arrayexpress/experiments/E-GEUV-1/samples/?full=true&s_page=25&s_pagesize=25&s_sortby=col_23&s_sortorder=ascending
EMBL_metadata = read.delim("input/EMBL_geuvadis_metadata.txt")

metadata = read.csv("input/geuvadis_SRA_metadata.csv")
vcf_individuals = readLines("input/individual_id_from_geuvadis_VCF.txt")
vcf_individuals = unlist(strsplit(vcf_individuals, "\t"))

found_vcf_individuals = match(vcf_individuals, EMBL_metadata$Source.Name)
found_vcf_individuals = found_vcf_individuals[!is.na(found_vcf_individuals)]
EMBL_metadata = EMBL_metadata[found_vcf_individuals ,]

metadata$individual_id = EMBL_metadata$Source.Name[match(metadata$sample_id, EMBL_metadata$Comment.ENA_RUN.)]

write.csv(metadata, "geuvadis_metadata.csv", row.names=F, quote=F)
