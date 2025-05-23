library(optparse)

option_list <- list(
	make_option(c("--sample_id_rep"), type = "character",
    help = "Input: String: sample_id replicate"),
    make_option(c("--SNPs"), type = "character",
    help = "Input: R object: GRanges of SNPs."),
    make_option(c("--total"), type = "character",
    help = "Input: bigWig file path."),
    make_option(c("--alt"), type = "character",
    help = "Input: alt file path."),
    make_option(c("--cutoff"), type = "numeric",
    help = "Input numeric: Cutoff for coverage filtering."),
    make_option(c("--modelLattice"), type = "character",
    help = "Input: R object: Mean and variance model in lattice format."),
	  make_option(c("--model"), type = "character",
	  help = "Input: R object: Path to our genotyping model"),
    make_option(c("--prior"), type = "character",
    help = "Input: R object: Estimation of Pi prior for model."),
    make_option(c("--tempFolder"), type = "character",
    help = "Input: temp folder path."),
    make_option(c("--result"), type = "character", 
    help = "Output csv: chr, start, AF, M, S, pred_genotype"),
    make_option(c("--outOfLattice"), type = "character", 
    help = "Output csv: Values outside of lattice: chr, start, AF, M, S, pred_genotype=NA"))

opt <- parse_args(OptionParser(option_list = option_list))

if (length(opt$sample_id_rep) == 0 |
	length(opt$SNPs) == 0 |
	length(opt$cutoff) == 0 |
	length(opt$total) == 0 |
	length(opt$alt) == 0 |
    length(opt$modelLattice) == 0 | 
    length(opt$model) == 0 |
    length(opt$prior) == 0 |
    length(opt$tempFolder) == 0 |
    length(opt$outOfLattice) == 0 |
    length(opt$result) == 0 ) {
        stop("Not all arguments provided. Check --help for description.")
}

totalPtm <- proc.time()
library(data.table)
library(dplyr)
library(GenomicRanges)
library(rtracklayer)
library(caret)

recount3_chr_mapping <- "/dcs04/hansen/data/recount3/recount3.alts.chromosome_mappings.tsv"


get_coverage_count <- function(snps_path, bigWig_path, coverage_cutoff) {
	#Load in bigWig file to get `coverage_count` and `filtered_snps_gr`. 
	snps_gr <- readRDS(snps_path)
	cat("Loading in: ", bigWig_path, "\n")
	
	 bigwig <- tryCatch(
        {
            import(bigWig_path, format = "bigwig")
        },
        error=function(cond) {
            message(paste("Error loading bigwig file: ", bigWig_path))
            message(cond)
            return(NA)
        },
        warning=function(cond) {
        	message(paste("Warning loading bigwig file: ", bigWig_path))
            message(cond)
            return(NA)
        },
        finally={}
    )    
	if(all(is.na(bigwig))) {
		return(NA)
	}
	overlap_loci <- findOverlaps(snps_gr, bigwig)
	bigwig_count <- bigwig$score[subjectHits(overlap_loci)]
	filter_idx <- which(bigwig_count >= coverage_cutoff)
	bigwig_count <- bigwig_count[filter_idx]
	filtered_snps_gr <- snps_gr[queryHits(overlap_loci)[filter_idx]]
	return(list(coverage_count = bigwig_count,
				filtered_snps_gr = filtered_snps_gr))
}

get_alt_count <- function(alt_path, sample_id_rep, filtered_snps_gr, temp_folder) {
	#Load in alt file to construct `alt_count`:
	#1. Decompress .zst via zstdcat in command line, save output as *.alt.temp.csv
	#2. Read in *.alt.temp.csv via vroom (fast), convert to data.table in order 
	#to collapase alt counts. 
	#3. Collapse to `alt_count`, and create GRanges object for intersection.
	cat("Loading in: ", alt_path, "\n")
	temp_altFile <- paste0(temp_folder, sample_id_rep, ".alt.temp.csv")
	system(paste0("zstdcat ", alt_path, " > ", temp_altFile))
	alt <- fread(temp_altFile, select = c(1, 2, 4))
	system(paste0("rm ", temp_altFile))
	if(nrow(alt) == 0) {
		return(NA)
	}
	colnames(alt) <- c("chr", "pos", "alt")
	#collapse alt counts. 
	alt_count <- alt[, .N, by = .(chr, pos, alt)] 
	rm(alt)
	alt_count <- alt_count[!is.na(alt_count$alt) ,]
	#`alt_count` is 0-based, so we add 1 to positions. 
	alt_count$pos <- alt_count$pos + 1
	#use chromosome names that are used in recount3. 
	chr_mapping <- fread(recount3_chr_mapping, header = FALSE)
	alt_count$chr <- chr_mapping$V2[match(alt_count$chr, chr_mapping$V1)] 
	alt_count_gr <- GRanges(seqnames = alt_count$chr,
	                        ranges = IRanges(alt_count$pos,
	                                         alt_count$pos))
	ov <- findOverlaps(alt_count_gr, filtered_snps_gr)
	#Subset `alt_count_gr`, `alt_count` to the positions from SNPs, and compute 
	#its own unique key.
	alt_count_gr <- alt_count_gr[queryHits(ov)]
	alt_count <- alt_count[queryHits(ov)]
	alt_key <- paste(as.character(seqnames(alt_count_gr)), start(alt_count_gr), 
	    filtered_snps_gr$ref_seq[subjectHits(ov)], alt_count$alt, sep = "_")
	#Even though `alt_key` is in the same positions as `filtered_snps_key`, many of the 
	#`alt_key` entries refer to alternate alleles that we are not tracking. 
	#We need to match a second time, this time using the entire key.
	filtered_snps_key <- paste(as.character(seqnames(filtered_snps_gr)), start(filtered_snps_gr), 
		filtered_snps_gr$ref_seq, filtered_snps_gr$alt_seq, sep = "_")
	idx <- match(filtered_snps_key, alt_key)
	snps_key_idx <- which(!is.na(idx))
	alt_key_idx <- idx[!is.na(idx)]
	#Construct final `alt_count` relative to `filtered_snps_gr`
	final_alt_count <- rep(0, length(filtered_snps_gr))
	final_alt_count[snps_key_idx] <- alt_count[alt_key_idx]$N

	
	#calculate the sequencing counts aligned to non alt/ref sequence
	#this should be subtracted from the total later
	not_alt_count<-alt_count[-alt_key_idx,]
	not_alt_count<-not_alt_count %>% group_by(chr,pos) %>% summarize(error=sum(N))
	error_gr<-makeGRangesFromDataFrame(not_alt_count,seqnames="chr",start.field ="pos",end.field = "pos", keep.extra.columns = T)
	
	ov <- findOverlaps(error_gr, filtered_snps_gr)
	final_error_count <- rep(0, length(filtered_snps_gr))
	final_error_count[subjectHits(ov)] <- error_gr$error[queryHits(ov)]
	
	return(list(final_alt_count = final_alt_count,
	            final_error_count = final_error_count))
	
	
	
	
}

predict_genotype <- function(model_path, model_lattice_path, prior_path, M, S) {
	lattice_max <- 8.5 
	#Load in model information.
	model_MS_lattice <- readRDS(model_lattice_path)
	model_MS <- readRDS(model_path)
	prior <- readRDS(prior_path)
	if(all(is.na(prior)) == T) {
	  prior = c(.94, .04, .02) 
	}
	#Figure out which values of predictor S can be used for lattice lookup table.
	withinLattice_idx <- which(S <= lattice_max)
	outsideLattice_idx <- which(S > lattice_max)
	cat("Proportion of S values that can be inferred using lattice: ", length(withinLattice_idx)/length(S), "\n")
	cat("Number of S values to infer with model later: ", length(outsideLattice_idx), "\n")
	#Compute likelihood * prior for each genotype.
	likelihood_times_prior <- sapply(1:3, function(k){
	  #Inference on mu_M and sd_M
	  mu_M <- rep(NA, length(withinLattice_idx))
	  sd_M <- rep(NA, length(withinLattice_idx))
	  if(k == 1) {
	    #model_MS_lattice is a data.table with key set to its S values, 
	    #so we use the key to get lattice values.
	    mu_M <- model_MS_lattice[J(S[withinLattice_idx])]$mu_M_prediction_1
	    sd_M <- model_MS_lattice[J(S[withinLattice_idx])]$sd_M_prediction_1
	    indx<-which(is.na(mu_M))
	    mu_M[indx] <- predict(model_MS[[1]][[1]], data.frame(mu_S = S[withinLattice_idx][indx]))
	    sd_M[indx] <- predict(model_MS[[1]][[2]], data.frame(mu_S = S[withinLattice_idx][indx]))
	  }else if(k == 2) {
	    mu_M <- model_MS_lattice[J(S[withinLattice_idx])]$mu_M_prediction_2
	    sd_M <- model_MS_lattice[J(S[withinLattice_idx])]$sd_M_prediction_2
	    indx<-which(is.na(mu_M))
	    mu_M[indx] <- predict(model_MS[[2]][[1]], data.frame(mu_S = S[withinLattice_idx][indx]))
	    sd_M[indx] <- predict(model_MS[[2]][[2]], data.frame(mu_S = S[withinLattice_idx][indx]))
	  }else if(k == 3) {
	    mu_M <- model_MS_lattice[J(S[withinLattice_idx])]$mu_M_prediction_3
	    sd_M <- model_MS_lattice[J(S[withinLattice_idx])]$sd_M_prediction_3
	    indx<-which(is.na(mu_M))
	    mu_M[indx] <- predict(model_MS[[3]][[1]], data.frame(mu_S = S[withinLattice_idx][indx]))
	    sd_M[indx] <- predict(model_MS[[3]][[2]], data.frame(mu_S = S[withinLattice_idx][indx]))
	  }
	  sd_M[sd_M <= 0] = 0 #if we predict (extrapolate) any variance to be <= 0, set it to 0. 
	  return(prior[k] * dnorm(x = M[withinLattice_idx], mean = mu_M, sd = sd_M, log = FALSE))
	})
	


  	posterior <- sapply(1:3, function(k){
        if(is.null(dim(likelihood_times_prior))) { #if we only have one SNP that needs to be genotyped
            return(likelihood_times_prior[k] / sum(likelihood_times_prior))
        }else { #else, we have many SNPs in a dataframe
            return(likelihood_times_prior[, k] / rowSums(likelihood_times_prior))
        }
    })

    #Pick genotype with the max posterior genotype as our prediction.
    if(is.null(dim(posterior))) { #if we only have one SNP that needs to be genotyped
        predicted_genotype <- which(posterior == max(posterior))
    }else { #else, we have many SNPs in a dataframe
        predicted_genotype <- apply(posterior, MARGIN = 1, FUN = function(x) which(x == max(x)))
    }

  	
	return(list(predicted_genotype = predicted_genotype,
				withinLattice_idx = withinLattice_idx,
				outsideLattice_idx = outsideLattice_idx))
}

writeNullOutput = function(resultFile, latticeFile) {
	#Write a blank output when there is no data to genotype. We do not throw an error
	#because it allows the pipeline to keep running, and it will be apparent when the user
	#see that there's nothing that was genotyped.
	cat("Error loading in bigwig or no data to genotype! Writing blank results.\n")
	result <- data.table(chr = character(), start = numeric(), AF = numeric(), M = numeric(), S = numeric(), coverage = numeric(),ref_count= numeric(), pred_genotype = numeric())
	fwrite(result, file = resultFile)
	fwrite(result, file = latticeFile)
}


cat("Loading in SNPs and bigWig.\n")
ptm <- proc.time()

coverage_count_result <- get_coverage_count(snps_path = opt$SNPs, 
											bigWig_path = opt$total, 
											coverage_cutoff = as.numeric(opt$cutoff))

print(gc())
print(proc.time() - ptm)
cat("Finished loading in SNPs and bigWig.\n")

if(any(is.na(coverage_count_result)) | length(coverage_count_result$coverage_count) == 0) {
	writeNullOutput(opt$result, opt$outOfLattice)
	q()
}

cat("Loading in alt counts.\n")
ptm <- proc.time()

alt_count_df <- get_alt_count(alt_path = opt$alt, 
						   sample_id_rep = opt$sample_id_rep, 
						   filtered_snps_gr = coverage_count_result$filtered_snps_gr,
						   temp_folder = paste0(opt$tempFolder, "/"))


print(gc())
print(proc.time() - ptm)
cat("Finished loading in alt counts.\n")

if(all(is.na(alt_count_df))) {
	writeNullOutput(opt$result, opt$outOfLattice)
	q()
}
alt_count <- alt_count_df$final_alt_count
total_cov <- coverage_count_result$coverage_count- alt_count_df$final_error_count
filtered_snps_gr <- coverage_count_result$filtered_snps_gr

#At some locations, the alt count is larger than the total count: in these SNPs replace the alt count with total count
id<-which(alt_count>total_cov)
alt_count[id]<-total_cov[id]
#Compute M and S values. 
ref_count <- total_cov - alt_count
#We sometimes have cases where the alt counts > coverage counts (bigWig): 
#the alt reads were processed to keep overlapping pair-end reads
#whereas the coverage counts (bigWig) did not keep overlapping pair-end reads. 
#this is an ad hoc way to deal with negative counts in `ref_count`.
error_indx<-which(ref_count < 0 | alt_count < 0 | total_cov< opt$cutoff)
if(length(error_indx)>0){
  alt_count<- alt_count[-error_indx] 
  total_cov<- total_cov[-error_indx] 
  ref_count<- ref_count[-error_indx] 
  filtered_snps_gr<- filtered_snps_gr[-error_indx]
} 
#To prevent divde by 0 in M, S calculation, we add psuedocount of 1. 
M <- log2((ref_count + 1) / (alt_count + 1))
S <- log2(sqrt((ref_count + 1) * (alt_count + 1)))

cat("Predicting genotype.\n")
ptm <- proc.time()

pred_genotype_result <- predict_genotype(model_path = opt$model, 
										 model_lattice_path = opt$modelLattice, 
										 prior_path = opt$prior,
										 M = M, 
										 S = S)

print(gc())
print(proc.time() - ptm)
cat("Finished predicting genotype. Saving results now.\n")

#Predicted genotype result.
result <- data.table(chr = as.character(seqnames(filtered_snps_gr[pred_genotype_result$withinLattice_idx])),
					start = start(filtered_snps_gr[pred_genotype_result$withinLattice_idx]),
					AF = round(filtered_snps_gr$allele_freq[pred_genotype_result$withinLattice_idx], 4),
					M = M[pred_genotype_result$withinLattice_idx],
					S = S[pred_genotype_result$withinLattice_idx],
					coverage = ref_count[pred_genotype_result$withinLattice_idx] + alt_count[pred_genotype_result$withinLattice_idx],
					ref_count= ref_count[pred_genotype_result$withinLattice_idx],
					pred_genotype = pred_genotype_result$predicted_genotype)
fwrite(result, file = opt$result)

#M, S values that fall outside of our lattice, to be predicted with model file later.
outOfLattice_result <- data.table(chr = as.character(seqnames(filtered_snps_gr[pred_genotype_result$outsideLattice_idx])),
								start = start(filtered_snps_gr[pred_genotype_result$outsideLattice_idx]),
								AF = round(filtered_snps_gr$allele_freq[pred_genotype_result$outsideLattice_idx], 4),
								M = M[pred_genotype_result$outsideLattice_idx],
								S = S[pred_genotype_result$outsideLattice_idx],
								coverage = ref_count[pred_genotype_result$outsideLattice_idx] + alt_count[pred_genotype_result$outsideLattice_idx],
								ref_count= ref_count[pred_genotype_result$outsideLattice_idx],
								pred_genotype = NA)
fwrite(outOfLattice_result, opt$outOfLattice)


cat("Total time elapsed:\n")
print(proc.time() - totalPtm)


