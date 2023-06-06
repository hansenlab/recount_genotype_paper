library(optparse)

option_list <- list(
  make_option(c("-g", "--genotypedSampleComplete"), type = "character",
              help = "Input: CSV file of genotyped sample."),
  make_option(c("-m", "--accuracyModelLattice"), type = "character",
              help = "Input: Accuracy model lattice, in .rds format"),
  make_option(c("-o", "--output"), type = "character",
              help = "Output: genotyped sample with accuracy predictions, in CSV format."))

opt <- parse_args(OptionParser(option_list = option_list))

if (length(opt$genotypedSampleComplete) == 0 | length(opt$accuracyModelLattice) == 0 |
    length(opt$output) == 0) {
  stop("Not all arguments provided. Check --help for description.")
}

library(tidyverse)
library(data.table)
library(rms) #restricted cubic splines

expit <- function(x) {
  return(1/(1+exp(-x)))
}

get_low_major_AF_quantile <- function(accuracyModelLattice) {
  major_AF_info <- unique(accuracyModelLattice$majorAF_bin)
  major_AF_info <- lapply(str_split(major_AF_info, ","), function(x) x[1]) #string manipulations to get it into a numeric vector
  major_AF_info <-  as.numeric(substr(major_AF_info, 2, 999))
  return(major_AF_info)
}

totalPtm <- proc.time()


accuracyModelLattice <- readRDS(opt$accuracyModelLattice)

eval_data <- fread(opt$genotypedSampleComplete)
eval_data$major_AF <- case_when(eval_data$AF <= .5  ~ 1 - eval_data$AF,
                                eval_data$AF > .5  ~ eval_data$AF)

#prediction for low major AF
eval_low_majorAF <- eval_data %>% filter(major_AF < .95)
eval_low_majorAF$majorAF_bin <- cut(eval_low_majorAF$major_AF, 
                                   breaks = get_low_major_AF_quantile(accuracyModelLattice), 
                                   include.lowest=TRUE)
eval_low_majorAF <- left_join(eval_low_majorAF, accuracyModelLattice, by = c("coverage", "majorAF_bin"))
eval_low_majorAF <- eval_low_majorAF %>% select(c(names(eval_data), "predicted_accuracy"))

#prediction for high major AF
eval_high_majorAF <- eval_data %>% filter(major_AF >= .95)
eval_high_majorAF$majorAF_bin <- "[0.95, 1]"
eval_high_majorAF <- left_join(eval_high_majorAF, accuracyModelLattice, by = c("coverage", "majorAF_bin"))
eval_high_majorAF <- eval_high_majorAF %>% select(c(names(eval_data), "predicted_accuracy"))

#put the two together
eval_data <- rbind(eval_low_majorAF, eval_high_majorAF)
eval_data$predicted_accuracy[is.na(eval_data$predicted_accuracy)] <- 1 #for coverage > 100 outside of lattice, give it perfect prediction
eval_data <- eval_data %>% select(-major_AF)

#save
fwrite(eval_data, file = opt$output)

cat("Total time elapsed:\n")
print(proc.time() - totalPtm)
print(gc())
