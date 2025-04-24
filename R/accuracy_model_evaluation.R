library(optparse)

option_list <- list(
  make_option(c("-a", "--allGenotypedSamplesAgg"), type = "character",
              help = "Input: Each SNP from samples in long format, with true and predicted genotypes 
              and alternative allele fraction, in .rds format."),
  make_option(c("-m", "--accuracyModel"), type = "character",
              help = "Input: Accuracy model, in .rds format"),
  make_option(c("-o", "--accuracyModelEvaluation"), type = "character",
              help = "Output: Evaluation metrics of accuracy model in coverage and MAF bins, in .rds format."),
  make_option(c("-e", "--accuracyModelEvaluationError"), type = "character",
              help = "Output: Evaluation metrics of accuracy model in coverage and MAF bins by mean absolute error, in .rds format."),
  make_option(c("-d", "--downSampling"), type = "numeric", default = 1,
              help = "Option: Sampling proportion for predition. Default: 1 (no sampling)."))

opt <- parse_args(OptionParser(option_list = option_list))

if (length(opt$allGenotypedSamplesAgg) == 0 | length(opt$accuracyModel) == 0 |
    length(opt$accuracyModelEvaluation) == 0 |  length(opt$accuracyModelEvaluationError) == 0 ) {
  stop("Not all arguments provided. Check --help for description.")
}

library(tidyverse)
library(data.table)
library(rms) #restricted cubic splines

accuracyModel <- readRDS(opt$accuracyModel)
low_major_AF_quantile <- accuracyModel$low_majorAF$Design$parms$majorAF_bin #the AF bins used in the accuracy model
low_major_AF_quantile <- lapply(str_split(low_major_AF_quantile, ","), function(x) x[1])
low_major_AF_quantile <-  as.numeric(substr(low_major_AF_quantile, 2, 999))
low_major_AF_quantile <- c(low_major_AF_quantile, .95)

eval_data <- fread(opt$allGenotypedSamplesAgg)
downSampling <- as.numeric(opt$downSampling)
eval_data = eval_data %>% filter(true_genotype != "./." & true_genotype != "." & true_genotype !="")
eval_data$correct = as.numeric(eval_data$pred_genotype == eval_data$true_genotype)
eval_data$major_AF = case_when(eval_data$AF <= .5  ~ 1 - eval_data$AF,
                               eval_data$AF > .5  ~ eval_data$AF)

eval_low_majorAF = eval_data %>% filter(major_AF < .95)
eval_low_majorAF$majorAF_bin = cut(eval_low_majorAF$major_AF, breaks = low_major_AF_quantile, include.lowest=TRUE)
eval_low_majorAF = eval_low_majorAF[sample(1:nrow(eval_low_majorAF), downSampling * nrow(eval_low_majorAF)) ,]
eval_low_majorAF$predicted.values.logit = predict(accuracyModel[["low_majorAF"]], newdata = as.data.frame(eval_low_majorAF))
eval_low_majorAF$predicted.values.prob = 1/(1+exp(-eval_low_majorAF$predicted.values.logit))

eval_high_majorAF = eval_data %>% filter(major_AF >= .95)
eval_high_majorAF$majorAF_bin = "[0.95, 1]"
eval_high_majorAF = eval_high_majorAF[sample(1:nrow(eval_high_majorAF), downSampling * nrow(eval_high_majorAF)) ,]
eval_high_majorAF$predicted.values.logit = predict(accuracyModel[["high_majorAF"]], newdata = as.data.frame(eval_high_majorAF))
eval_high_majorAF$predicted.values.prob = 1/(1+exp(-eval_high_majorAF$predicted.values.logit))

eval_data = rbind(eval_low_majorAF, eval_high_majorAF)
eval_data$majorAF_bin <- factor(eval_data$majorAF_bin, ordered = TRUE)

eval_data_binned <- eval_data %>% 
  group_by(coverage, majorAF_bin) %>%
  summarise(empirical.prob = mean(correct),
            empirical.sd = sd(correct),
            predicted.prob = mean(predicted.values.prob),
            predicted.sd = sqrt(sum(predicted.values.prob * (1 - predicted.values.prob)) / n()),
            n = n()) %>%
  pivot_longer(cols = -c("coverage", "majorAF_bin", "n"), names_to = c("method", ".value"), names_sep = "\\.") 


eval_data_binned_error <- eval_data_binned %>% 
  select(coverage, majorAF_bin, n, prob, method) %>%
  filter(n >= 50 & coverage <= 30) %>%
  pivot_wider(names_from = "method", values_from = "prob") %>%
  summarise(absolute_error = mean(abs(empirical - predicted)))


fwrite(eval_data_binned, file = opt$accuracyModelEvaluation)
fwrite(eval_data_binned_error, file = opt$accuracyModelEvaluationError)
#fwrite(eval_data, file = opt$allGenotypedSamplesAggUpdated)

