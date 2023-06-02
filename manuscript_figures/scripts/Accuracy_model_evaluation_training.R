library(tidyverse)
library(data.table)
library(rms) #restricted cubic splines
library(cowplot)
theme_set(theme_cowplot())


plot_file = "../ready_to_plot/Accuracy_model_evaluation_training.rds"

if(!file.exists(plot_file)) {
	GTEx_training_metadata = read.csv("/dcs04/hansen/data/recount_genotype/pipeline/AggregateFiles/GTEx_Blended_Tissue_Training_metadata.csv")
	model = readRDS(GTEx_training_metadata$accuracy_model[1])
	training_lowAF = model$low_majorAF$data
	training_lowAF$predicted.values.logit = predict(model$low_majorAF, newdata = as.data.frame(training_lowAF))
	training_lowAF$predicted.values.prob = 1/(1+exp(-training_lowAF$predicted.values.logit))

	training_highAF = model$high_majorAF$data
	training_highAF$predicted.values.logit = predict(model$high_majorAF, newdata = as.data.frame(training_highAF))
	training_highAF$predicted.values.prob = 1/(1+exp(-training_highAF$predicted.values.logit))

	training = rbind(training_highAF, training_lowAF)

	train_metadata = read.csv("../../training_snakemake/GTEx_Blended_Tissue_Training_metadata.csv")
	training$tissue = train_metadata$tissue[match(training$sample_id_rep, train_metadata$sample_id_rep)]

	training_bins = list()
	for(i in 1:length(unique(training$tissue))) {
	  t = unique(training$tissue)[i]
	  cat(t, "\n")
	  training_bin_i <- training %>% filter(tissue == t) %>%
		    group_by(coverage, majorAF_bin) %>%
		    summarise(empirical.prob = mean(correct),
		              empirical.sd = sd(correct),
		              predicted.prob = mean(predicted.values.prob),
		              predicted.sd = sqrt(sum(predicted.values.prob * (1 - predicted.values.prob)) / n()), #n()^2??
		              n = n()) %>%
		    pivot_longer(cols = -c("coverage", "majorAF_bin", "n"), names_to = c("method", ".value"), names_sep = "\\.") %>%
		    filter(n >= 50)
	  training_bin_i$tissue = t
	  training_bins[[i]] = training_bin_i
	}

	training_bins = do.call(rbind, training_bins)
	saveRDS(training_bins, plot_file)
} else {
	training_bins = readRDS(plot_file)
}


pdf("../figures/Accuracy_model_evaluation_train_set.pdf", width = 7, height = 5)

ggplot(training_bins %>% filter(method == "empirical"), aes(x = coverage, y = prob)) + 
  geom_line(aes(group = tissue), alpha = .2, linewidth = .85) +
  facet_wrap(~majorAF_bin) +
  labs(x = "Coverage", y = "P(Correct Genotype)") +
  scale_x_continuous(limits = c(4, 30)) + 
  theme(legend.position = "none")

dev.off()
