library(tidyverse)
library(cowplot)
library(RColorBrewer)

theme_set(theme_cowplot())

plot_file = "../ready_to_plot/Decision_boundary_simulation.rds"
if(!file.exists(plot_file)) {
	priors = list(c(.33, .33, .33),
				  c(.94, .04, .02),
				  c(.02, .04, .94))

	model_MS = readRDS("/dcs04/hansen/data/recount_genotype/pipeline/GTEx_Blended_Tissue_Training/train_genotyping_model/GTEx_Blended_Tissue_Training_genotyping_model.rds")


	lattice_all = data.frame()

	for(i in 1:length(priors)) {
		prior = priors[[i]]
		#construct a grid to find our decision boundary. 
		lattice = expand.grid(S = seq(0, 8, .1), M = seq(-10, 10, .1))
		lattice$mu_homRef <- predict(model_MS[[1]][[1]], data.frame(mu_S = lattice$S))
		lattice$sd_homRef <- predict(model_MS[[1]][[2]], data.frame(mu_S = lattice$S))
		lattice$mu_het <- predict(model_MS[[2]][[1]], data.frame(mu_S = lattice$S))
		lattice$sd_het <- predict(model_MS[[2]][[2]], data.frame(mu_S = lattice$S))
		lattice$mu_homAlt <- predict(model_MS[[3]][[1]], data.frame(mu_S = lattice$S))
		lattice$sd_homAlt <- predict(model_MS[[3]][[2]], data.frame(mu_S = lattice$S))

		#calculate log(prior * likelihood)
		lattice$dens_1 <- log(prior[1]) + dnorm(x = lattice$M, mean = lattice$mu_homRef, sd = lattice$sd_homRef, log = T)
		lattice$dens_2 <- log(prior[2]) + dnorm(x = lattice$M, mean = lattice$mu_het, sd = lattice$sd_het, log = T)
		lattice$dens_3 <- log(prior[3]) + dnorm(x = lattice$M, mean = lattice$mu_homAlt, sd = lattice$sd_homAlt, log = T)

		#calculate the MAP estimate from the posterior (3 classes), 
		#via log-sum-exp numerical stability computation: 
		#https://gregorygundersen.com/blog/2020/02/09/log-sum-exp/
		lattice$genotype = mapply(function(x, y, z) {
				log_prior_times_likelihood = c(x, y, z)
				log_sum_exp = log(sum(exp(log_prior_times_likelihood)))
				posterior = exp(log_prior_times_likelihood - log_sum_exp)
				return(which(posterior == max(posterior)))}, 
			lattice$dens_1, lattice$dens_2, lattice$dens_3)

		lattice$genotype = as.factor(lattice$genotype)
		levels(lattice$genotype) = c("AA", "AB", "BB")
		#lattice$priorName = paste0(round(prior, 2), collapse=", ")
		lattice$priorName = paste0("AA: ", prior[1], " AB: ", prior[2], " BB: ", prior[3])
		lattice_all = rbind(lattice_all, lattice)
	}
	saveRDS(lattice_all, file = plot_file)
}else {
	lattice_all = readRDS(plot_file)
}


pdf("../figures/Decision_boundary_simulation.pdf", width = 8.5, height = 3)
ggplot(lattice_all) + geom_point(aes(x = S, y = M, color = genotype), alpha = .2) + 
	scale_x_continuous(limits = c(1.161, 8)) +
	facet_wrap(~priorName) 
	#scale_color_brewer(palette = "Set2")
dev.off()


pdf("../figures/Decision_boundary_simulation_simple.pdf", width = 4, height = 3)
ggplot(lattice_all %>% filter(priorName == "AA: 0.33 AB: 0.33 BB: 0.33")) + 
	geom_point(aes(x = S, y = M, color = genotype)) + 
	scale_x_continuous(limits = c(1.161, 8)) +
	theme(legend.position = "none")
	#scale_color_brewer(palette = "Set2")
dev.off()