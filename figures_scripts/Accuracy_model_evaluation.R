library(tidyverse)
library(data.table)
library(rms) #restricted cubic splines
library(cowplot)
theme_set(theme_cowplot())

setwd("~/recount_genotype/redo_manuscript_figures")

plot_file = "ready_to_plot/Accuracy_model_evaluation.rds"

if(!file.exists(plot_file)) {
  GTEx_metadata = read.csv("/dcs07/hansen/data/recount_genotype/new_count_pipeline/new_count_pipeline/AggregateFiles/GTEx_testing.csv")
  #Geuvadis_metadata = read.csv("/dcs04/hansen/data/recount_genotype/pipeline/AggregateFiles/Geuvadis_testing_metadata.csv")
  #dcs07/hansen/data/recount_genotype/new_count_pipeline/new_count_pipeline/output_model_eval/
    # GTEx_metadata$accuracyModelEvaluation<-paste0("/dcs07/hansen/data/recount_genotype/new_count_pipeline/new_count_pipeline/output_model_eval/",GTEx_metadata$tissue,"_accuracy_model_evaluation.rds" )
    # GTEx_metadata$accuracyModelEvaluationError<-paste0("/dcs07/hansen/data/recount_genotype/new_count_pipeline/new_count_pipeline/output_model_eval/",GTEx_metadata$tissue,"_accuracy_model_evaluation_err.rds" )
    # fwrite(GTEx_metadata, "/dcs07/hansen/data/recount_genotype/new_count_pipeline/new_count_pipeline/AggregateFiles/GTEx_testing.csv")    
     all_GTEx = list()
  for(i in 1:nrow(GTEx_metadata)) {
    print(i)
    GTEx_tissue = fread(GTEx_metadata$accuracyModelEvaluation[i])
    GTEx_tissue$tissue = GTEx_metadata$tissue[i]
    all_GTEx[[i]] = GTEx_tissue
  }
  all_GTEx = do.call(rbind, all_GTEx)
  all_GTEx$study = "GTEx Test"
  
  geuvadis = fread("/dcs07/hansen/data/recount_genotype/new_count_pipeline/new_count_pipeline/output_model_eval/geuvadis_accuracyModelEvaluation.csv.gz")
  geuvadis$tissue = "Geuvadis LCL"
  geuvadis$study = "Geuvadis OOS"
  GTEx_vs_Geuvadis = rbind(all_GTEx, geuvadis)
  
  GTEx_vs_Geuvadis$study <- factor(GTEx_vs_Geuvadis$study, ordered = TRUE, levels = c("GTEx Test", "Geuvadis OOS"))
  
  GTEx_vs_Geuvadis$majorAF_bin = as.character(GTEx_vs_Geuvadis$majorAF_bin)
  GTEx_vs_Geuvadis$majorAF_bin = paste0("MAF: ", GTEx_vs_Geuvadis$majorAF_bin)
  GTEx_vs_Geuvadis$majorAF_bin[which(GTEx_vs_Geuvadis$majorAF_bin == "MAF: [0.5,0.632]")] = "MAF: [.5, .63]"
  GTEx_vs_Geuvadis$majorAF_bin[which(GTEx_vs_Geuvadis$majorAF_bin == "MAF: (0.632,0.749]")] = "MAF: (.63, .75]"
  GTEx_vs_Geuvadis$majorAF_bin[which(GTEx_vs_Geuvadis$majorAF_bin == "MAF: (0.749,0.839]")] = "MAF: (.75, .84]"
  GTEx_vs_Geuvadis$majorAF_bin[which(GTEx_vs_Geuvadis$majorAF_bin == "MAF: (0.839,0.908]")] = "MAF: (.84, .91]"
  GTEx_vs_Geuvadis$majorAF_bin[which(GTEx_vs_Geuvadis$majorAF_bin == "MAF: (0.908,0.95]")] = "MAF: (.91, .95]"
  GTEx_vs_Geuvadis$majorAF_bin[which(GTEx_vs_Geuvadis$majorAF_bin == "MAF: [0.95, 1]")] = "MAF: (.95, 1]"
  GTEx_vs_Geuvadis$majorAF_bin <- factor(GTEx_vs_Geuvadis$majorAF_bin, ordered = TRUE, levels = c("MAF: [.5, .63]",
                                                                                                  "MAF: (.63, .75]", "MAF: (.75, .84]", "MAF: (.84, .91]", "MAF: (.91, .95]", "MAF: (.95, 1]"))
  saveRDS(GTEx_vs_Geuvadis, plot_file)
}else {
  GTEx_vs_Geuvadis = readRDS(plot_file)
}

pdf("figure/Accuracy_model_evaluation_test_set.pdf", width = 8, height = 3.5)
ggplot(GTEx_vs_Geuvadis %>% filter(method == "empirical"), aes(x = coverage, y = prob)) + 
  geom_point(alpha = .85, size = 1, color = "darkslateblue") +
  geom_line(data = GTEx_vs_Geuvadis %>% filter(method == "predicted"), aes(x = coverage, y = prob, group = study), size = 1.1, color = "cornflowerblue") +
  facet_grid(study~majorAF_bin) +
  labs(x = "Coverage", y = "P(Correct Genotype)") +
  scale_x_continuous(limits = c(4, 30)) + 
  theme(legend.position = "none") +
  ylim(0.7,1)
dev.off()

tissue_comparison = rbind(GTEx_vs_Geuvadis %>% filter(tissue == "Cells_EBV-transformed_lymphocytes"),
                          GTEx_vs_Geuvadis %>% filter(tissue == "Thyroid"),
                          GTEx_vs_Geuvadis %>% filter(tissue == "Whole_Blood"),
                          GTEx_vs_Geuvadis %>% filter(tissue == "Geuvadis LCL"))
tissue_comparison$tissue = gsub("Cells_EBV-transformed_lymphocytes", "GTEx LCL", tissue_comparison$tissue)
tissue_comparison$tissue = gsub("Whole_Blood", "GTEx Whole Blood", tissue_comparison$tissue)
tissue_comparison$tissue = gsub("Thyroid", "GTEx Thyroid", tissue_comparison$tissue)
tissue_comparison$tissue <- factor(tissue_comparison$tissue, ordered = TRUE, levels = c("GTEx Thyroid", "GTEx Whole Blood", "GTEx LCL", "Geuvadis LCL"))

pdf("figure/Accuracy_model_evaluation_test_set_specificTissues.pdf", width = 8.5, height = 7)
ggplot(tissue_comparison %>% filter(method == "empirical"), aes(x = coverage, y = prob)) + 
  geom_point(alpha = .8, size = 1, color = "darkslateblue") +
  geom_line(data = tissue_comparison %>% filter(method == "predicted"), aes(x = coverage, y = prob), size = 1.1, color = "cornflowerblue") +
  labs(x = "Coverage", y = "P(Correct Genotype)") +
  facet_grid(tissue~majorAF_bin) +
  scale_x_continuous(limits = c(4, 30)) +
  theme(legend.position = "none")+
  ylim(0.7,1)
dev.off()


