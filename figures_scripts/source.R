conditional1<-function(geno_df){ 
  geno_df %>% 
    filter(AF <= .5) %>%
    mutate(
      major_accuracy = ifelse(true_genotype == "1" & pred_genotype == "1", 1, 0) +
        ifelse(true_genotype == "1" & pred_genotype == "2", .5, 0) +
        ifelse(true_genotype == "1" & pred_genotype == "3", 0, 0) +
        ifelse(true_genotype == "2" & pred_genotype == "1", 1, 0) +
        ifelse(true_genotype == "2" & pred_genotype == "2", 1, 0) +
        ifelse(true_genotype == "2" & pred_genotype == "3", 0, 0) +
        ifelse(true_genotype == "3" & pred_genotype == "1", NA, 0) +
        ifelse(true_genotype == "3" & pred_genotype == "2", NA, 0) +
        ifelse(true_genotype == "3" & pred_genotype == "3", NA, 0),
      minor_accuracy = ifelse(true_genotype == "1" & pred_genotype == "1", NA, 0) +
        ifelse(true_genotype == "1" & pred_genotype == "2", NA, 0) +
        ifelse(true_genotype == "1" & pred_genotype == "3", NA, 0) +
        ifelse(true_genotype == "2" & pred_genotype == "1", 0, 0) +
        ifelse(true_genotype == "2" & pred_genotype == "2", 1, 0) +
        ifelse(true_genotype == "2" & pred_genotype == "3", 1, 0) +
        ifelse(true_genotype == "3" & pred_genotype == "1", 0, 0) +
        ifelse(true_genotype == "3" & pred_genotype == "2", .5, 0) +
        ifelse(true_genotype == "3" & pred_genotype == "3", 1, 0),
      overall_accuracy = ifelse(true_genotype == "1" & pred_genotype == "1", 1, 0) +
        ifelse(true_genotype == "1" & pred_genotype == "2", .5, 0) +
        ifelse(true_genotype == "1" & pred_genotype == "3", 0, 0) +
        ifelse(true_genotype == "2" & pred_genotype == "1", .5, 0) +
        ifelse(true_genotype == "2" & pred_genotype == "2", 1, 0) +
        ifelse(true_genotype == "2" & pred_genotype == "3", .5, 0) +
        ifelse(true_genotype == "3" & pred_genotype == "1", 0, 0) +
        ifelse(true_genotype == "3" & pred_genotype == "2", .5, 0) +
        ifelse(true_genotype == "3" & pred_genotype == "3", 1, 0)) %>%
    summarise(sum_overall_accuracy = sum(overall_accuracy, na.rm=T),
              sum_major_accuracy = sum(major_accuracy, na.rm=T),
              sum_minor_accuracy = sum(minor_accuracy, na.rm=T),
              n = sum(!is.na(true_genotype)),
              n_major = sum(!is.na(major_accuracy)),
              n_minor = sum(!is.na(minor_accuracy))
    ) %>%
    collect()
}

cat("Computing minor accuracy.\n")


conditional2<-function(geno_df){ 
  geno_df %>% 
    filter(AF > .5) %>%
    mutate(
      minor_accuracy = ifelse(true_genotype == "1" & pred_genotype == "1", 1, 0) +
        ifelse(true_genotype == "1" & pred_genotype == "2", .5, 0) +
        ifelse(true_genotype == "1" & pred_genotype == "3", 0, 0) +
        ifelse(true_genotype == "2" & pred_genotype == "1", 1, 0) +
        ifelse(true_genotype == "2" & pred_genotype == "2", 1, 0) +
        ifelse(true_genotype == "2" & pred_genotype == "3", 0, 0) +
        ifelse(true_genotype == "3" & pred_genotype == "1", NA, 0) +
        ifelse(true_genotype == "3" & pred_genotype == "2", NA, 0) +
        ifelse(true_genotype == "3" & pred_genotype == "3", NA, 0),
      major_accuracy = ifelse(true_genotype == "1" & pred_genotype == "1", NA, 0) +
        ifelse(true_genotype == "1" & pred_genotype == "2", NA, 0) +
        ifelse(true_genotype == "1" & pred_genotype == "3", NA, 0) +
        ifelse(true_genotype == "2" & pred_genotype == "1", 0, 0) +
        ifelse(true_genotype == "2" & pred_genotype == "2", 1, 0) +
        ifelse(true_genotype == "2" & pred_genotype == "3", 1, 0) +
        ifelse(true_genotype == "3" & pred_genotype == "1", 0, 0) +
        ifelse(true_genotype == "3" & pred_genotype == "2", .5, 0) +
        ifelse(true_genotype == "3" & pred_genotype == "3", 1, 0),
      overall_accuracy = ifelse(true_genotype == "1" & pred_genotype == "1", 1, 0) +
        ifelse(true_genotype == "1" & pred_genotype == "2", .5, 0) +
        ifelse(true_genotype == "1" & pred_genotype == "3", 0, 0) +
        ifelse(true_genotype == "2" & pred_genotype == "1", .5, 0) +
        ifelse(true_genotype == "2" & pred_genotype == "2", 1, 0) +
        ifelse(true_genotype == "2" & pred_genotype == "3", .5, 0) +
        ifelse(true_genotype == "3" & pred_genotype == "1", 0, 0) +
        ifelse(true_genotype == "3" & pred_genotype == "2", .5, 0) +
        ifelse(true_genotype == "3" & pred_genotype == "3", 1, 0)) %>%
    summarise(sum_overall_accuracy = sum(overall_accuracy,na.rm=T),
              sum_major_accuracy = sum(major_accuracy, na.rm=T),
              sum_minor_accuracy = sum(minor_accuracy, na.rm=T),
              n = sum(!is.na(true_genotype)),
              n_major = sum(!is.na(major_accuracy)),
              n_minor = sum(!is.na(minor_accuracy))
    ) %>%
    collect()
}

