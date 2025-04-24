library(tidyverse)
library(data.table)
recount3_metadata<-read_tsv("Recount3_metadata.tsv") #this metadata was downloaded from the Recount3 website

#Only select bulk samples
met<-recount3_metadata[,c(colnames(recount3_metadata)[1:6],"run_published","library_layout","seq_type")]
met<-met[met$seq_type=="bulk",]
#Count the number of runs in each experiment
pp<-met %>% group_by(experiment_acc) %>% summarise(n=n())
 
#Make the metadata to rerun the genotyping pipeline
#Here we want to only rerun the pipeline for experiments with more than 1 run
#Make sure that all the runs were submitted on the same date! Remove any experiment that have runs with variable submission date bc we cannot trust those runs

id<-pp[pp$n>1,]

re_run<-met[met$experiment_acc %in% id$experiment_acc,]

unique_date<-re_run %>% mutate(date=date(run_published)) %>% group_by(experiment_acc) %>% summarise(num_date=length(unique(date)))

id2<-unique_date[unique_date$num_date==1,]
exp_metadata<-met[met$experiment_acc %in% unique_date$experiment_acc,]

#if runs do not have the same published date, do not trust those samples
exp_metadata$warning<-NA
exp_metadata$warning[!(exp_metadata$experiment_acc %in% id2$experiment_acc)]<-"Possible bad sample:Runs published in multiple days"
dim(exp_metadata) #55684
fwrite(exp_metadata, file="multi_run_exp.csv")


#add the alt and total bigwig path:
sra<-read.csv("all_SRA.csv")
sra_exp<-sra[sra$sample_id %in% exp_metadata$external_id,]

sra_exp$experiment_acc<-exp_metadata$experiment_acc[match(sra_exp$sample_id,exp_metadata$external_id)]
sra_exp$warning<-exp_metadata$warning[match(sra_exp$sample_id,exp_metadata$external_id)]

fwrite(sra_exp, file="experiment_metadata.csv")
