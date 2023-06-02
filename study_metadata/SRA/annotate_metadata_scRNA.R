# Load required packages
library(recount3)
library(tidyverse)

# Load the spreadsheet containing the "study" column
project_data <- read.csv("/dcs04/hansen/data/recount_genotype/pipeline/AggregateFiles/all_SRA.csv")
study_info = data.frame(study = unique(project_data$study))

# Define a function to retrieve metadata for a given project and return the first element of the "pred.type" column
get_pred_type <- function(project_id) {
  cat(project_id, "\n")

  # Locate the metadata URL for the project
  url <- recount3::locate_url(project = project_id, project_home = 'data_sources/sra', type = 'metadata')
  
  # Retrieve the metadata files
  metadata_files <- recount3::file_retrieve(url)

  # Check that the last element of the URL contains the substring "sra.recount_pred"
  if(!grepl("sra\\.recount_pred", tail(metadata_files, 1))) {
    cat("The metadata file does not contain 'sra.recount_pred' in the name.\n")
    return(NA)
  }
  
  # Load the metadata file as a data frame
  metadata_df <- read.table(tail(metadata_files, 1), header = TRUE, sep = "\t")
  
  # Return the first element of the "pred.type" column
  return(metadata_df$pred.type[1])
}

# Apply the function to the "study" column of the project_data data frame and store the results in a new column "assay"

study_info$assay <- sapply(study_info$study, get_pred_type)

project_data$assay = NA
project_data$assay = study_info$assay[match(project_data$study, study_info$study)]

# Save the updated data frame to a new CSV file
write.csv(project_data, "/dcs04/hansen/data/recount_genotype/pipeline/AggregateFiles/all_SRA.scRNA_annotated.csv", row.names = FALSE)


