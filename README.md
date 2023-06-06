# recount_genotype_paper

### Model training 

To train the genotyping model and accuracy model using our GTEx training set, run the snakemake pipeline in `training_snakemake/`. This was created using snakemake version 7.3.7.

### Model evaluation

We evaluated our genotyping model and accuracy model for the following datasets:

- GTEx testing set, reporting accuracy metrics for each tissue type: `testing_snakemake/GTEx_tissue_specific`

- GTEx testing set, reporting accuracy metrics for all tissues together: `testing_snakemake/GTEx_all_tissues`

- Geuvadis out-of-study test set, reporting accuracy metrics for the entire study: `testing_snakemake/geuvadis`

### Model application

- We applied our genotyping and accuracy model for all of SRA studies, which has no gold-standard genotypes to compute our model accuracy: `SRA_snakemake/`.

- We also applied our model for TCGA normal samples: `TCGA_snakemake/`.

### Manuscript figures

To generate our manuscript figures, see: `manuscript_figures/scripts`. 

