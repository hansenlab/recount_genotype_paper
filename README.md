# Genotype calling from expression data in Recount3

Afrooz Razi, Christopher C. Lo, Siruo Wang, Jeffrey T. Leek, Kasper D. Hansen

see https://github.com/raziafrooz/RecountGenotyper for R package.

### Overview:
In this paper, we used raw reference and alternative read counts from Recount3 expression data to predict the genotype at biallelic SNPs. Using our machine learning model we were able to predict the genotype for all the recount3 samples and used that to generate a PCA plot containing the underlying population structure of the samples. 

## Candidate list of biallelic SNPs from GTEx

GTEx whole genome DNA-seq was used to select biallelic SNPs using "bcftools view -m2 -M2 -v snps".

Only SNPs in the protein coding regions were kept using ensemble V85
  
## Model training 

We trained the model using samples in various tissues in GTEx
To train the genotyping model and accuracy model using our GTEx training set, run the snakemake pipeline in `training_snakemake/`. This was created using snakemake version 7.3.7.

## Model evaluation

We evaluated our genotyping model and accuracy model for the following datasets:

- GTEx testing set, reporting accuracy metrics for each tissue type: `testing_snakemake/GTEx_tissue_specific`

- GTEx testing set, reporting accuracy metrics for all tissues together: `testing_snakemake/GTEx_all_tissues`

- Geuvadis out-of-study test set, reporting accuracy metrics for the entire study: `testing_snakemake/geuvadis`

## Model application

- We applied our genotyping and accuracy model for all of SRA studies, which has no gold-standard genotypes to compute our model accuracy: `SRA_snakemake/`.

- We also applied our model for TCGA normal samples: `TCGA_snakemake/`.

## Manuscript figures

To generate our manuscript figures, see: `figures_scripts`. 

