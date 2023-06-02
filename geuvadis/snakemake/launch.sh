snakemake --cluster-sync "qsub -sync yes -l mem_free=50G -l h_vmem=50G -l h_fsize=100G -e /users/clo/recount_genotype/geuvadis -o /users/clo/recount_genotype/geuvadis" --jobs 100 --latency-wait 60
