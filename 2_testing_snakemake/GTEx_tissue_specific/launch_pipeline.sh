snakemake --configfile GTEx_config.yaml --cluster-sync "qsub -l mem_free={resources.memory} -l h_vmem={resources.memory} -l h_fsize=800G -e logs -o logs" --jobs 1000 --latency-wait 10 -k
