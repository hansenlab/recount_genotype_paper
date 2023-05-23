snakemake --cluster-sync "qsub -sync yes -l mem_free=50G -l h_vmem=50G -l h_fsize=700G -e ~/error -o ~/error" --jobs 100 --latency-wait 60
