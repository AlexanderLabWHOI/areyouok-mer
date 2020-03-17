#!/bin/bash
#SBATCH --partition=compute

snakemake --jobs 10 --rerun-incomplete --use-conda --cluster-config cluster.yaml --cluster "sbatch --parsable --qos=unlim --partition={cluster.queue} --job-name=areyouokmer.{rule}.{wildcards} --mem={cluster.mem}gb --time={cluster.time} --ntasks={cluster.threads} --nodes={cluster.nodes}"


