#!/bin/bash
#SBATCH --partition=compute
#SBATCH --qos=unlim
#SBATCH -n 8

accessionNums="SRR4068200 SRR4068201 SRR4068202 SRR1202305 SRR1202307 SRR1202309 SRR1202311 SRR1202313 SRR1202315"

for f in $accessionNums; do fasterq-dump "$f" -O data/ncbi_thaps_fastqs; done
