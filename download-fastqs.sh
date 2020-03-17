#!/bin/bash
#SBATCH --partition=compute
#SBATCH --qos=unlim
#SBATCH -n 8

#accessionNums="SRR4068200 SRR4068201 SRR4068202 SRR1202305 SRR1202307 SRR1202309 SRR1202311 SRR1202313 SRR1202315"
accessionNums="SRR1202305 SRR1202307 SRR1202309 SRR1202311 SRR1202313 SRR1202315 SRR1202317 SRR1202319 SRR1202321 SRR1202323 SRR1202325 SRR1202327 SRR1202329 SRR1202331 SRR1202333 SRR1202335 SRR1202337 SRR1202339 SRR4068200 SRR4068201 SRR4068202 SRR4068203 SRR4068204 SRR4068205 SRR4068206 SRR4068207 SRR4068208"

accessionNums_242525="SRR6296274 SRR6296278 SRR6296279 SRR6296282 SRR6296286 SRR6296288 SRR6296294 SRR6296295 SRR6296296 SRR6296298 SRR6296276 SRR6296277 SRR6296280 SRR6296283 SRR6296284 SRR6296289 SRR6296290 SRR6296293 SRR6296297 SRR6296273 SRR6296275 SRR6296281 SRR6296285 SRR6296287 SRR6296291 SRR6296292 SRR6296299"

for f in $accessionNums; do fasterq-dump "$f" --split-files -O data/ncbi_thaps_fastqs; done
for f in $accessionNums; do awk -F"/| " '{ if (NR%2==1) { print $1 "/1 " $4 " " $5} else { print } }'  "data/ncbi_thaps_fastqs/$f"_1.fastq > "data/ncbi_thaps_fastqs/$f"_3.fastq; done
for f in $accessionNums; do awk -F"/| " '{ if (NR%2==1) { print $1 "/2 " $4 " " $5} else { print } }'  "data/ncbi_thaps_fastqs/$f"_2.fastq > "data/ncbi_thaps_fastqs/$f"_4.fastq; done

#for f in $accessionNums_242525; do sed -i -r 's/(^[\@\+]SRR\S+)/\1\/1/' "data/ncbi_thaps_fastqs/$f".fastq; done
#for f in $accessionNums_242525; do sed -i -r 's/(^[\@\+]SRR\S+)/\1\/2/' "data/ncbi_thaps_fastqs/$f".fastq; done
#for f in $accessionNums; do sed -i -r 's/(^[\@\+]SRR\S+)/\1\/1/' "data/ncbi_thaps_fastqs/$f"_1.fastq; done
#for f in $accessionNums; do sed -i -r 's/(^[\@\+]SRR\S+)/\1\/2/' "data/ncbi_thaps_fastqs/$f"_1.fastq; done
#for f in $accessionNums; do sed -i -r 's/(^[\@\+]SRR\S+)/\1\/1/' "data/ncbi_thaps_fastqs/$f"_2.fastq; done
#for f in $accessionNums; do sed -i -r 's/(^[\@\+]SRR\S+)/\1\/2/' "data/ncbi_thaps_fastqs/$f"_2.fastq; done

#SRA_TOOLKIT/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files file.sra
