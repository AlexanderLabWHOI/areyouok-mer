# areyouok-mer
eukhashing project!

Directory for the Snakemake workflow for checking sourmash similarities of _Thaps_ experimental sequences and assembled genomes

## Steps

Branch 1: Raw Sequences (experiments pulled from online)
1. Interleave/concatenate
2. Do error trimming (coverage 5-10%) - good tool available on Conda = https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-017-1469-3
https://anaconda.org/bioconda/afterqc
3. Sourmash with k = 21, 31, 51, track abundance 

Branch 2: pre-assembled
1. Sourmash step 
2. Run on whole genome as well as just the coding region of the genome
3. Just work with cleaned assemblies, not combined assembles
