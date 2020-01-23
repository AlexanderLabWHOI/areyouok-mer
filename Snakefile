configfile: "config.yaml"

import io
import os
import yaml
import pandas as pd
import numpy as np
import pathlib
from snakemake.exceptions import print_exception, WorkflowError 

FASTQDUMP = config["fastqdumpdir"]
FASTQFILES = config["fastqdir"]
TRIMMEDFASTQFILES = config["trimmedfastqdir"]
TRIMMEDINTERMED = config["trimmedfastqdirintermed"]
OUTPUTDIR = config["outputsourmash"]
COMPAREDIR = config["outputcomparison"]

# Before `rule all`: run `download-fastqs.sh` with output into FASTQDUMP directory.

# how many unique accession numbers do we have from our fasterqdump?
fastqfiles = os.listdir(FASTQDUMP)
uniqueaccessions = set([(curr.split("_")[0]).split(".fastq")[0] for curr in fastqfiles])
print(uniqueaccessions)
kmers = [21,31,51]

include: "modules/concatinterleave"
include: "modules/trimup"
include: "modules/cutlowabundance"
include: "modules/sourmash"

rule all:
    input:
        # interleave the output files from the NCBI fasterqdump
        fastqs_interleaved = expand(os.path.join(FASTQFILES, "{accession}.fastq"), accession = uniqueaccessions),
        # next, we need to do error trimming. we use seqtk
        trimmed_intermed = expand(os.path.join(TRIMMEDINTERMED, "{accession}.fastq"), accession = uniqueaccessions),
        # next, cut low abundance kmers. we use khmer
        trimmed_fastqs = expand(os.path.join(TRIMMEDFASTQFILES, "{accession}.fastq"), accession = uniqueaccessions),
        # the next step is computing the sourmash signatures 
        signatures = expand(os.path.join(OUTPUTDIR, "{accession}.fastq.sig"), accession = uniqueaccessions),
        # the final step is to compute the sourmash comparisons at different k-mer lengths
        comparisons = expand(os.path.join(COMPAREDIR, "thaps_k{kmer}.cmp.csv"), kmer = kmers)
    
