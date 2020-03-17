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
TRINITYFASTADIR = config["trinityfastadir"]
TRINITYOUTDIR = config["trinityoutputdir"]
OUTPUTDIR = config["outputsourmash"]
COMPAREDIR = config["outputcomparison"]
TRANSDECODERDIR = config["transdecoder"]
BACKTRANSLATEDDIR = config["backtrans"]
BACKTRANSLATEDSIG = config["backtranssig"]

# Before `rule all`: run `download-fastqs.sh` with output into FASTQDUMP directory.

# how many unique accession numbers do we have from our fasterqdump?
fastqfiles = os.listdir(FASTQDUMP)
uniqueaccessions = set([(curr.split("_")[0]).split(".fastq")[0] for curr in fastqfiles])
uniqueaccessionsplusref = uniqueaccessions.union(set(["thapsreference","thapsreferencecoding","pseudreference","pseudreferencecoding","fragreference","fragreferencecoding"]))

toassemble = [curr for curr in fastqfiles if (("_1.fastq" in curr) | ("_2.fastq" in curr))]
toassemble = set([(curr.split("_")[0]).split(".fastq")[0] for curr in toassemble])
print(toassemble)
print("assemblers above")
uniqueaccessionsplusrefandtrinity = uniqueaccessionsplusref.union(toassemble)

trinityoutputs = os.listdir(TRINITYFASTADIR)
print(uniqueaccessions)
kmers = [21,31,51]

#trinity_assemblies = expand(os.path.join(TRINITYOUTDIR, "trinity_results_sample_{accession}", "Trinity.fasta"), accession = toassemble),

include: "modules/concatinterleave"
include: "modules/trimup"
include: "modules/cutlowabundance"
include: "modules/sourmash"
include: "modules/trinity"
include: "modules/transdecoder"
include: "modules/convertgff"

rule all:
    input:
        # interleave the output files from the NCBI fasterqdump
        #fastqs_interleaved = expand(os.path.join(FASTQFILES, "{accession}.fastq"), accession = uniqueaccessions),
        # next, we need to do error trimming. we use seqtk
        #trimmed_intermed = expand(os.path.join(TRIMMEDINTERMED, "{accession}.fastq"), accession = uniqueaccessionsplusref),
        # next, cut low abundance kmers. we use khmer
        #trimmed_fastqs = expand(os.path.join(TRIMMEDFASTQFILES, "{accession}.fastq"), accession = uniqueaccessionsplusref),
        # next, do a Trinity assembly for the paired end reads
        #trinity_assemblies = expand(os.path.join(TRINITYFASTADIR, "{accession}.fasta"), accession = toassemble),
        # the next step is computing the sourmash signatures for all except the Trinity assemblies
        signatures = expand(os.path.join(OUTPUTDIR, "{accession}.fastq.sig"), accession = uniqueaccessionsplusref),
        # now compute the signatures for the Trinity assemblies
        signaturestrinity = expand(os.path.join(OUTPUTDIR, "{accession}.fasta.sig"), accession = toassemble),
        # now do the transdecoder files!
        #transdecoder = expand(os.path.join(TRANSDECODERDIR, "{accession}.fasta.transdecoder.pep"), accession = toassemble),
        # now convert the transdecoder files to AA!
        backtrans = expand(os.path.join(BACKTRANSLATEDDIR, "{accession}.fasta"), accession = toassemble),
        # now do sourmash scores!
        signaturesbacktrans = expand(os.path.join(BACKTRANSLATEDSIG, "{accession}.fasta.sig"), accession = toassemble),
        # now do sourmash scores!
        signaturesprotein = expand(os.path.join(TRANSDECODERDIR, "{accession}.fasta.transdecoder.pep.sig"), accession = toassemble),
        # the final step is to compute the sourmash comparisons at different k-mer lengths
        comparisons = expand(os.path.join(COMPAREDIR, "thaps_k{kmer}.cmp.csv"), kmer = kmers)
    
