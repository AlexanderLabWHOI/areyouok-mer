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
STARTDIR = config["startingdir"] # do we start with assemblies or with something else? 
STARTEXTENSION = config["startingextension"] # extension for our starting point files
OUTPUTDIR = config["outputsourmash"]
COMPAREDIR = config["outputcomparison"]
TRANSDECODERDIR = config["transdecoder"]
BACKTRANSLATEDDIR = config["backtrans"]
BACKTRANSLATEDSIG = config["backtranssig"]
OUTPUTNAME = config["outputname"]
MMETSPDIR = config["mmetspdir"]


# Before `rule all`: run `download-fastqs.sh` with output into FASTQDUMP directory.

# how many unique accession numbers do we have from our fasterqdump?
fastqfiles = os.listdir(STARTDIR)
uniqueaccessions = set([(curr.split("_")[0]).split(STARTEXTENSION)[0] for curr in fastqfiles if STARTEXTENSION in curr])
uniqueaccessionsplusref = uniqueaccessions.union(set(config["referencenames"]))

# get all the files from the MMETSP dir
fastqfiles = os.listdir(MMETSPDIR)
mmetspnames = set([(curr.split("_")[0]).split(".fasta")[0] for curr in fastqfiles if ".fasta" in curr])
mmetspnames = mmetspnames.union(set(config["referencenames"]))
print("MMETSP names are: ")
print(mmetspnames)

# get rid of this later - more hard-coded for NCBI fastq dump
#fastqfiles = os.listdir(FASTQDUMP)
#uniqueaccessions = set([(curr.split("_")[0]).split(".fastq")[0] for curr in fastqfiles])
#uniqueaccessionsplusref = #uniqueaccessions.union(set(["thapsreference","thapsreferencecoding","pseudreference","pseudreferencecoding","fragreference","fragreferencecoding"]))

# currently set to only assemble paired end reads
toassemble = [curr for curr in fastqfiles if (("_1.fastq" in curr) | ("_2.fastq" in curr))]
toassemble = set([(curr.split("_")[0]).split(".fastq")[0] for curr in toassemble])
print(toassemble)
print("assemblers above")
uniqueaccessionsplusrefandtrinity = uniqueaccessionsplusref.union(toassemble)

trinityoutputs = os.listdir(TRINITYFASTADIR)
print(uniqueaccessions)
kmers = [21,33,51]
kmersprot = [12,15,18,21,33,51]

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
        #signatures = expand(os.path.join(OUTPUTDIR, "{accession}.fastq.sig"), accession = uniqueaccessionsplusref),
        # now compute the signatures for the Trinity assemblies
        # signaturestrinity = expand(os.path.join(OUTPUTDIR, "{accession}.fasta.sig"), accession = trinityoutputs),
        # now compute the signatures for the MMETSP/reference files
        signaturesmmetsp = expand(os.path.join(MMETSPDIR, "{accession}.fasta.sig"), accession = mmetspnames),
        # now do the transdecoder files!
        transdecoder = expand(os.path.join(TRANSDECODERDIR, "{accession}.fasta.transdecoder.pep"), accession = mmetspnames),
        # now convert the transdecoder files to AA!
        backtrans = expand(os.path.join(BACKTRANSLATEDDIR, "{accession}.fasta"), accession = mmetspnames),
        # now do sourmash scores!
        signaturesbacktrans = expand(os.path.join(BACKTRANSLATEDSIG, "{accession}.fasta.sig"), accession = mmetspnames),
        # now do sourmash scores!
        signaturesprotein = expand(os.path.join(TRANSDECODERDIR, "{accession}.fasta.transdecoder.pep.sig"), accession = mmetspnames),
        # the final step is to compute the sourmash comparisons at different k-mer lengths
        comparisons = expand(os.path.join(COMPAREDIR, OUTPUTNAME + "_k{kmer}.cmp.csv"), kmer = kmers),
        # now do protein comparisons
        comparisonsprotein = expand(os.path.join(COMPAREDIR, OUTPUTNAME + "_k{kmer}_protein.cmp.csv"), kmer = kmersprot)
    
