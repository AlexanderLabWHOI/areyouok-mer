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
REFERENCEDIR = config["referencedir"]
BLASTOUTDIR = config["blastoutdir"]


# Before `rule all`: run `download-fastqs.sh` with output into FASTQDUMP directory. 
# Also, run `organize_commands.ipynb` to create the reference file, and specify this file in the config.
indexfile = pd.read_csv(config["indexfile"], sep = "\t", index_col = 2)
referencefiles = (indexfile.loc["refgen"])
refname = list(referencefiles["FileName"])
refext = list(referencefiles["FileExtension"])

rawfiles = (indexfile.loc[["onesid","twosid"]])
mmetspfiles = (indexfile.loc[["reftrans"]])

kmers = [21,33,51]
kmersprot = [12,15,18,21,33,51]

#trinity_assemblies = expand(os.path.join(TRINITYOUTDIR, "trinity_results_sample_{accession}", "Trinity.fasta"), accession = toassemble),

include: "modules/concatinterleave"
include: "modules/trimup"
#include: "modules/cutlowabundance"
include: "modules/sourmash"
include: "modules/trinity"
include: "modules/transdecoder"
include: "modules/convertgff"
include: "modules/blastassemblies"
    
print(list(rawfiles["FileExtension"]))

rule all:
    input:
        # interleave the output files from the NCBI fasterqdump
        fastqs_interleaved = expand(os.path.join(FASTQFILES, "{accession}.{ext}"), zip, accession = list(rawfiles["FileName"]), ext = list(rawfiles["FileExtension"])),
        # next, we need to do error trimming. we use seqtk
        trimmed_intermed = expand(os.path.join(TRIMMEDINTERMED, "{accession}.{ext}"), zip, accession = list(rawfiles["FileName"]), ext = list(rawfiles["FileExtension"])),
        # next, cut low abundance kmers. we use khmer
        #trimmed_fastqs = expand(os.path.join(TRIMMEDFASTQFILES, "{accession}.fastq"), accession = uniqueaccessionsplusref),
        # next, do a Trinity assembly for the paired end reads (need to edit Trinity module)
        #trinity_assemblies = expand(os.path.join(TRINITYFASTADIR, "{accession}_twoside.fasta"), accession = toassemble),
        # next, do a Trinity assembly for all the reads
        trinity_assemblies = expand(os.path.join(TRINITYFASTADIR, "{accession}.fasta"), accession = list(rawfiles["FileName"])),
        # the next step is computing the sourmash signatures for all except the Trinity assemblies
        signatures = expand(os.path.join(OUTPUTDIR, "{accession}.{ext}.sig"), zip, accession = list(rawfiles["FileName"]), ext = list(rawfiles["FileExtension"])),
        # now compute the signatures for the Trinity assemblies
        signaturestrinity = expand(os.path.join(OUTPUTDIR, "{accession}.fasta.sig"), accession = list(rawfiles["FileName"])),
        # now compute the signatures for the MMETSP files
        signaturesmmetsp = expand(os.path.join(MMETSPDIR, "{accession}.fasta.sig"), accession = list(mmetspfiles["FileName"])),
        # now compute the signatures for the reference files
        signaturesref = expand(os.path.join(OUTPUTDIR, "{accession}.{ext}.sig"), accession = refname, ext = refext),
        # now do the transdecoder files!
        transdecoder = expand(os.path.join(TRANSDECODERDIR, "{accession}.fasta.transdecoder.pep"), accession = list(rawfiles["FileName"])),
        # now convert the transdecoder files to AA!
        backtrans = expand(os.path.join(BACKTRANSLATEDDIR, "{accession}.fasta"), accession = list(rawfiles["FileName"])),
        # now do sourmash scores!
        signaturesbacktrans = expand(os.path.join(BACKTRANSLATEDSIG, "{accession}.fasta.sig"), accession = list(rawfiles["FileName"])),
        # now do sourmash scores!
        signaturesprotein = expand(os.path.join(TRANSDECODERDIR, "{accession}.fasta.transdecoder.pep.sig"), accession = list(rawfiles["FileName"])),
        # the final step is to compute the sourmash comparisons at different k-mer lengths
        comparisons = expand(os.path.join(COMPAREDIR, OUTPUTNAME + "_k{kmer}.cmp.csv"), kmer = kmers),
        # now do protein comparisons
        comparisonsprotein = expand(os.path.join(COMPAREDIR, OUTPUTNAME + "_k{kmer}_protein.cmp.csv"), kmer = kmersprot),
        # now do a blast of the Trinity assemblies vs the raw for _1 files
        blastcomparisons = expand(os.path.join(BLASTOUTDIR, "{sample}.{ext}.txt"), zip, sample = list(rawfiles["FileName"]), ext = list(rawfiles["FileExtension"]))
    
