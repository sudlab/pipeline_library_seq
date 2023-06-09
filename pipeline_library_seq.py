"""===========================
pipeline_library_seq.py
===========================

Overview
========

This pipeline merges paired end reads, deduplicates UMIs, profiles
motif content in syntehtic UTRs and counts them.

files :file:``pipeline.yml` and :file:`conf.py`.

Usage
=====

See :ref:`PipelineSettingUp` and :ref:`PipelineRunning` on general
information how to use CGAT pipelines.

Configuration
-------------

The pipeline requires a configured :file:`pipeline.yml` file.
CGATReport report requires a :file:`conf.py` and optionally a
:file:`cgatreport.ini` file (see :ref:`PipelineReporting`).

Default configuration files can be generated by executing:

   python <srcdir>/pipeline_library_seq.py config

Input files
-----------

Trimmed fastq files with UMIs on read1)

Input reads must be the lone present in the cwd with suffix ".fastq.[1|2].gz"

Requirements
------------

On top of the default CGAT setup, the pipeline requires the following
software to be in the path:
    - NGmerge module
    - UMI tools
Required R modules:
   - tidyverse

Pipeline output
===============

- number total reads before merge [PREFIX].nreads
- fastq with extracted UMIs in extractedUMIs.dir
- merged paired reads in mergedReads.dir
  see MGmerge docs
- motifContent.dir
    - [PREFIX]_failed.tsv
    List of reads failing UTR profiling, contains read name and log of fail
    Failed sequencing:
        1. The two RE cutting sites are back to back (directly connected)
        2. The RE cut sites are present (likely empty vector).
        3. The sequence upstream of the insert was not detected.
        4. The sequence downstream of the insert was not detected.
        6. Ambiguity in linker of motif sequence during profiling.
    - [PREFIX]_sUTR_counts.tsv
    - [PREFIX]_stats.txt

Code
====

"""

import sys
import os
import sqlite3
import csv
import re
import glob

from cgatcore import pipeline as P
import cgat.GTF as GTF
import cgatcore.iotools as IOTools
from ruffus import *

#Load config file options
PARAMS = P.get_parameters(
    ["%s/pipeline.yml" % os.path.splitext(__file__)[0],
     "../pipeline.yml",
     "pipeline.yml"])
########################


@transform("*.fastq.1.gz",
           regex("(.+).fastq.1.gz"),
           r"\1.nreads")
def countReads(infile, outfile):
    '''Count total number of reads'''
    statement = '''
    expr $(zcat %(infile)s | wc -l) / 4 > %(outfile)s
    '''
    P.run(statement,
          job_memory="2G")

@transform("*.fastq.1.gz",
           regex("(.+).fastq.1.gz"),
           add_inputs(r"\1.fastq.2.gz"),
           r"extractedUMIs.dir/\1")
def extractUMIs(infiles, outfile):
    '''umi extract'''
    read1, read2 = infiles
    out2 = P.snip(outfile, ".1.gz")+".2.gz"
    logfile = outfile+".log"
    statement = '''
    umi_tools extract -I %(read1)s --read2-in=%(read2)s
                      -S %(outfile)s --read2-out=%(out2)s
                       --bc-pattern=NNNNNNNN 
                       --log=%(logfile)s
    '''
    P.run(statement,
          job_memory="8G")


@transform(extractUMIs,
           regex("(.+).fastq.1.gz"),
           add_inputs(r"\1.fastq.2.gz"),
           r"mergedReads.dir/\1.fastq.gz")
def mergePairedReads(infiles, outfile):
    '''Merge paired end reads'''
    read1, read2 = infiles
    logfile1 = P.snip(outfile, ".fastq.gz")+".log"
    logfile2 = P.snip(outfile, ".fastq.gz")+".align.log"
    fail_reads = P.snip(outfile, ".fastq.gz")+"_fail.fastq"
    n_threads = PARAMS["number_of_threads"]
    n_align = PARAMS["min_alignment"]
    mismatch = PARAMS["fraction_mismatch_align"]
    statement = '''
    module load bio/NGmerge &&
    NGmerge -1 %(read1)s -2 %(read2)s -m %(n_align)s
    -o %(outfile)s
    -f %(fail_reads)s
    -l %(logfile1)s 
    -j %(logfile2)s
    -z 
    -n 2
    -p %(mismatch)s
    '''
    P.run(statement,
          job_memory="8G")

@transform(mergePairedReads,
           regex("(.+)_merged.fastq.gz"),
           r"motifContent.dir/\1_SUTRs_counts.tsv")
def motifContent(infile, outfile):
    '''Infer motif content in sequenced sUTRs'''
    pre = P.snip(outfile, "_SUTRs_motifs_content.tsv")
    logfile = P.snip(outfile, ".tsv")+".log"
    link = PARAMS["linker"]
    up = PARAMS["upstream"]
    dw = PARAMS["downstream"]
    RE1 = PARAMS["5primeRE"]
    RE2 = PARAMS["3primeRE"]
    motifs = PARAMS["motifs_list"]
    script_path = os.path.join((os.path.dirname(__file__)),
                                "scripts",
                                "motifContentFromFastQ.py")
    statement = '''
    python %(script_path)s 
           -i %(infile)s
           -o ./ 
           -op %(pre)s
           -l %(link)s
           -1 %(RE1)s
           -2 %(RE2)s
           -u %(up)s
           -d %(dw)s
           -m %(motifs)s
           -L %(logfile2)s
    '''
    P.run(statement,
          job_memory="8G")



@follows(motifContent)
def full():
    '''Later alligator'''
    pass

P.main()
