
This pipeline merges paired end reads, deduplicates UMIs, profiles
motif content in syntehtic UTRs and counts them.

# Inputs

- Trimmed fastq files with UMIs on read1
Input reads must be the lone present in the cwd with suffix ".fastq.[1|2].gz"

- see pipeline.yml for more options/requirements.

# Outputs

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