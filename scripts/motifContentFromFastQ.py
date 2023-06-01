'''
motifContentFromSequence.py - template for cgat scripts
====================================================

:Author:
:Tags: Python

Purpose
-------

From a fastq file containing sequencing results for synthetic UTRs library building,
determine motifs content in the sUTRs, deduplicate reads and counts reads per sUTR.

Outputs:
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


Allowing 1 mismatch in linker/motif sequence.

Usage
-----

.. Example use case

Example::

   python motifContentFromSequence.py -i input.fastq(.gz) -o /path/to/out/directory

Type::

   python motifContentFromSequence.py --help

for command line help.

Command line options
--------------------

"-i", "--input-fasta": Fastq file of sequencing, with IDs
"-l", "--linker-sequence":Table of lowstab motifs to remove from highstab
"-1", "--RE1-site": Sequence of Restriction Enzyme 1 cut site
"-2", "--RE2-site": Sequence of Restriction Enzyme 1 cut site
"-o", "--output-directory: path to directory where outputs are written
"-p", "--output-prefix": Output prefix
"-u", "--upstream-sequence: Excepected sequence upstream of insert
"-d", "--downstream-sequence: Excepected sequence downstream of insert
"-m", "--list-motifs": List of motifs expected in insert



'''
import sys
import cgatcore.experiment as E
import cgatcore.iotools as iotools
from Bio import SeqIO
import pandas as pd
import regex
from umi_tools import UMIClusterer
import os

def main(argv=None):
    """script main.
    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $1.0$",
                            usage=globals()["__doc__"])

    parser.add_option("-i", "--input-fasta", dest="fastq", type=str,
                        help="Fasta files of sequencing, with IDs")
    parser.add_option("-l", "--linker-sequence", dest="linker", type=str,
                        help="Table of lowstab motifs to remove from highstab")
    parser.add_option("-1", "--RE1-site", dest="RE1site", type=str,
                        help="Sequence of Restriction Enzyme 1 cut site")
    parser.add_option("-2", "--RE2-site", dest="RE2site", type=str,
                        help="Sequence of Restriction Enzyme 1 cut site")
    parser.add_option("-o", "--output-directory", dest="out_dir", type=str,
                        help="Full path directory where outputs are written")
    parser.add_option("-p", "--output-prefix", dest="out_pre", type=str,
                        help="Output prefix")
    parser.add_option("-u", "--upstream-sequence", dest="upstream", type=str,
                        help="Excepected sequence upstream of insert")
    parser.add_option("-d", "--downstream-sequence", dest="downstream", type=str,
                            help="Excepected sequence downstream of insert")
    parser.add_option("-m", "--list-motifs", dest="motifs", type=str,
                        help="List of motifs expected in insert")
    parser.set_defaults(
        motifs="/shared/sudlab1/General/projects/SynthUTR_hepG2_a549/motif_merge_linkers/correct_strandedness/highstab_list_filtered_oligos.list",
        linker="GGCGAT",
        RE1site="GGCGCC", #NarI
        RE2site="ATCGAT", #ClaI
        out_dir=None,
        out_pre="out",
        upstream="GCGGCCGC", #NotI
        downstream="TCGTACG" #T+BsiWI
    )

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.start(parser, argv=argv)
    
    linker = options.linker #sequence linker
    #List of possible motifs
    motifs = open(options.motifs).readlines()
    motifs = [s.strip() for s in motifs]
    #Sorting it so it checks for matches with longer motifs first (less ambiguity)
    motifs = sorted(motifs, key=len, reverse = True)
    before = (options.upstream)
    after = (options.downstream)
    siteRE1 = (options.RE1site)
    siteRE2 = (options.RE2site)

    #Wrong sequences
    wrong_upstream = before + siteRE1
    wrong_downstream = siteRE2 + after
    no_insert = before + linker + after

    def check_motifs(sequence):
        for m in motifs:
            where_is_motif = regex.search("^(%s){s<=1}" % m, sequence)
            #Check if start with it
            if where_is_motif is None:
                continue
            else:
                return m
            
        
    def check_linker(sequence):
        where_is_linker = regex.search("^(%s){s<=1}" % linker, sequence)
        #Check if start with it
        if where_is_linker is None:
            return False
        else:
            #return end position of linker
            end_l = where_is_linker.span()[1]
            return end_l

    #Checking linker and motifs ambiguity 
    linker_ambiguity = []
    for m in motifs:
        matched = regex.search("(%s){s<=1}" % linker, m)
        if matched is None:
            continue
        else:
            linker_ambiguity.append(m)
    
    ###Motifs Content###
    SUTRs_content_ID_SUTR = []
    SUTRs_content_UMI_SUTR = []
    SUTRs_content_length = []
    SUTRs_content_NumberOfMotifs = []
    SUTRs_content_MotifsOrder = []

    SUTRs_failed_ID_SUTR = []
    SUTR_failed_problem = []

    SUTR_fastq =  SeqIO.parse(iotools.open_file(options.fastq), "fastq")

    #Loop through sequences
    n = 0
    ambiguous = 0
    for sequence in SUTR_fastq:
        n += 1
        E.debug("Processing read %i" % n)
        #Get sutr id and extract UMI
        seq_id = sequence.id
        umi_id = seq_id.split("_")[1]
        #Get sequence
        seq = str(sequence.seq)

        #Detect if any RE sites back to back, should not happen but who knows
        if regex.search("(%s%s){s<=1}" % (siteRE1,siteRE2), seq) is not None:
            #empty
            SUTRs_failed_ID_SUTR.append(seq_id)
            SUTR_failed_problem.append("Back to back RE sites")
            continue

        #Detect if wrongly inserted
        if (regex.search("(%s){s<=1}" % wrong_upstream, seq) is not None or
            regex.search("(%s){s<=1}" % wrong_downstream, seq) is not None):
            SUTRs_failed_ID_SUTR.append(seq_id)
            SUTR_failed_problem.append("RE sites present: wrong direction insertion")
            continue

        #Detect if nothing inserted
        if (regex.search("(%s){s<=1}" % no_insert, seq) is not None):
            SUTRs_failed_ID_SUTR.append(seq_id)
            SUTR_failed_problem.append("RE sites present: empty vector")
            continue

        #Detect sequence before 1st linker
        try:
            #Search first instance
            start = regex.search("(%s){s<=1}" % before, seq)
            #Start point of first match instance
            start_span = start.span()[1]
        except AttributeError:
            SUTRs_failed_ID_SUTR.append(seq_id)
            SUTR_failed_problem.append("Didn't detect sequence upstream of 1st linker")
            continue

        #Get sequence from 1st linker to last linker
        try:
            #Finding best match
            end = regex.search("(?b)(%s){s<=1}" % after, seq)
            #End point of end sequences' best match 
            end_span = end.span()[0]
            seq = seq[start_span:end_span]
        except AttributeError:
            SUTRs_failed_ID_SUTR.append(seq_id)
            SUTR_failed_problem.append("Didn't detect sequence downstream of insert")
            continue

        #Now loop in that sequence
        matched_motifs_seq = []
        number_motifs = 0
        length_utr = len(seq)
        while seq is not None:
            #print(seq)
            #First should be a linker
            end_linker = check_linker(seq) 
            if end_linker == False:
                 SUTRs_failed_ID_SUTR.append(seq_id)
                 SUTR_failed_problem.append("Linker ambiguity at %s" % seq)
                 E.debug((seq_id, " failed at %s" % seq))
                 ambiguous += 1
                 #Trashing sequence
                 seq = None
                 continue
            
            #Check if last linker
            if len(seq) == len(linker):
                SUTRs_content_length.append(length_utr)
                SUTRs_content_ID_SUTR.append(seq_id)
                SUTRs_content_UMI_SUTR.append(umi_id.encode())
                SUTRs_content_NumberOfMotifs.append(number_motifs)
                SUTRs_content_MotifsOrder.append(matched_motifs_seq)
                seq = None
                continue

            #Else continue loop
            else:
                seq = seq[end_linker:]
                #matched_motifs_seq.append("-")
            
            #Now check the motif
            motif_check = check_motifs(sequence = seq)
            if motif_check ==  None:
                SUTRs_failed_ID_SUTR.append(seq_id)
                SUTR_failed_problem.append("Motif ambiguity at %s" % seq)
                seq = None
                E.debug((seq_id, " failed at %s" % seq))
                ambiguous += 1
                continue
            else:
                matched_motifs_seq.append(str(motif_check))
                number_motifs += 1
                #print("Motif after linker", str(motif_check))
                seq = seq[len(motif_check):]

    data = {
        'ID_SUTR': SUTRs_content_ID_SUTR,
        'UMI_SUTR': SUTRs_content_UMI_SUTR,
        'length':SUTRs_content_length,
        'NumberOfMotifs':SUTRs_content_NumberOfMotifs,
        'MotifsOrder':SUTRs_content_MotifsOrder
    }
    SUTRs_content = pd.DataFrame(data)  

    data = {
        'read_name': SUTRs_failed_ID_SUTR,
        'problem':SUTR_failed_problem
    }

    SUTRs_failed = pd.DataFrame(data)

    SUTRs_failed.to_csv(os.path.join(options.out_dir, options.out_pre + "_failed.tsv"),
                        sep="\t",
                        index = False)

    E.debug("There was %i reads in fastq file, %i failed, %i for ambiguity in linker or motif sequence" % (n,len(SUTRs_failed), ambiguous))


    ###UMI deduplication###

    #UMI function
    clusterer = UMIClusterer(cluster_method="directional")

    #Motif as string sequence to group
    motif_sequences = []
    for i in SUTRs_content.index:
        motifs_seq = ''.join(SUTRs_content["MotifsOrder"][i])
        motif_sequences.append(motifs_seq)
    SUTRs_content.loc[:,"Sequence"] = motif_sequences

    #UMI clustering, dedup and count
    unique_sutrs = set(SUTRs_content.loc[:,("Sequence")])
    grouped = SUTRs_content.groupby(["Sequence"])
    avg_umi_per_utr = []
    sutr_counts = []
    counter = 0
    for s in unique_sutrs:
        counter += 1
        #Clustering
        this_group = grouped.get_group(s)
        group_dict = this_group.groupby(["UMI_SUTR"]).count().to_dict()['Sequence']
        E.debug("Group %i contains %i UMIs" %(counter, len(group_dict)))
        clustered_umis = clusterer(group_dict, threshold=1)
        clustered_umis = [umi[0] for umi in clustered_umis]
        avg_umi_per_utr.append(len(clustered_umis))
        #Dedup and count
        if len(clustered_umis) == 1:
            counts = len(this_group) 
        else:
            counts = len(this_group[~this_group["UMI_SUTR"].isin(clustered_umis)])+1 
        sutr_counts.append(["sUTR_%i" % counter,
                            ','.join(this_group.iloc[0]["MotifsOrder"]),
                            this_group.iloc[0]["length"],
                            this_group.iloc[0]["NumberOfMotifs"],
                            counts])


    avg_umi_per_utr = sum(avg_umi_per_utr) / len(avg_umi_per_utr)
    E.debug("On average, %i UMIs per UTR" % avg_umi_per_utr)
    

    sutr_counts_df = pd.DataFrame(sutr_counts, columns = ("sUTR_ID",
                                                          "MotifsOrder",
                                                          "length",
                                                          "NumberOfMotifs",
                                                          "counts"))
    
    sutr_counts_df.to_csv(os.path.join(options.out_dir, options.out_pre + "_sUTR_counts.tsv"),
                         sep="\t",
                         index = False)

    #Final log file 
    logging = ["Total number of reads: %i" % n,
               "Failed reads: %i" % len(SUTRs_failed),
               "Failed because of ambiguity in linker/motif: %i" % ambiguous,
               "Number of unique UTRs: %i" % len(sutr_counts_df),
               "Average UMIs per UTR: %i" % avg_umi_per_utr,
               '''List of motifs matching the linker 
                  with 1 mismatch: %s''' % ','.join(linker_ambiguity)]

    logfile = open(os.path.join(options.out_dir, options.out_pre + "_stats.txt"), "w")
    for l in logging:
        logfile.write(l + "\n")
    logfile.close()

    # write footer and output benchmark information.
    E.stop()


if __name__ == "__main__":
    sys.exit(main(sys.argv))
