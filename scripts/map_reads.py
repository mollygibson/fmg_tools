#!/usr/bin/env python                                                                                                                                                       
__author__ = "Molly Gibson"
__email__ = "molly.gibson@wustl.edu"

# Python imports                                                                                                                                                                  
import argparse, subprocess, os, itertools, operator, random, sys
from os.path import basename
from Bio import SeqIO 
import pandas

# Other imports                                                                                                                                                             
import parse_config, parse_mapping

def main():
    parser = argparse.ArgumentParser(description="This script takes in functional metagenomic contigs (or ORFs), and forward and reverse reads and then maps the reads back to the contigs or ORFs and creates a final pileup, finally calculating an RPKM of sample by contig or ORF")
    parser.add_argument('-contigs', dest="contig_fp", help="Path to contig file to map reads")
    parser.add_argument('-orfs', dest="orf_fp", help="Path to nucleotides ORF file to map reads")
    parser.add_argument('-for_pair', dest="forward_fp", help="Path to quality filtered forward reads paired")
    parser.add_argument('-rev_pair', dest="reverse_fp", help="Path to quality filetered reverse reads paired")
    parser.add_argument('-unpair', dest="unpaired_fp", help="Path to quality filtered reads not paired")
    parser.add_argument('-o', dest="output_fp", help="Output directory")
    args = parser.parse_args()

    # Make output directory
    if not os.path.exists(args.output_fp):
        cmd = "mkdir " + args.output_fp
        run_command(cmd)

    # Align the reads
    [pileup_file_paired, pileup_file_unpaired, sort_bam_unpaired_fp, sort_bam_paired_fp, alignment_fp] = align_reads(args)

    # get total reads
    unpaired_stats = open(args.output_fp + "/unpaired_alignment_stats.txt", 'r')
    unpaired_count = unpaired_stats.readline().split()[0]
    paired_stats = open(args.output_fp + "/paired_alignment_stats.txt", 'r')
    paired_count = paired_stats.readline().split()[0]
    total_reads = int(unpaired_count) + int(paired_count)

    # Calcualte RPKM
    calculate_rpkm(args, sort_bam_paired_fp, sort_bam_unpaired_fp, total_reads)

def align_reads(args):
    # Determine which file we will align rarefied reads to                                                                                                                                                                                                        
    if args.contig_fp:
        alignment_fp = args.contig_fp
    elif args.orf_fp:
        alignment_fp = args.orf_fp

    # Align all of your reads to the contigs                                                                                                                                                                                                                      
    index = os.path.splitext(alignment_fp)[0]
    sam_unpaired_fp = args.output_fp + "/unpaired.sam"
    if not os.path.exists(sam_unpaired_fp):
        if args.contig_fp:
            cmd = "bowtie2 -x " + index + " -U " + args.unpaired_fp + " -S " + sam_unpaired_fp + " --end-to-end 2>&1 | tee " + args.output_fp + "/unpaired_alignment_stats.txt"
        elif args.orf_fp:
            cmd = "bowtie2 -x " + index + " -U " + args.unpaired_fp + " -S " + sam_unpaired_fp + " --end-to-end 2>&1 | tee " + args.output_fp + "/unpaired_alignment_stats.txt"
        run_command(cmd)
    sam_paired_fp = args.output_fp + "/paired.sam"
    if not os.path.exists(sam_paired_fp):
        if args.contig_fp:
            cmd = "bowtie2 -x " + index + " -1 " + args.forward_fp + " -2 " + args.reverse_fp + " -S " + sam_paired_fp + " --end-to-end -I 250 -X 500 --no-discordant --no-contain 2>&1 | tee " + args.output_fp + "/paired_alignment_stats.txt"
        elif args.orf_fp:
            cmd = "bowtie2 -x " + index + " -1 " + args.forward_fp + " -2 " + args.reverse_fp + " -S " + sam_paired_fp + " --end-to-end -I 250 -X 500 --no-contain 2>&1 | tee " + args.output_fp + "/paired_alignment_stats.txt"
        run_command(cmd)

    # Create pileups using samtools                                                                                                                                                                                                                               
    sort_bam_unpaired_fp = args.output_fp + "/unpaired_sorted"
    if not os.path.exists(sort_bam_unpaired_fp + ".bam"):
        cmd = "samtools view -buS " + args.output_fp + "/unpaired.sam | samtools sort -m 500000000 - " + sort_bam_unpaired_fp
        run_command(cmd)
        cmd = "samtools index " + sort_bam_unpaired_fp + ".bam"
        run_command(cmd)

    sort_bam_paired_fp = args.output_fp + "/paired_sorted"
    if not os.path.exists(sort_bam_paired_fp + ".bam"):
        cmd = "samtools view -buS " + args.output_fp + "/paired.sam | samtools sort -m 500000000 - " + sort_bam_paired_fp
        run_command(cmd)
        cmd = "samtools index " + sort_bam_paired_fp + ".bam"
        run_command(cmd)

    pileup_file_unpaired = args.output_fp + "/unpaired.pileup"
    if not os.path.exists(pileup_file_unpaired):
        cmd = "samtools mpileup -f " + alignment_fp + " " + sort_bam_unpaired_fp + ".bam > " + pileup_file_unpaired
        run_command(cmd)

    pileup_file_paired = args.output_fp + "/paired.pileup"
    if not os.path.exists(pileup_file_paired):
        cmd = "samtools mpileup -f " + alignment_fp + " " + sort_bam_paired_fp + ".bam > " + pileup_file_paired
        run_command(cmd)

    return [pileup_file_paired, pileup_file_unpaired, sort_bam_unpaired_fp, sort_bam_paired_fp, alignment_fp]


def calculate_rpkm(args, sort_bam_paired_fp, sort_bam_unpaired_fp, total_reads):
    
    paired_idxstats = args.output_fp + "/paired_idxstats.txt"
    unpaired_idxstats = args.output_fp + "/unpaired_idxstats.txt"

    cmd = "samtools idxstats " + sort_bam_paired_fp + ".bam > " + paired_idxstats
    run_command(cmd)
    cmd = "samtools idxstats " + sort_bam_unpaired_fp + ".bam > " + unpaired_idxstats
    run_command(cmd)

    paired = pandas.io.parsers.read_table(paired_idxstats, header = None)
    unpaired = pandas.io.parsers.read_table(unpaired_idxstats, header = None)

    index = [paired[0], unpaired[0]]
    tuples = list(zip(*index))
    paired.index = tuples

    index = [unpaired[0], paired[0]]
    tuples = list(zip(*index))
    unpaired.index = tuples

    paired.columns = ['orf1', 'len1', 'mapped1', 'unmapped1']
    unpaired.columns = ['orf2', 'len2', 'mapped2', 'unmapped2']

    joined_df = unpaired.join(paired, how="outer")
    num_mapped = joined_df.loc[:,['orf1', 'len1', 'mapped1', 'mapped2']]

    num_mapped['total_mapped'] = num_mapped['mapped1'] + num_mapped['mapped2']

    num_mapped['RPKM'] = (num_mapped['total_mapped'] * 10000000000)/(total_reads*num_mapped['len1'])

    final_df = num_mapped.loc[:,['orf1', 'len1', 'total_mapped', 'RPKM']]

    final_df.to_csv(args.output_fp + "/rpkm.txt", sep="\t", index=False)



def run_command(command):
    runCmd = subprocess.Popen(command,shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output,errors = runCmd.communicate()
    print command
    print output
    print errors

if __name__ == "__main__":
    main()
