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
    parser = argparse.ArgumentParser(description="This script takes in functional metagenomic contigs (or proteins), and forward and reverse reads and rarifies them to a certain read depth, and then maps the reads back to the contigs or ORFs and creates a final pileup.")
    parser.add_argument('-contigs', dest="contig_fp", help="Path to collapsed contig file")
    parser.add_argument('-orfs', dest="orf_fp", help="Path to nucleotides file")
    parser.add_argument('-for_pair', dest="forward_fp", help="Path to quality filtered forward reads paired")
    parser.add_argument('-rev_pair', dest="reverse_fp", help="Path to quality filetered reverse reads paired")
    parser.add_argument('-unpair', dest="unpaired_fp", help="Path to quality filtered reads not paired")
    parser.add_argument('-rare_paired', dest="paired_levl", help="Level to rarify paired reads")
    parser.add_argument('-rare_unpaired', dest="unpaired_levl", help="Level to rarify unpaired reads")
    parser.add_argument('-o', dest="output_fp", help="Output directory")
    args = parser.parse_args()

    # Make output directory
    if not os.path.exists(args.output_fp):
        cmd = "mkdir " + args.output_fp
        run_command(cmd)

    # Rareify our reads
    [unpaired_q_fp, forward_paired_q_fp, reverse_paired_q_fp] = rarefy_reads(args)

    # Align the rareified reads
    [pileup_file_paired, pileup_file_unpaired, alignment_fp] = align_reads(args, unpaired_q_fp, forward_paired_q_fp, reverse_paired_q_fp)

    [orf_cov_dict, orf_len_dict] = generate_final_pileup(pileup_file_paired, pileup_file_unpaired, alignment_fp)

    coverage_fp = open(args.output_fp + "/coverage.txt", 'w')
    coverage_fp.write("orf_name\ttotal_coverage\torf_len\tave_coverage\n")
    for orf in orf_len_dict.keys():
        if orf in orf_cov_dict.keys():
            coverage_fp.write(orf + "\t" + str(orf_cov_dict[orf]) + "\t" + str(orf_len_dict[orf]) + "\t" + str(float(orf_cov_dict[orf])/float(orf_len_dict[orf])) + "\n")
        else:
            coverage_fp.write(orf +"\t0\t"+ str(orf_len_dict[orf]) + "\t0\n")


def align_reads(args, unpaired_q_fp, forward_paired_q_fp, reverse_paired_q_fp):
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
            cmd = "bowtie2 -x " + index + " -U " + unpaired_q_fp + " -S " + sam_unpaired_fp + " --score-min L,0,-0.17 -a"
        elif args.orf_fp:
            cmd = "bowtie2 -x " + index + " -U " + unpaired_q_fp + " -S " + sam_unpaired_fp + " --score-min L,0,-0.17"
        run_command(cmd)
    sam_paired_fp = args.output_fp + "/paired.sam"
    if not os.path.exists(sam_paired_fp):
        if args.contig_fp:
            cmd = "bowtie2 -x " + index + " -1 " + forward_paired_q_fp + " -2 " + reverse_paired_q_fp + " -S " + sam_paired_fp + " --score-min L,0,-0.17 -I 250 -X 500 --no-discordant --no-contain -a"
        elif args.orf_fp:
            cmd = "bowtie2 -x " + index + " -1 " + forward_paired_q_fp + " -2 " + reverse_paired_q_fp + " -S " + sam_paired_fp + " --score-min L,0,-0.17 -I 250 -X 500 --no-discordant --no-contain"
        run_command(cmd)

    # Create pileups using samtools                                                                                                                                                                                                                               
    sort_bam_unpaired_fp = args.output_fp + "/unpaired_sorted"
    if not os.path.exists(sort_bam_unpaired_fp + ".bam"):
        cmd = "samtools view -buS " + args.output_fp + "/unpaired.sam | samtools sort -m 4000000000 - " + sort_bam_unpaired_fp
        run_command(cmd)
        cmd = "samtools index " + sort_bam_unpaired_fp + ".bam"
        run_command(cmd)

    sort_bam_paired_fp = args.output_fp + "/paired_sorted"
    if not os.path.exists(sort_bam_paired_fp + ".bam"):
        cmd = "samtools view -buS " + args.output_fp + "/paired.sam | samtools sort -m 4000000000 - " + sort_bam_paired_fp
        run_command(cmd)
        cmd = "samtools index " + sort_bam_paired_fp + ".bam"
        run_command(cmd)

    pileup_file_unpaired = args.output_fp + "/unpaired_" + str(args.unpaired_levl) +".pileup"
    if not os.path.exists(pileup_file_unpaired):
        cmd = "samtools mpileup -f " + alignment_fp + " " + sort_bam_unpaired_fp + ".bam > " + pileup_file_unpaired
        run_command(cmd)

    pileup_file_paired = args.output_fp + "/paired_" + str(args.paired_levl) +".pileup"
    if not os.path.exists(pileup_file_paired):
        cmd = "samtools mpileup -f " + alignment_fp + " " + sort_bam_paired_fp + ".bam > " + pileup_file_paired
        run_command(cmd)

    return [pileup_file_paired, pileup_file_unpaired, alignment_fp]


def generate_final_pileup(paired_pileup_fp, unpaired_pileup_fp, reference):

    # Import both paire and unpaired pileups
    paired = pandas.io.parsers.read_table(paired_pileup_fp, header=None)
    paired_coverage = paired.loc[:,[0,1,3]]

    unpaired = pandas.io.parsers.read_table(unpaired_pileup_fp, header=None)
    unpaired_coverage = unpaired.loc[:,[0,1,3]]

    # Join the two dataframes together on both the contig and the position
    index = [unpaired_coverage[0], unpaired_coverage[1]]
    tuples = list(zip(*index))
    unpaired_coverage.index = tuples

    index = [paired_coverage[0], paired_coverage[1]]
    tuples = list(zip(*index))
    paired_coverage.index = tuples

    paired_coverage.columns = ['contig1', 'position1', 'coverage1']
    unpaired_coverage.columns = ['contig2', 'position2', 'coverage2']

    joined_df = unpaired_coverage.join(paired_coverage, how="outer")
    coverage_df = joined_df.loc[:,['contig1', 'coverage1', 'coverage2']]

    # sum the coverage for paired and unpaired
    coverage_df['sum'] = coverage_df.sum(axis=1)

    orfs = coverage_df.groupby('contig1')
    
    cov_per_orf = orfs.sum()
    cov_per_orf_dict = cov_per_orf['sum'].to_dict()
            
    # Calcualte the total number of bases in the functional metagenomics file                                                                                                                                                                                     
    contig_metadata = pandas.io.parsers.read_table(reference + ".fai", header=None)
    contig_metadata.columns = ['contig_name', 'length', 'trash', 'trash2', 'trash3']
    orf_length_dict = contig_metadata.set_index('contig_name')['length'].to_dict()

    return [cov_per_orf_dict, orf_length_dict]

def rarefy_reads(args):
    # Rareify paired reads                                                                                                                                                                                                                                        
    forward_paired_q_fp = args.output_fp + "/forward_quality_paired_" + str(args.paired_levl) + ".fastq"
    reverse_paired_q_fp= args.output_fp + "/reverse_quality_paired_" + str(args.paired_levl) +".fastq"
    if not os.path.exists(forward_paired_q_fp) or not os.path.exists(reverse_paired_q_fp):
        # Read in forward and reverse reads                                                                                                                                                                                                                        
        for_file = open(args.forward_fp, "rU") # FORWARD                                                                                                                                                                                                           
        forward_paired_reads = list(SeqIO.parse(for_file, "fastq"))
        print "Forward reads read."
        sys.stdout.flush()
        rev_file = open(args.reverse_fp, "rU") # REVERSE                                                                                                                                                                                                           
        reverse_paired_reads = list(SeqIO.parse(rev_file, "fastq"))
        print "Reverse reads read."
        sys.stdout.flush()

        # pick indecies to keep                                                                                                                                                                                                                                    
        randIndex = random.sample(range(len(forward_paired_reads)), int(args.paired_levl))
        randIndex.sort()

        # forward reads rarefy                                                                                                                                                                                                                                     
        randIndex = random.sample(range(len(forward_paired_reads)), int(args.paired_levl))
        randIndex.sort()
        rare_forward_reads = [forward_paired_reads[i] for i in randIndex]
        SeqIO.write(rare_forward_reads, forward_paired_q_fp, "fastq")

        # reverse reads rarefy                                                                                                                                                                                                                                     
        rare_reverse_reads = [reverse_paired_reads[i] for i in randIndex]
        SeqIO.write(rare_reverse_reads, reverse_paired_q_fp, "fastq")

    # Rareify unpaired reads                                                                                                                                                                                                                                       
    unpaired_q_fp= args.output_fp + "/quality_unpaired_" + str(args.unpaired_levl) +".fastq"
    if not os.path.exists(unpaired_q_fp):
        unpaired_file = open(args.unpaired_fp, "rU")
        unpaired_reads = list(SeqIO.parse(unpaired_file, "fastq"))
        rare_unpaired_reads = random.sample(unpaired_reads, int(args.unpaired_levl))
        print "Unpaired\t" + str(len(unpaired_reads)) + "\t" + str(len(rare_unpaired_reads))
        SeqIO.write(rare_unpaired_reads, unpaired_q_fp, "fastq")

    return [unpaired_q_fp, forward_paired_q_fp, reverse_paired_q_fp]
        
def run_command(command):
    runCmd = subprocess.Popen(command,shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output,errors = runCmd.communicate()
    print command
    print output
    print errors



if __name__ == "__main__":
    main()
