#!/usr/bin/env python                                                                                                                                                       
__author__ = "Molly Gibson"
__email__ = "molly.gibson@wustl.edu"


# Python imports                                                                                                                                                                  
import argparse, subprocess, os, itertools, operator, random
from os.path import basename
from Bio import SeqIO

# Other imports                                                                                                                                                             
import parse_config, parse_mapping

def main():
    parser = argparse.ArgumentParser(description="This script takes in functional metagenomic contigs, and forward and reverse reads and rarifies them to a certain read depth, and then maps the reads back to the contigs and calcualtes average coverage")
    parser.add_argument('-contigs', dest="contig_fp", help="Path to collapsed contig file")
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

    # Rarify Reads
    print "Rarefying Reads"
    print "Dataset\tTotal\tRare"

    seed_number = random.randint(0,100)  # Use a seed so that the paired reads pick the same random subset

    forward_paired_q_fp = args.output_fp + "/forward_quality_paired_" + str(args.paired_levl) + ".fastq"
    if not os.path.exists(forward_paired_q_fp):
        for_file = open(args.forward_fp, "rU")
        forward_paired_reads = list(SeqIO.parse(for_file, "fastq"))
        random.seed(seed_number)
        rare_forward_reads = random.sample(forward_paired_reads, int(args.paired_levl))
        print "Forward Paired\t" + str(len(forward_paired_reads)) + "\t" + str(len(rare_forward_reads))
        SeqIO.write(rare_forward_reads, forward_paired_q_fp, "fastq")
    
    reverse_paired_q_fp= args.output_fp + "/reverse_quality_paired_" + str(args.paired_levl) +".fastq"
    if not os.path.exists(reverse_paired_q_fp):
        rev_file = open(args.reverse_fp, "rU")
        reverse_paired_reads = list(SeqIO.parse(rev_file, "fastq"))
        random.seed(seed_number)
        rare_reverse_reads = random.sample(reverse_paired_reads, int(args.paired_levl))
        print "Reverse Paired\t" + str(len(reverse_paired_reads)) + "\t" + str(len(rare_reverse_reads))
        SeqIO.write(rare_reverse_reads, reverse_paired_q_fp, "fastq")

    unpaired_q_fp= args.output_fp + "/quality_unpaired_" + str(args.unpaired_levl) +".fastq"
    if not os.path.exists(unpaired_q_fp):
        unpaired_file = open(args.unpaired_fp, "rU")
        unpaired_reads = list(SeqIO.parse(unpaired_file, "fastq"))
        rare_unpaired_reads = random.sample(unpaired_reads, int(args.unpaired_levl))
        print "Unpaired\t" + str(len(unpaired_reads)) + "\t" + str(len(rare_unpaired_reads))
        SeqIO.write(rare_unpaired_reads, unpaired_q_fp, "fastq")
                           
    # Align all of your reads to the contigs
    index = os.path.splitext(args.contig_fp)[0]
    sam_unpaired_fp = args.output_fp + "/unpaired.sam"
    if not os.path.exists(sam_unpaired_fp):
        cmd = "bowtie2 -x " + index + " -U " + unpaired_q_fp + " -S " + sam_unpaired_fp + " --score-min L,0,-0.17"
        run_command(cmd)
    sam_paired_fp = args.output_fp + "/paired.sam"
    if not os.path.exists(sam_paired_fp):
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
        cmd = "samtools mpileup -f " + args.contig_fp + " " + sort_bam_unpaired_fp + ".bam > " + pileup_file_unpaired
        run_command(cmd)

    pileup_file_paired = args.output_fp + "/paired_" + str(args.paired_levl) +".pileup"
    if not os.path.exists(pileup_file_paired):
        cmd = "samtools mpileup -f " + args.contig_fp + " " + sort_bam_paired_fp + ".bam > " + pileup_file_paired
        run_command(cmd)

def run_command(command):
    runCmd = subprocess.Popen(command,shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output,errors = runCmd.communicate()
    print command
    print output
    print errors



if __name__ == "__main__":
    main()
