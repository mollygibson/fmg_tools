#!/usr/bin/env python                                                                                                                                                       
__author__ = "Molly Gibson"
__email__ = "molly.gibson@wustl.edu"

# Python imports                                                                                                                                                                  
import argparse, subprocess, os, os.path, itertools, operator, random, sys
from os.path import basename
from Bio import SeqIO 
import pandas

# Other imports                                                                                                                                                             
import parse_config, parse_mapping

def main():
    parser = argparse.ArgumentParser(description="This script takes in functional metagenomic contigs, and forward and reverse reads and rarifies them to a certain read depth, and then maps the reads back to the contigs and calcualtes average coverage")
    parser.add_argument('-i', dest="fastq_fp", help="Path to fastq file to rarefy")
    parser.add_argument('-dir', dest="directory", help="Path to list of fastq files")
    parser.add_argument('-level', dest="level", help="Level to rarify reads")
    parser.add_argument('-o', dest="out_fp", help="output")
    args = parser.parse_args()


    if args.fastq_fp:
        if args.out_fp:
            rareified_fp = args.out_fp
        else:
            rareified_fp = args.fastq_fp.split('.')[0] + "_" + str(args.level) + ".fasta"
        rarefy(args.fastq_fp, rareified_fp, args.level)
    elif args.directory:
        num_files = len([name for name in os.listdir(args.directory) if os.path.isfile(name)])
        each_level = int(args.level)/int(num_files)
        print args.level-(each_level*num_files)
        for filename in os.listdir(args.directory):
            rarefy(filename, filename + "_rare", each_level)

def rarefy(in_fp, out_fp, level):
    fastq_file = open(in_fp, "rU") 
    reads = list(SeqIO.parse(fastq_file, "fastq"))

    randIndex = random.sample(range(len(reads)), int(level))
    randIndex.sort()

    rare_reads = [reads[i] for i in randIndex]
    SeqIO.write(rare_reads, out_fp, "fasta")


if __name__ == "__main__":
    main()
