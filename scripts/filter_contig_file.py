#!/usr/bin/env python                                                                                                                                                                                                                                                                                              

__author__ = "Molly Gibson"
__email__ = "molly.gibson@wustl.edu"

# Python imports                                                                                                                                                                                                                 
import argparse, subprocess, os, itertools, operator
from Bio import SeqIO
# Other imports                                                                                                                                                                                                                  
import parse_config, parse_mapping

def main():
    parser = argparse.ArgumentParser(description="This script takes in a contig file from PARFuMS and filters it by length and if it contains an annotated resistance gene or not.")
    parser.add_argument('-contigs', dest="contig_fp", help="Path to contig fasta file")
    parser.add_argument('-min_len', dest="min_len", help="Minimum length gene to retain") 
    parser.add_argument('-annotations', dest="annotation_fp", help="Path to annotations file (only necessary if filtering by resfams)")
    parser.add_argument('--res_filter', dest="res_filter", help="Flag to determine if keeping only contigs that contain resistance genes or all contigs.", action='store_true', default=False)
    parser.add_argument('-o', dest="output_fp", help="Path to filtered fasta")
    args = parser.parse_args()

    # Determine genes that have resfams annotation
    resfam_contigs = find_res_genes(args.annotation_fp)

    # Filter fasta file and output
    SeqIO.write(filter_records(args.contig_fp, args.min_len, resfam_contigs), args.output_fp, "fasta")

def find_res_genes(annotations):
    resfam_ids = []
    resfam = False
    for line in open(annotations, 'r'):
        if line.startswith(">"):
            if resfam:
                resfam_ids.append(contig_name)
            resfam = False
            contig_name = line.split()[0].strip('>') + "_Contig_" + line.split("Contig:")[1].split()[0]
        elif 'ResFam' in line:
                resfam = True
    if resfam:
        resfam_ids.append(contig_name)

    print resfam_ids
    return resfam_ids

def filter_records(fasta_fp, min_len,resfams):
    sequences = SeqIO.parse(fasta_fp, 'fasta')
    filtered_seqs = []
    for record in sequences:
        contig_id = record.id.split("_Mean")[0].strip('>')
        if len(record.seq) >= int(min_len) and contig_id in resfams:
            filtered_seqs.append(record)
    return filtered_seqs

if __name__ == "__main__":
    main()


