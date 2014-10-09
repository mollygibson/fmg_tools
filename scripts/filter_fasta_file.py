#!/usr/bin/env python                                                                                                                                                                                                                                                                                              

__author__ = "Molly Gibson"
__email__ = "molly.gibson@wustl.edu"


# Python imports                                                                                                                                                                                                                 
import argparse, subprocess, os, itertools, operator
from Bio import SeqIO
# Other imports                                                                                                                                                                     
import helper as h
import parse_config, parse_mapping

def main():
    parser = argparse.ArgumentParser(description="This script takes in a nucleotide file and filters it by both length and if it is annotated as a resistance gene or not.")
    parser.add_argument('-fasta', dest="nucl_fp", help="Path to nucleotide or protein fasta file")
    parser.add_argument('-min_len', dest="min_len", help="Minimum length gene to retain") 
    parser.add_argument('-annotations', dest="annotation_fp", help="Path to annotations file (only necessary if filtering by resfams) - this is the *.final.txt annotation file")
    parser.add_argument('--res_filter', dest="res_filter", help="Flag to determine if keeping only resistance genes or all genes.", action='store_true', default=False)
    parser.add_argument('-ids', dest="ids", help="gene ids to keep")
    parser.add_argument('-o', dest="output_fp", help="Path to filtered fasta output directory")
    args = parser.parse_args()

    if args.min_len:
        min_len = args.min_len
    else:
        min_len = 0

    if not os.path.exists(args.output_fp):
        h.run_command("mkdir " + args.output_fp)

    filtered_file = args.output_fp.rstrip('/') + "/filtered.fasta"
    other_file = args.output_fp.rstrip('/') + "/other.fasta"

    # Determine genes that have resfams annotation
    if args.res_filter:
        resfam_genes = find_res_genes(args.annotation_fp)
        # Filter fasta file and output
        SeqIO.write(filter_records(args.nucl_fp, min_len, resfams=resfam_genes)[0], filtered_file, "fasta")
        SeqIO.write(filter_records(args.nucl_fp, min_len, resfams=resfam_genes)[1], other_file, "fasta")
    elif args.ids:
        SeqIO.write(filter_records(args.nucl_fp, min_len, ids=[line.rstrip() for line in open(args.ids)])[0], filtered_file, "fasta")
        SeqIO.write(filter_records(args.nucl_fp, min_len, ids=[line.rstrip() for line in open(args.ids)])[1], other_file, "fasta")
    else:
        SeqIO.write(filter_records(args.nucl_fp, min_len)[0], filtered_file, "fasta")
        SeqIO.write(filter_records(args.nucl_fp, min_len)[1], other_file, "fasta")
        

def find_res_genes(annotations):
    resfam_ids = []
    resfam = False
    for line in open(annotations, 'r'):
        if line.startswith('>'):
            if resfam:
                resfam_ids.append(gene)
            resfam = False
            contig_name = line.split()[0].strip('>') + "." + line.split("Contig:")[1].split()[0]
        elif not line.startswith("\t\t"):
            if resfam:
                resfam_ids.append(gene)
            resfam = False
            gene = contig_name + ":" + line.split("\t")[1] + "-" + line.split("\t")[2]
        else:
            if 'ResFam' in line:
                resfam = True
    if resfam:
        resfam_ids.append(gene)

    return resfam_ids

def filter_records(fasta_fp, min_len,resfams="",ids=""):
    sequences = SeqIO.parse(fasta_fp, 'fasta')
    filtered_seqs = []
    other_seqs = []
    for record in sequences:
        if len(record.seq) >= int(min_len):
            if ids:
                print record.id
                if record.id in ids:
                    filtered_seqs.append(record)
                else:
                    other_seqs.append(record)
            elif not resfams or record.id in resfams:
                filtered_seqs.append(record)
            else:
                other_seqs.append(record)
    return [filtered_seqs,other_seqs]

if __name__ == "__main__":
    main()


