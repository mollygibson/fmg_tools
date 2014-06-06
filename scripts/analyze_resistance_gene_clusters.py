#!/usr/bin/env python                                                                                                                                                                                                                        
__author__ = "Molly Gibson"
__email__ = "molly.gibson@wustl.edu"

# Python imports                                                                                                                                                                                                                             
import argparse, subprocess, os, itertools, operator
# Other imports                                                                                                                                                                                                                              
from Bio import SeqIO
import parse_config, parse_mapping
import helper as h

def main():
    parser = argparse.ArgumentParser(description="Analyze clusters")

    # arguments                                                     
    parser.add_argument('-clusters', dest="cluster_fp", help="Path to .clstr data")
    parser.add_argument('-fasta', dest="fasta_fp", help="Path to fasta of genes or proteins we care about")
    parser.add_argument('-o', dest="output_fp", help="Path to output directory")
    parser.add_argument('-f', dest="override", help="Force override otuput directory", action="store_true", default=False)
    args = parser.parse_args()

    res_genes = []
    for line in open(args.fasta_fp, 'r'):
        if line.startswith(">"):
            res_genes.append(line.split()[0].strip(">"))

    clusters = {}
    first = True
    save_cluster = False
    for line in open(args.cluster_fp, 'r'):
        if line.startswith(">"):
            if not first and save_cluster:
                clusters[cluster_num] = cluster_data
            # Reset data
            cluster_data = []
            cluster_num = line.rstrip().strip(">")
            # Reset flags
            first = False
            save_cluster = False
        else:
            gene_member = line.split(">")[1].split("...")[0]
            if gene_member in res_genes:
                save_cluster = True
            cluster_data.append(gene_member)
    if save_cluster:
        clusters[cluster_num] = cluster_data

    output = open(args.output_fp, 'w')
    for cluster_num in clusters:
        individuals = []
        for indv in clusters[cluster_num]:
            pieces = indv.split("_")
            individual = "_".join(pieces[0:2])
            if not individual in individuals:
                individuals.append(individual)

        output.write(cluster_num + "\t" + str(len(clusters[cluster_num])) + "\t" + str(len(individuals)) + "\t" +  ";".join(clusters[cluster_num]) + "\n")

if __name__ == "__main__":
    main()
