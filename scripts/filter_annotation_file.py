#!/usr/bin/env python                                                                                                                                                                                                                                                                                              

__author__ = "Molly Gibson"
__email__ = "molly.gibson@wustl.edu"


# Python imports                                                                                                                                                                                                                 
import argparse, subprocess, os, itertools, operator, pandas
# Other imports                                                                                                                                                                                                                  
import parse_config, parse_mapping

def main():
    parser = argparse.ArgumentParser(description="This script takes in a filtered contig file and filters annotations")
    parser.add_argument('-proteins', dest="protein_fp", help="Path to collapsed protein file")
    parser.add_argument('-anno_tab', dest="ano_tab_fp", help="Path to annotation file")
    parser.add_argument('-o', dest="output_fp", help="Path to filtered annotations")
    args = parser.parse_args()

    proteins_to_keep = []
    for line in open(args.protein_fp, 'r'):
        if line.startswith('>'):
            proteins_to_keep.append(line.rstrip(">").split()[0])

    filtered_df = pandas.DataFrame()
    if args.ano_tab_fp:
        annotations = pandas.io.parsers.read_csv(args.ano_tab_fp, sep="\t")
        for index, row in annotations.iterrows():
            contig_pieces = row['contig'].split()
            protein = contig_pieces[0] + "." + contig_pieces[2].split(":")[1] + ":" + str(row['start']) + "-" + str(row['stop'])
            
            if protein in proteins_to_keep:
                filtered_df = filtered_df.append(row)
            else:
                print row['contig'] + "\t" + str(row['start']) + "\t" + str(row['stop'])

    filtered_df.to_csv(args.output_fp, sep="\t")

if __name__ == "__main__":
    main()


