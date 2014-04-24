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
    parser = argparse.ArgumentParser(description="This script takes the coverage calculated from reads aligned to resistance gene and the profile of selections identified in and calculates the relative resistance potential of a microbiome.")
    parser.add_argument('-m', dest="resmap_fp", help="Path to gene to resistance profile map.")
    parser.add_argument('-coverage', dest="coverage_fp", help="Path to coverage file from rarefy_and_map_reads.py")
    parser.add_argument('-o', dest="output_fp", help="Path to output file")
    args = parser.parse_args()

    # import resistance map
    resistance_map = pandas.io.parsers.read_table(args.resmap_fp)
    res_profile_dictionary = resistance_map.set_index('orf_name')['res_profile'].to_dict()

    # import coverage
    coverage = pandas.io.parsers.read_table(args.coverage_fp)
    orf_coverage_dictionary = coverage.set_index('orf_name')['ave_coverage'].to_dict()

    # calc abx coverage
    abx_coverage = {}
    for orf in orf_coverage_dictionary.keys():
        resistance_list = res_profile_dictionary[orf].split(',')
        num_resistance_genes = len(resistance_list)
        for abx in resistance_list:
            if abx in abx_coverage.keys():
                abx_coverage[abx] = abx_coverage[abx] + (orf_coverage_dictionary[orf]/num_resistance_genes)
            else:
                abx_coverage[abx] = (orf_coverage_dictionary[orf]/num_resistance_genes)
    total = sum(abx_coverage.values())

    #output
    out = open(args.output_fp, 'w')
    out.write('abx\tcoverage\tnormalize_coverage\n')
    for abx in abx_coverage.keys():
        out.write(abx + "\t" + str(abx_coverage[abx]) + "\t" + str(abx_coverage[abx]/total) + "\n")
    


if __name__ == "__main__":
    main()
