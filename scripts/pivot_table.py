#!/usr/bin/env python                                                                                                                                                          

__author__ = "Molly Gibson"
__email__ = "molly.gibson@wustl.edu"


# Python imports
import argparse, subprocess, os, itertools, operator
# Other imports
from Bio import SeqIO
import parse_config, parse_mapping
import helper as h
import pandas

def main():
    parser = argparse.ArgumentParser(description="This script collapses all redundant proteins \
                        or ORFs by percent identity, can be done within either library or antibiotic.")
    # arguments
    parser.add_argument('-i', dest="input_fp", help="Input matrix")
    parser.add_argument('-o', dest="output_fp", help="Output matrix")
    parser.add_argument('-reverse', dest="reverse", action="store_true", help="Reverse pivot table")
    args = parser.parse_args()
    
    if args.reverse:
        in_table = pandas.io.parsers.read_table(args.input_fp, index_col=0)
        unstaked = in_table.unstack()
        unstaked.to_csv(args.output_fp, sep="\t")
    else:
        in_table = pandas.io.parsers.read_table(args.input_fp, header=False)
        in_table.columns = ("rows", "columns", "values")
        stacked = in_table.pivot(index='rows', columns='columns', values='values')
        final = stacked.fillna(value=0)
        final.to_csv(args.output_fp, sep="\t")
        

if __name__ == "__main__":
    main()



