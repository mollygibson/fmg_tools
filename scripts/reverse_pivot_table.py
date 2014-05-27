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
    args = parser.parse_args()
    
    in_table = pandas.io.parsers.read_table(args.input_fp, index_col=0)
    unstaked = in_table.unstack()
    unstaked.to_csv(args.output_fp, sep="\t")

if __name__ == "__main__":
    main()



