#!/usr/bin/env python                                                                                                                                                          

__author__ = "Molly Gibson"
__email__ = "molly.gibson@wustl.edu"


# Python imports
import argparse, subprocess, os, operator, sys, time
# Other imports
import helper as h
import pandas as pd


def main():
    parser = argparse.ArgumentParser(description="This script just takes blast output format 6 and creates file with pid and query coverage.")
    # arguments
    parser.add_argument('-outfmt6', dest="blast_fp", help="Path to outfmt6 file")
    parser.add_argument('-o', dest="output_fp", help="Path to output file")
    args = parser.parse_args()

#    out = open(args.output_fp, 'a')

  
    blast = pd.io.parsers.read_table(args.blast_fp, header=False)
    blast.columns = ['Query', 'Target', 'Pid', 'Alignment_len', 'Mismatches', 'Gaps', 'Query_start', 'Query_end', 'Target_start', 'Target_end', 'Eval', 'Bit']
 
    blast['Coverage'] = (blast['Query_end']-blast['Query_start']+1) / blast['Alignment_len']
 
    blast.to_csv(args.output_fp, sep="\t", columns=['Query', 'Target', 'Pid', 'Coverage'])
    


if __name__ == "__main__":
    main()


