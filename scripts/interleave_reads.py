#!/usr/bin/env python                                                                                                                                                          

__author__ = "Molly Gibson"
__email__ = "molly.gibson@wustl.edu"


# Python imports
import argparse, subprocess, os, operator, sys, time
# Other imports
import helper as h
import pandas as pd


def main():
    parser = argparse.ArgumentParser(description="This script is to be used in conjunction with demultiplex_reads.py \
                                      and will interleave forward and reverse reads.")
    # arguments
    parser.add_argument('-for', dest="forward_fp", help="Path to forward reads")
    parser.add_argument('-rev', dest="reverse_fp", help="Path to reverse reads")
    parser.add_argument('-o', dest="output_fp", help="Path to output file")
    args = parser.parse_args()

    start_time = time.time()
    
    out = open(args.output_fp, 'a')

    # Loop through forward and reverse reads at the same time
    total_reads = 0
    with open(args.forward_fp) as fwd, open(args.reverse_fp) as rev:
        while True:
            try:
                # iterate total reads
                total_reads = total_reads + 1

                # Read name
                fwd_name = next(fwd).rstrip()
                rev_name = next(rev).rstrip()

                # Sequence
                fwd_seq = next(fwd).rstrip()
                rev_seq = next(rev).rstrip()
                
                # throw this line away
                trash = next(fwd).rstrip()
                trash = next(rev).rstrip()
                
                # quality scores
                fwd_qual = next(fwd).rstrip()
                rev_qual = next(rev).rstrip()
                

                out.write(fwd_name + "\n" + fwd_seq + "\n+\n" + fwd_qual + "\n" + rev_name + "\n" + rev_seq + "\n+\n" + rev_qual + "\n")

                if total_reads % 1000000 == 0:
                    current_time = time.time() - start_time
                    print(str(total_reads) + " processed - " + str(current_time) + " seconds")
                    sys.stdout.flush()
                
            except StopIteration:
                break


if __name__ == "__main__":
    main()


