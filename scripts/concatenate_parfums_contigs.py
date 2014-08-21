#!/usr/bin/env python                                                                                                                                                               

__author__ = "Molly Gibson"
__email__ = "molly.gibson@wustl.edu"


# Python imports                                                                                                                                                                   
import argparse, subprocess, os, itertools, operator
# Other imports                                                                                                                                                                    
import parse_config, parse_mapping

def main():
    parser = argparse.ArgumentParser(description="This script will concatenate all assembled contigs from PARFuMS and add sample name to each contig.")
 
    # arguments                                                                                                                                                                    
    parser.add_argument('-parfums_fp', dest="parfums_fp", help="Path to directory where PARFuMS assembly is run")
    parser.add_argument('-sample_ids', dest="sample_ids", help="If provided, will only concatenate these samples")
    parser.add_argument('-o', dest="output", type=argparse.FileType('w'),  help="Output file path")
    args = parser.parse_args()

    if args.sample_ids:
        samples_to_keep = [line.rstrip() for line in open(args.sample_ids)]

    for file_fp in os.listdir(args.parfums_fp):
        if os.path.isdir(args.parfums_fp.rstrip('/') + "/" + file_fp) and (not args.sample_ids or (args.sample_ids and file_fp in samples_to_keep)):
            sample_name = file_fp

            try:
                for line in open(args.parfums_fp.rstrip('/') + "/" + file_fp + '/' + sample_name +'.lastContigs.fna', 'r'):
                    if line.startswith(">"):
                        header = ">" + sample_name + "_" + line.rstrip().strip(">")
                        args.output.write(header + "\n")
                    else:
                        args.output.write(line.rstrip() + "\n")
            except Exception, e:
                print e
            

if __name__ == "__main__":
    main()
