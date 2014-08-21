#!/usr/bin/env python                                                                                                                                                          

__author__ = "Molly Gibson"
__email__ = "molly.gibson@wustl.edu"


# Python imports
import argparse, subprocess, os, operator, sys, time
# Other imports
import helper as h
import pandas as pd


def main():
    parser = argparse.ArgumentParser(description="This script collapses all redundant proteins \
                        or ORFs by percent identity, can be done within either library or antibiotic.")
    # arguments
    parser.add_argument('-for', dest="forward_fp", help="Path to forward reads")
    parser.add_argument('-rev', dest="reverse_fp", help="Path to reverse reads")
    parser.add_argument('-m', dest="mapping_fp", help="Path to mapping file from barcode to sequence")
    parser.add_argument('-o', dest="output_fp", help="Path to output directory")
    args = parser.parse_args()

    start_time = time.time()
    
    if not os.path.exists(args.output_fp):
        h.run_command("mkdir " + args.output_fp)

    # make a mapping file
    mapping = pd.io.parsers.read_csv(args.mapping_fp, sep="\t",header=False)
    map_dict = mapping.set_index('barcode')['id'].to_dict()

    # middle 6 
    map_dict_mid6 = {}
    for barcode in map_dict.keys():
        barcode_mid6 = barcode[1:7]
        map_dict_mid6[barcode_mid6] = map_dict[barcode]

    # make dictionary of files
    forward_files = {}
    reverse_files = {}
    for sample_id in map_dict_mid6.values():
        forward_files[sample_id] = open(args.output_fp + "/Forward_" + sample_id + ".fastq", 'a')
        reverse_files[sample_id] = open(args.output_fp + "/Reverse_" + sample_id + ".fastq", 'a')

    forward_files["unknown"] = open(args.output_fp + "/Forward_unknown.fastq", 'a')
    reverse_files["unknown"] = open(args.output_fp + "/Reverse_unknown.fastq", 'a')

    # read counters
    total_reads = 0
    barcode_reads = {}

    # Loop through forward and reverse reads at the same time
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
                
                # check to see if it matches a barcode                                                                                                                                                                                       
                barcode_match = "unknown"
                if fwd_seq[1:7] in map_dict_mid6.keys() and rev_seq[1:7] in map_dict_mid6.keys():
                    # found a match                                                                                                                                                                                                          
                    barcode_match = fwd_seq[1:7]
                if not barcode_match in barcode_reads.keys():
                    barcode_reads[barcode_match] = 1
                else:
                    barcode_reads[barcode_match] = barcode_reads[barcode_match] + 1
                        
                # throw this line away
                trash = next(fwd)
                trash = next(rev)
                
                # quality scores
                fwd_qual = next(fwd).rstrip()
                rev_qual = next(rev).rstrip()
                
                # get the sample id for the barcode
                sample_id = "unknown"
                if barcode_match in map_dict_mid6.keys():
                    sample_id = map_dict_mid6[barcode_match]

                # Print out read to forward file
                seq = fwd_seq[8:len(fwd_seq)]
                qual = fwd_qual[8:len(fwd_qual)]
                if len(seq) == len(qual):
                    forward_files[sample_id].write(fwd_name + "\n" + fwd_seq[8:len(fwd_seq)] + "\n+\n" + fwd_qual[8:len(fwd_qual)] + "\n")
                else:
                    print("Sequence and quality not same length for forward: " + fwd_name)
                    
                # Print out the read tot her reverse file
                seq = rev_seq[8:len(fwd_seq)]
                qual = rev_qual[8:len(fwd_qual)]
                if len(seq) == len(qual):
                    reverse_files[sample_id].write(rev_name + "\n" + rev_seq[8:len(rev_seq)] + "\n+\n" + rev_qual[8:len(rev_qual)] + "\n")
                else:
                    print("Sequence and quality not same length for reverse: " + rev_name) 

                if total_reads % 1000000 == 0:
                    current_time = time.time() - start_time
                    print(str(total_reads) + " processed - " + str(current_time) + " seconds")
                    sys.stdout.flush()
                
            except StopIteration:
                break

    # Ourput stats
    with open(args.output_fp + "/distribution_stats.txt", 'w') as dist_out:
        for barcode in map_dict_mid6.keys():
            if barcode in barcode_reads.keys():
                dist_out.write(barcode + "\t" + map_dict_mid6[barcode] + "\t" + str(barcode_reads[barcode]) + "\n")
            else:
                dist_out.write(barcode + "\t" + map_dict_mid6[barcode] + "\t0\n")

        dist_out.write("total_reads\t\t" + str(total_reads) + "\n")

if __name__ == "__main__":
    main()


