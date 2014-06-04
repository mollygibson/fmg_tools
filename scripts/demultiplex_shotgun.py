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

    total_reads = 0
    barcode_reads = {}

    # Keep track of what line we are on
    found = False
    print_reads = False
    # Loop through forward and reverse reads at the same time
    with open(args.forward_fp) as fwd, open(args.reverse_fp) as rev:
        while True:
            try:
                fwd_name = next(fwd).rstrip()
                rev_name = next(rev).rstrip()
                
                total_reads = total_reads + 1
                
                fwd_seq = next(fwd).rstrip()
                rev_seq = next(rev).rstrip()
                
                # check to see if it matches a barcode                                                                                                                                                                                       
                if fwd_seq[1:7] in map_dict_mid6.keys() and rev_seq[1:7] in map_dict_mid6.keys():
                    # found a match                                                                                                                                                                                                          
                    found = True
                    barcode_match = fwd_seq[1:7]
                    if not barcode_match in barcode_reads.keys():
                        barcode_reads[barcode_match] = 1
                    else:
                        barcode_reads[barcode_match] = barcode_reads[barcode_match] + 1
                else:
                    if not "unknown" in barcode_reads.keys():
                        barcode_reads["unknown"] = 1
                    else:
                        barcode_reads["unknown"] = barcode_reads["unknown"] + 1
                        
                trash = next(fwd)
                trash = next(rev)
                
                fwd_qual = next(fwd).rstrip()
                rev_qual = next(rev).rstrip()
                
                sample_id = "unknown"
                if found:
                    sample_id = map_dict_mid6[barcode_match]
                    # Print out read to forwad and reverse files
                    with open(args.output_fp + "/Forward_" + sample_id + ".fastq", 'a') as forward_out:
                        seq = fwd_seq[8:len(fwd_seq)]
                        qual = fwd_qual[8:len(fwd_qual)]
                        if len(seq) == len(qual):
                            forward_out.write(fwd_name + "\n" + fwd_seq[8:len(fwd_seq)] + "\n+\n" + fwd_qual[8:len(fwd_qual)] + "\n")
                        else:
                            print("Sequence and quality not same length for forward: " + fwd_name)
                    with open(args.output_fp + "/Reverse_" + sample_id + ".fastq", 'a') as reverse_out:
                        seq = rev_seq[8:len(fwd_seq)]
                        qual = rev_qual[8:len(fwd_qual)]
                        if len(seq) == len(qual):
                            reverse_out.write(rev_name + "\n" + rev_seq[8:len(rev_seq)] + "\n+\n" + rev_qual[8:len(rev_qual)] + "\n")
                        else:
                            print("Sequence and quality not same length for reverse: " + rev_name) 

                if total_reads % 100000 == 0:
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

if __name__ == "__main__":
    main()


