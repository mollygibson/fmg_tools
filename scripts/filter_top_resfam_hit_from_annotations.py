#!/usr/bin/env python                                                                                                                                                                                                                        
__author__ = "Molly Gibson"
__email__ = "molly.gibson@wustl.edu"

# Python imports                                                                                                                                                                                                                             
import argparse, subprocess, os, itertools, operator
# Other imports                                                                                                                                                                                                                              
import helper as h
import pandas

def main():
    parser = argparse.ArgumentParser(description="Analyze clusters for gene sharing and co-resistance")

    # arguments                                                     
    parser.add_argument('-i', dest="input_fp", help="Path to annotation tab file")
    parser.add_argument('-o', dest="output_fp", help="Path to resfam filtered annotation tab file")
    parser.add_argument('-anno_map', dest="map_out_fp", help="Path to map output")
    args = parser.parse_args()

    # Filter annotations
    annotations = pandas.io.parsers.read_csv(args.input_fp, sep="\t")
    resfams = annotations[annotations.database=='ResFam']
    top_resfams = resfams.groupby(['contig','start']).first().reset_index()
    top_resfams.to_csv(args.output_fp, sep="\t")

    # Create map
    top_resfams['contig_num'] = top_resfams['contig_num'].astype(str)
    top_resfams['start'] = top_resfams['start'].astype(str)
    top_resfams['stop'] = top_resfams['stop'].astype(str)
    
    top_resfams['gene_id'] = top_resfams['sample'] + '.' + top_resfams['contig_num'] + ':' + top_resfams['start'] + '-' + top_resfams['stop']
    top_resfams.to_csv(args.map_out_fp, sep="\t", cols=["gene_id", "anno_id", "description"], header=False, index=False)


if __name__ == "__main__":
    main()
