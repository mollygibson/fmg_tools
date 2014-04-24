#!/usr/bin/env python                                                                                                                                                               

__author__ = "Molly Gibson"
__email__ = "molly.gibson@wustl.edu"


# Python imports                                                                                                                                                                   
import argparse, subprocess, os, itertools, operator
# Other imports                                                                                                                                                                    
import parse_config, parse_mapping

def main():
    parser = argparse.ArgumentParser(description="This script summarizes output from PARFuMS annotations")
 
    # arguments                                                                                                                                                                    
    parser.add_argument('-annotations', dest="annotation_fp", help="Path to annotation file from functional metagenomic selections")
    parser.add_argument('-min_prot_length', dest="min_plen", help="Minimum protein length to consider for functional annotation", default=0)
    args = parser.parse_args()


    resfam_ids = []
    function_by_library = {}
    for line in open(annotation_fp, 'r'):
        if line.startswith("ID:"):
            current_status = {}
            elements = line.rstrip().split()
            for item in elements:
                item_s = item.split(":")
                current_status[item_s[0]] = current_status[item_s[1]]

    
