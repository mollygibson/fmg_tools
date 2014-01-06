#!/usr/bin/env python                                                                                                                                                          

__author__ = "Molly Gibson"
__email__ = "molly.gibson@wustl.edu"


# Python imports
import argparse
import subprocess
import os

# fmg_script imports
import annotate_contigs
import parse_config

parser = argparse.ArgumentParser(description="This script analyzes output from \
antibiotic functional metagenomic selections assembled using PARFuMS.")

# Arguments
parser.add_argument('-contigs', dest="contig_fp", help="Path to contig file from PARFuMS output")
parser.add_argument('-o', dest="output_fp", help="Path to output directory for analysis")

args = parser.parse_args()

# Call ORFs and Annotate Proteins using Pfam, TIGRFam, and Resfams
prefix = os.path.splitext(args.contig_fp)[0]
protein_fp = annotate_contigs.call_orfs(args.contig_fp, prefix)
annotation_fp = annotate_contigs.annotate_proteins(protein_fp, prefix)
