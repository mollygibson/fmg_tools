#!/usr/bin/env python                                                                                                                                                               

__author__ = "Molly Gibson"
__email__ = "molly.gibson@wustl.edu"


# Python imports                                                                                                                                                                   
import argparse, subprocess, os, itertools, operator
import pandas
import numpy as np
# Other imports                                                                                                                                                                    
import parse_config, parse_mapping

def main():
    parser = argparse.ArgumentParser(description="This script summarizes output from PARFuMS annotations")
 
    # arguments                                                                                                                                                                    
    parser.add_argument('-annotations', dest="annotation_fp", help="Path to annotation file from functional metagenomic selections")
    parser.add_argument('-o', dest="output_fp", help="Path to output directory")
    args = parser.parse_args()

    if not os.path.isdir(args.output_fp):
        subprocess.call('mkdir ' + args.output_fp, shell=True)    

    # import annotations and extract library and antibiotic information
    annotations = pandas.io.parsers.read_table(args.annotation_fp)
    extraction = annotations['contig'].str.extract('.*ID:(?P<library>[^ ]*) Contig:(?P<contig_num>\d*) .* Len:(?P<contig_len>\d*).*abx:(?P<antibiotic>[^ ]*).*')
    annotations = pandas.merge(annotations,extraction,on=annotations.index,how='outer')
    annotations['gene_len'] = annotations['stop'].astype(float) - annotations['start'].astype(float)


    # analyze contig information
    [contigs_by_library, contigs_by_abx, contigs_by_lib_abx, contig_dist] = count_unique_contigs(annotations)
    contigs_by_library.to_csv(args.output_fp + "/contigs_by_library.txt", sep="\t")
    contigs_by_abx.to_csv(args.output_fp + "/contigs_by_abx.txt", sep="\t")
    contigs_by_lib_abx.to_csv(args.output_fp + "/contigs_by_lib_abx.txt", sep="\t")
    contig_dist.to_csv(args.output_fp + "/contig_dist.txt", sep="\t")

    # analyze gene information (all greater than 300bp)
    [gene_by_lib, gene_by_abx,gene_by_lib_abx, res_gene_by_lib, res_gene_by_abx, res_gene_by_lib_abx] = count_unique_genes(annotations)
    gene_by_lib.to_csv(args.output_fp + "/genes_by_lib.txt", sep="\t")
    gene_by_abx.to_csv(args.output_fp + "/genes_by_abx.txt", sep="\t")
    gene_by_lib_abx.to_csv(args.output_fp + "/genes_by_lib_abx.txt", sep="\t")
    res_gene_by_lib.to_csv(args.output_fp + "/res_gene_by_lib.txt", sep="\t")
    res_gene_by_abx.to_csv(args.output_fp + "/res_gene_by_abx.txt", sep="\t")
    res_gene_by_lib_abx.to_csv(args.output_fp + "/res_gene_by_lib_abx.txt", sep="\t")
    
   # make annotation table
    annotation_table = make_annotation_table(annotations)
    annotation_table.to_csv(args.output_fp + "/annotation_table.txt", sep="\t")

    # analyze co-localized resistance genes
    num_res_genes_per_contig = analyze_colocalized_resistance(annotations)
    num_res_genes_per_contig.to_csv(args.output_fp + "/num_res_genes_per_contig.txt", sep="\t")

    # analyze mobilization elements
    mobile_element_key_words = ['transposase', 'transposon', 'conjugative', 'integrase', 'integron', 'recombinase', 'conjugal', 'mobilization', 'recombination', 'plasmid']
    annotations['mobilization'] = annotations.apply(lambda row: is_mobile(row['description'], mobile_element_key_words), axis=1)
    res_and_mobile_per_contig = analyze_mobilization_elements(annotations)
    res_and_mobile_per_contig.to_csv(args.output_fp + "/res_and_mobile_per_contig.txt", sep="\t")


#------------------------------------------------------------
def is_mobile(description, key_words):
    mobile=False
    for word in key_words:
        if word in description.lower():
            mobile=True
    return mobile

def analyze_mobilization_elements(annotations):
    top_hits = annotations[annotations.gene_len>300].groupby(['contig', 'start', 'stop'], as_index=False).first()

#    res_or_mobil_top_hits = top_hits[((top_hits.database == "ResFam") | (top_hits.mobilization == True))]

    res_top_hits = top_hits[top_hits.database == "ResFam"]
    mob_top_hits = top_hits[top_hits.mobilization == True]

    num_res_genes_per_contig = pandas.DataFrame(res_top_hits.groupby('contig').size())
    num_res_genes_per_contig.columns = ['num_res']
    num_mob_genes_per_contig = pandas.DataFrame(mob_top_hits.groupby('contig').size())
    num_mob_genes_per_contig.columns = ['num_mob']

    res_or_mobil_top_hits = num_res_genes_per_contig.join(num_mob_genes_per_contig, on=num_res_genes_per_contig.index, how='outer').fillna(0)
    res_and_mobil_top_hits = res_or_mobil_top_hits[((res_or_mobil_top_hits.num_mob>0) & (res_or_mobil_top_hits.num_res>0))]

    return res_and_mobil_top_hits



def analyze_colocalized_resistance(annotations):
    # Top hits and only resistance hits
    top_hits = annotations[annotations.gene_len>300].groupby(['contig', 'start', 'stop'], as_index=False).first()
    res_top_hits = top_hits[top_hits.database == "ResFam"]

    num_res_genes_per_contig = pandas.DataFrame(res_top_hits.groupby('contig').size())
    return num_res_genes_per_contig

def make_annotation_table(annotations):
    top_hits = annotations[annotations.gene_len>300].groupby(['contig', 'start', 'stop'], as_index=False).first()

    anno_counts = top_hits.groupby(['library', 'anno_id'], as_index=False).size().reset_index()
    anno_counts.columns = ['library', 'anno_id', 'count']

    anno_table = anno_counts.pivot_table('count', 'library', 'anno_id').fillna(0)
    return anno_table

def count_unique_genes(annotations):
    # filter to top hits
    top_hits = annotations.groupby(['contig', 'start', 'stop'], as_index=False).first()

    unique_by_lib = pandas.DataFrame(top_hits.groupby('library').size())
    unique_by_abx = pandas.DataFrame(top_hits.groupby('antibiotic').size())
    unique_by_lib_abx = pandas.DataFrame(top_hits.groupby(['library','antibiotic']).size())

    res_genes = annotations[annotations.database == "ResFam"]
    res_top_hits = res_genes.groupby(['contig', 'start', 'stop'], as_index=False).first()

    res_unique_by_lib = pandas.DataFrame(res_top_hits.groupby('library').size())
    res_unique_by_abx = pandas.DataFrame(res_top_hits.groupby('antibiotic').size())
    res_unique_by_lib_abx = pandas.DataFrame(res_top_hits.groupby(['library','antibiotic']).size())
    
    return [unique_by_lib, unique_by_abx, unique_by_lib_abx, res_unique_by_lib, res_unique_by_abx, res_unique_by_lib_abx]


def count_unique_contigs(annotations):
    # unique contigs by library
    unique_by_lib = pandas.DataFrame(annotations.groupby('library').contig.nunique())

    annotations['contig_len'] = annotations['contig_len'].astype(float)
    large_contigs = annotations[annotations.contig_len>500]
    unique_by_lib_large = pandas.DataFrame(large_contigs.groupby('library').contig.nunique())

    contigs_by_library = pandas.merge(unique_by_lib,unique_by_lib_large, on=unique_by_lib.index, how='outer')
    contigs_by_library.columns= ['library', 'total', 'total_greater_500']

    # unique contigs by library and antibiotic
    unique_by_lib_abx = pandas.DataFrame(annotations.groupby(['library','antibiotic']).contig.nunique())
    unique_by_lib_abx = unique_by_lib_abx.reset_index()

    annotations['contig_len'] = annotations['contig_len'].astype(float)
    large_contigs = annotations[annotations.contig_len>500]
    unique_by_lib_abx_large = pandas.DataFrame(large_contigs.groupby(['library','antibiotic']).contig.nunique())
    unique_by_lib_abx_large = unique_by_lib_abx_large.reset_index()

    contigs_by_lib_abx = pandas.merge(unique_by_lib_abx,unique_by_lib_abx_large, on=unique_by_lib_abx.index, how='outer')
    contigs_by_lib_abx = contigs_by_lib_abx[['library_x', 'antibiotic_x', '0_x', '0_y']]
    contigs_by_lib_abx.columns= ['library', 'antibiotic','total', 'total_greater_500']

    # unique contigs by antibiotic
    unique_by_abx = pandas.DataFrame(annotations.groupby('antibiotic').contig.nunique())

    annotations['contig_len'] =annotations['contig_len'].astype(float)
    large_contigs = annotations[annotations.contig_len>500]
    unique_by_abx_large = pandas.DataFrame(large_contigs.groupby('antibiotic').contig.nunique())

    contigs_by_abx = pandas.merge(unique_by_abx,unique_by_abx_large, on=unique_by_abx.index, how='outer')
    contigs_by_abx.columns= ['antibitoic', 'total', 'total_greater_500']    
    
    grouped  = annotations.groupby('library')
    contig_stats = grouped['contig_len'].agg([np.mean, np.std])

    return [contigs_by_library, contigs_by_abx, contigs_by_lib_abx, contig_stats]




if __name__ == "__main__":
    main()

