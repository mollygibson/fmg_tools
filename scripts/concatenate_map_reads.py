#!/usr/bin/env python                                                                                                                                                                                                                                                                                                                                                       

__author__ = "Molly Gibson"
__email__ = "molly.gibson@wustl.edu"


# Python imports                                                                                                                                                                                                                                                                                                                                                            
import argparse, subprocess, os, itertools, operator, pandas
# Other imports                                                                                                                                                                                                                                                                                                                                                             
import parse_config, parse_mapping

def main():
    parser = argparse.ArgumentParser(description="This script will concatenate all stats from read mapping using map_reads.py.")

    # arguments                                                                                                                                                                                                                                                                                                                                                             
    parser.add_argument('-map_fp', dest="map_fp", help="Path to directory where map_reads.py was run")
    parser.add_argument('-sample_ids', dest="sample_ids", help="Concatenate these samples")
    parser.add_argument('-annotations', dest="annotations_fp", help=".tab annotations file")
    parser.add_argument('-o', dest="output", type=argparse.FileType('w'),  help="Output file path")
    args = parser.parse_args()

    if args.sample_ids:
        samples_to_keep = [line.rstrip() for line in open(args.sample_ids)]

    for file_fp in os.listdir(args.map_fp):

        if os.path.isdir(args.map_fp.rstrip('/') + "/" + file_fp) and ((args.sample_ids and file_fp in samples_to_keep) or not args.sample_ids):
            sample_name = file_fp

            # get total reads                                                                                                                                     
            unpaired_stats = open(args.map_fp.rstrip("/") + "/" + sample_name + "/unpaired_alignment_stats.txt", 'r')
            unpaired_count = unpaired_stats.readline().split()[0]
            paired_stats = open(args.map_fp.rstrip("/") + "/" + sample_name  + "/paired_alignment_stats.txt", 'r')
            paired_count = paired_stats.readline().split()[0]
            total_reads = int(unpaired_count) + int(paired_count)

            try:
                for line in open(args.map_fp.rstrip('/') + "/" + file_fp + "/rpkm.txt"):
                    if not line.startswith("orf"):
                        gene_id = line.split("\t")[0].rstrip()

                        anno_desc = find_res_genes_annotation(args.annotations_fp, gene_id)
                        
                        args.output.write(sample_name + "\t" + line.rstrip() + "\t" + str(total_reads) + "\t" + anno_desc + "\n")
#                        print sample_name + "\t" + line.rstrip() + "\t" + str(total_reads) + "\t" + anno_desc
            except Exception, e:
                print e


def find_res_genes_annotation(annotations, gene):
    index = gene.split('.')[0]
    contig_num = gene.split('.')[1].split(':')[0]
    start = gene.split(':')[1].split('-')[0]
    stop = gene.split(':')[1].split('-')[1]

    anno_df = pandas.io.parsers.read_csv(annotations, sep="\t")
    anno_df.set_index(['sample', 'contig_num', 'start', 'stop'], inplace=True)
    anno_df.sortlevel(inplace=True)

    filtered_anno_df = anno_df.ix[(index, int(contig_num), int(start), int(stop))]
    res_annos = filtered_anno_df[filtered_anno_df['database'] == 'ResFam']

    return res_annos.iloc[0,3]    

if __name__ == "__main__":
    main()
