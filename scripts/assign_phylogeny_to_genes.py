#!/usr/bin/env python
__author__ = "Molly Gibson"
__email__ = "molly.gibson@wustl.edu"


# Python imports
import argparse, subprocess, os, itertools, operator
# Other imports
import pandas as pd
import helper as h
import parse_mapping

def main():
    parser = argparse.ArgumentParser(description="This script takes a fasta file that comes out of the annotations script, the blast results from phymmbl and add the phylogeny to the fasta")
    # arguments  
    parser.add_argument('-i', dest="fasta_fp", help="Path to final fasta from the annotation script")
    parser.add_argument('-o', dest="output_fp", help="Path to output directory")
    parser.add_argument('-m', dest="mapping_fp", help="Mapping file")
    parser.add_argument('-phymmbl', dest="phymmbl_fp", help="Path to phymmbl blast output")
    args = parser.parse_args()

    if not os.path.exists(args.output_fp):
        h.run_command("mkdir " + args.output_fp)

    phymmbl_blast = pd.io.parsers.read_table(args.phymmbl_fp)

    query_to_best_match = phymmbl_blast.set_index('QUERY_ID')['BEST_MATCH'].to_dict()
    query_to_score = phymmbl_blast.set_index('QUERY_ID')['SCORE'].to_dict()
    query_to_g = phymmbl_blast.set_index('QUERY_ID')['GENUS'].to_dict()
    query_to_f = phymmbl_blast.set_index('QUERY_ID')['FAMILY'].to_dict()
    query_to_o = phymmbl_blast.set_index('QUERY_ID')['ORDER'].to_dict()
    query_to_c = phymmbl_blast.set_index('QUERY_ID')['CLASS'].to_dict()
    query_to_p = phymmbl_blast.set_index('QUERY_ID')['PHYLUM'].to_dict()

    output = open(args.output_fp.strip("/") + "/" + os.path.splitext(os.path.basename(args.fasta_fp))[0] + "_with_tax.fasta", 'w')
    output_bh = open(args.output_fp.strip("/") + "/" + os.path.splitext(os.path.basename(args.fasta_fp))[0] + "_best_hit_summary.txt", 'w')
    output_bg = open(args.output_fp.strip("/") + "/" + os.path.splitext(os.path.basename(args.fasta_fp))[0] + "_best_genus_summary.txt", 'w')
    output_bh_s = open(args.output_fp.strip("/") + "/" + os.path.splitext(os.path.basename(args.fasta_fp))[0] + "_best_hit_summary_by_sample.txt", 'w')
    output_bg_s = open(args.output_fp.strip("/") + "/" + os.path.splitext(os.path.basename(args.fasta_fp))[0] + "_best_genus_summary_by_sample.txt", 'w')
    output_nh = open(args.output_fp.strip("/") + "/" + os.path.splitext(os.path.basename(args.fasta_fp))[0] + "_no_hit.txt", 'w')

    best_hit = {}
    best_g = {}
    best_hit_by_sample = {}
    best_g_by_sample = {}
    for line in open(args.fasta_fp, 'r'):
        if line.startswith(">"):
            sample_id = line.split()[0].split('.')[0].strip('>')
            contig = line.split()[0].split('.')[1].split(':')[0]
            query = sample_id + "_" + contig
            
            indv_id = parse_mapping.main(args.mapping_fp).id[sample_id]

            if query in query_to_best_match and isinstance(query_to_g[query],basestring):
                output.write(line.rstrip() + " best_hit:" + query_to_best_match[query] + " taxonomy:p__" + query_to_p[query] + ";c__" + query_to_c[query] + ";o__" + query_to_o[query]+ ";f__" + query_to_f[query]+ ";g__" + query_to_g[query] + "\n")
                
                if indv_id in best_hit_by_sample:
                    if query_to_best_match[query] in best_hit_by_sample[indv_id]:
                        best_hit_by_sample[indv_id][query_to_best_match[query]] = best_hit_by_sample[indv_id][query_to_best_match[query]] + 1
                    else:
                        best_hit_by_sample[indv_id][query_to_best_match[query]] = 1
                else:
                    best_hit_by_sample[indv_id] = {}
                    best_hit_by_sample[indv_id][query_to_best_match[query]] = 1

                if indv_id in best_g_by_sample:
                    if query_to_g[query] in best_g_by_sample[indv_id]:
                        best_g_by_sample[indv_id][query_to_g[query]] = best_g_by_sample[indv_id][query_to_g[query]] + 1
                    else:
                        best_g_by_sample[indv_id][query_to_g[query]] = 1
                else:
                    best_g_by_sample[indv_id] = {}
                    best_g_by_sample[indv_id][query_to_g[query]] = 1

                if query_to_best_match[query] in best_hit:
                    best_hit[query_to_best_match[query]] = best_hit[query_to_best_match[query]] + 1
                else:
                    best_hit[query_to_best_match[query]] = 1

                if query_to_g[query] in best_g:
                    best_g[query_to_g[query]] = best_g[query_to_g[query]] + 1
                else:
                    best_g[query_to_g[query]] = 1

            else:
                output_nh.write(line.rstrip() + "\n")
                output.write(line.rstrip() + " best_hit:none taxonomy:none\n")
        else:
            output.write(line.rstrip() + "\n")

    for hit in best_hit:
        output_bh.write(hit + str(best_hit[hit]) + "\n")
    
    for g in best_g:
        output_bg.write(g + str(best_g[g]) + "\n")

    for sample_id in best_hit_by_sample:
        for hit in best_hit_by_sample[sample_id]:
            output_bh_s.write(sample_id + "\t" + hit + "\t" + str(best_hit_by_sample[sample_id][hit]) + "\t" + str(sum(best_hit_by_sample[sample_id].values())) + "\t" + str(best_hit_by_sample[sample_id][hit]/float(sum(best_hit_by_sample[sample_id].values()))) + "\n")

    for sample_id in best_g_by_sample:
        for g in best_g_by_sample[sample_id]:
            output_bg_s.write(sample_id + "\t" + g + "\t" + str(best_g_by_sample[sample_id][g]) + "\t" + str(sum(best_g_by_sample[sample_id].values())) + "\t" + str(best_g_by_sample[sample_id][g]/float(sum(best_g_by_sample[sample_id].values()))) +"\n")

if __name__ == "__main__":
    main()
