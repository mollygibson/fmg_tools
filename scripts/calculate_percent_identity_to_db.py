#!/usr/bin/env python                                                                                                                                                                                                                        
from __future__ import print_function

__author__ = "Molly Gibson"
__email__ = "molly.gibson@wustl.edu"

# Python imports                                                                                                                                                                                                                            
import argparse, subprocess, os, itertools, operator, sys
from Bio import SeqIO
from Bio import Entrez
import pandas

# Other imports                                                                                                                                                                                                                              
import helper as h

def main():
    parser = argparse.ArgumentParser(description="This script collapses all redundant contigs \
                        by percent identity, can be done within either library or antibiotic.")
    # arguments                                                                                                                                                                                                                             
    parser.add_argument('-prot', dest="protein_fp", help="Path to proteins")
    parser.add_argument('-prot_db', dest="protdb_fp", help="Path to database")
    parser.add_argument('-nucl', dest="nucl_fp", help="Path to nucleotides")
    parser.add_argument('-nucl_db', dest="nucldb_fp", help="Path to nucleotide database")
    parser.add_argument('-anno_tab', dest="anno_tab", help="Path to annotation file (tab)")
    parser.add_argument('-o', dest="output_fp", help="Path to output")
    parser.add_argument('-f', dest="force_flag", help="Force new output", action='store_true')
    args = parser.parse_args()

    # Make output dir
    if not os.path.exists(args.output_fp):
        h.run_command("mkdir " + args.output_fp)

    # BLAST Proteins
    if (args.protein_fp and args.protdb_fp) and (not os.path.exists(args.output_fp + "/blastp_results.txt") or args.force_flag):
        blast = "blastp -query " + args.protein_fp + " -db " + args.protdb_fp + " -out " + args.output_fp + "/blastp_results.txt -outfmt 6 -max_target_seqs 1"
        h.run_command(blast)

    if  (os.path.exists(args.output_fp + "/blastp_results.txt") and not os.path.exists(args.output_fp + "/final_protein_identity.txt")) or args.force_flag:
        global_identity = calculate_global_identity_blastp(args)

    if (args.nucl_fp and args.nucldb_fp) and (not os.path.exists(args.output_fp + "/blastn_results.txt") or args.force_flag):
        blast = "blastp -query " + args.protein_fp + " -db " + args.protdb_fp + " -out " + args.output_fp + "/blastn_results.txt -outfmt 6 -max_target_seqs 1"
        h.run_command(blast)

def calculate_global_identity_blastp(args):
    # set NCBI information
    Entrez.email = "molly.gibson@wustl.edu"

    # import proteins as fasta
    protein_dict = SeqIO.to_dict(SeqIO.parse(args.protein_fp, "fasta"))

    final = open(args.output_fp + "/final_protein_identity.txt", 'w')
    for result in open(args.output_fp + "/blastp_results.txt", 'r'):
        try:
            SeqIO.write(protein_dict[result.split()[0]], args.output_fp + "/temp_fasta1.fa", "fasta")
            id_number1 = result.split()[0]
            
            handle = Entrez.efetch(db='protein', id=result.split()[1].split('|')[3], rettype="text", retmode="fasta")
            id_number2 = result.split()[1].split('|')[3]
            output = open(args.output_fp + "/temp_fasta2.fa", 'w')
            output.write(handle.read())
            output.close()
        except:
            print("(1) Didn't work for " + id_number1 + " " + id_number2, file=sys.stderr)
            
        needle = "needle -outfile=" + args.output_fp + "/temp_identity.txt -asequence " + args.output_fp + "/temp_fasta1.fa -bsequence " + args.output_fp + "/temp_fasta2.fa -gapopen=10 -gapextend=0.5"
        h.run_command(needle)
        
        identity = "NA"
        for line in open(args.output_fp + "/temp_identity.txt", 'r'):
            if line.startswith("# Identity:"):
                identity = line.split("(")[1].rstrip("\n").rstrip(")").rstrip("%")
                    
        try:
            handle = Entrez.efetch(db='protein', id=result.split()[1].split('|')[3], rettype="text", retmode="gb")
            record = SeqIO.read(handle, "genbank")
            
            if "source" in record.annotations and "taxonomy" in record.annotations:
                final.write(id_number1 + "\t" + retrieve_annotations(id_number1, args.anno_tab) + "\t" + id_number2 + "\t" + identity + "\t" + record.description + "\t" + record.annotations["source"] + "\t" + ";".join(record.annotations["taxonomy"]) + "\n")
            elif "source" in record.annotations:
                final.write(id_number1 + "\t" + retrieve_annotations(id_number1, args.anno_tab) + "\t"+ id_number2 + "\t" + identity + "\t" + record.description + "\t" + record.annotations["source"]  + "\n")
            else:
                final.write(id_number1 + "\t" + retrieve_annotations(id_number1, args.anno_tab) + "\t"+ id_number2 + "\t" + identity + "\t" + record.description + "\t" + "\n")
        except:
            final.write(id_number1 + "\t" + retrieve_annotations(id_number1, args.anno_tab) + "\t"+ id_number2 + "\t" + identity + "\n")
            print("(2) Didn't work for " + id_number1 + " " + id_number2, file=sys.stderr)
        
    h.run_command("rm " + args.output_fp +  "/temp_*")

def retrieve_annotations(gene, anno_fp):
    gene_id = gene.split('.')[0].strip('>')
    contig = int(gene.split('.')[1].split(':')[0])
    gene_start = int(gene.split('.')[1].split(':')[1].split('-')[0])
    gene_stop = int(gene.split('.')[1].split(':')[1].split('-')[1])

    annotations = pandas.io.parsers.read_csv(anno_fp, sep="\t")

    gene_annotation = annotations[(annotations.index==gene_id)&(annotations.contig_num==contig)&(annotations.start==gene_start)&(annotations.stop==gene_stop)&(annotations.database=='ResFam')]
    return gene_annotation[:1][["anno_id","description"]].values[0][0] + "\t" + gene_annotation[:1][["anno_id","description"]].values[0][1]

if __name__ == "__main__":
    main()
