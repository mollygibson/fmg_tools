#!/usr/bin/env python                                                                                                                                                                                                                        
__author__ = "Molly Gibson"
__email__ = "molly.gibson@wustl.edu"

# Python imports                                                                                                                                                                                                                            
import argparse, subprocess, os, itertools, operator
from Bio import SeqIO
from Bio import Entrez
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
    parser.add_argument('-o', dest="output_fp", help="Path to output")
    args = parser.parse_args()

    # Make output dir
    if not os.path.exists(args.output_fp):
        h.run_command("mkdir " + args.output_fp)

    # BLAST Proteins
    if (args.protein_fp and args.protdb_fp) and not os.path.exists(args.output_fp + "/blastp_results.txt"):
        blast = "blastp -query " + args.protein_fp + " -db " + args.protdb_fp + " -out " + args.output_fp + "/blastp_results.txt -outfmt 6 -max_target_seqs 1"
        h.run_command(blast)

    if os.path.exists(args.output_fp + "/blastp_results.txt"):
        global_identity = calculate_global_identity_blastp(args)

    if args.nucl_fp and args.nucldb_fp:
        blast = "blastp -query " + args.protein_fp + " -db " + args.protdb_fp + " -out " + args.output_fp + "/blastn_results.txt -outfmt 6 -max_target_seqs 1"
        h.run_command(blast)

def calculate_global_identity_blastp(args):
    # set NCBI information
    Entrez.email = "molly.gibson@wustl.edu"

    # import proteins as fasta
    protein_dict = SeqIO.to_dict(SeqIO.parse(args.protein_fp, "fasta"))

    final = open(args.output_fp + "/final_protein_identity.txt", 'w')
    for result in open(args.output_fp + "/save_blastp.txt", 'r'):
        try:
            SeqIO.write(protein_dict[result.split()[0]], args.output_fp + "/temp_fasta1.fa", "fasta")
            id_number1 = result.split()[0]
            
            handle = Entrez.efetch(db='protein', id=result.split()[1].split('|')[3], rettype="text", retmode="fasta")
            id_number2 = result.split()[1].split('|')[3]
            output = open(args.output_fp + "/temp_fasta2.fa", 'w')
            output.write(handle.read())
            output.close()
            
            needle = "needle -outfile=" + args.output_fp + "/temp_identity.txt -asequence " + args.output_fp + "/temp_fasta1.fa -bsequence " + args.output_fp + "/temp_fasta2.fa -gapopen=10 -gapextend=0.5"
            h.run_command(needle)
            
            identity = "NA"
            for line in open(args.output_fp + "/temp_identity.txt", 'r'):
                if line.startswith("# Identity:"):
                    identity = line.split("(")[1].rstrip("\n").rstrip(")").rstrip("%")
                    
            final.write(id_number1 + "\t" + id_number2 + "\t" + identity + "\n")
        except:
            print "Didn't work for " + id_number1 + " " + id_number2
        

if __name__ == "__main__":
    main()
