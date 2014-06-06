#!/usr/bin/env python                                                                                                                                                          

__author__ = "Molly Gibson"
__email__ = "molly.gibson@wustl.edu"


# Python imports
import argparse, subprocess, os, itertools, operator
# Other imports
from Bio import SeqIO
import parse_config, parse_mapping
import helper as h

def main():
    parser = argparse.ArgumentParser(description="This script collapses all redundant proteins \
                        or ORFs by percent identity, can be done within either library or antibiotic.")
    # arguments
    parser.add_argument('-prot', dest="protein_fp", help="Path to updated proteins from annotation output")
    parser.add_argument('-nucl', dest="nucleotide_fp", help="Path to updated nucleotides from annotation output")
    parser.add_argument('-id', dest="percent_id", help="Percent identity to collapse genes on", default="1.0")
    parser.add_argument('-m', dest="mapping_fp", help="Mapping file path")
    parser.add_argument('--lib', dest="within_lib", help="Collapse within library", action="store_true", default=False)
    parser.add_argument('--abx', dest="within_abx", help="Collapse within antibiotic", action="store_true", default=False)
    parser.add_argument('-o', dest="output_fp", help="Path to output directory")
    parser.add_argument('-f', dest="override", help="Force override otuput directory", action="store_true", default=False)
    args = parser.parse_args()
    
    # check to see that either nucleotides or proteins are provided
    if not args.nucleotide_fp and not args.protein_fp:
        parser.exit(status=0, message="You must provide either nucleotides or proteins to collapse. \n Check usage with 'cluster_genes.py -h'.\n\n")

    if (args.within_lib or args.within_abx) and not args.mapping_fp:
        parser.exit(status=0, meassage="You need a mapping file if collaping within library or antibiotic. \n Check usage with 'cluster_functional_genes.py -h. \n\n")

    # determine the output directory if it isn't given
    if not args.output_fp:
        if args.within_lib and args.within_abx:
            if args.nucleotide_fp:
                output_fp = 'collapsed_nuc_lib-abx_' + str(args.percent_id)
            else:
                output_fp = 'collapse_prot_lib-abx_' + str(args.percent_id)
        if args.within_lib:
            if args.nucleotide_fp:
                output_fp = 'collapsed_nuc_lib_' + str(args.percent_id)
            else:
                output_fp = 'collapse_prot_lib_' + str(args.percent_id)
        elif  args.within_abx:
            if args.nucleotide_fp:
                output_fp = 'collapsed_nuc_abx_' + str(args.percent_id)
            else:
                output_fp = 'collapse_prot_abx_' + str(args.percent_id)
        else:
            if args.nucleotide_fp:
                output_fp = 'collapsed_nuc_' + str(args.percent_id)
            else:
                output_fp = 'collapse_prot_' + str(args.percent_id)
    else:
        output_fp = args.output_fp

    # make output directory
    if os.path.isdir(output_fp):
        if args.override:
            subprocess.call('rm -r ' + output_fp, shell=True)
            subprocess.call('mkdir ' + output_fp, shell=True)
        else:
            parser.exit(status=0, message="If you want to override the output directory use the -f flag. \n Check usage with \
                          'cluster_functional_genes.py -h'. \n\n")
    else:
        subprocess.call('mkdir ' + output_fp, shell=True)

    # split the files by library 
    split_fasta_file(output_fp, args)

    # run clustering
    cluster_fasta(output_fp, args)

    # create resistance profile for each gene and output new fasta
    if not args.within_lib and not args.within_abx and args.mapping_fp:
        create_resistance_profile(output_fp, args)

    summarize_output(output_fp, args)


# ----------------------------------------------------------------------------------------------
# methods
def split_fasta_file(output_fp, args):
    if args.within_lib and args.within_abx:
        if args.protein_fp:
            for line in open(args.protein_fp, 'r'):
                if line.startswith('>'):
                    library = line.split()[1].split(":")[1]
                    abx = line.split()[5].split(":")[1]
                    file_out = open(output_fp + "/" + library + "-" + abx + ".faa", 'a')
                file_out.write(line)
        else:
            for line in open(args.nucleotide_fp, 'r'):
                if line.startswith('>'):
                    library = line.split()[1].split(":")[1]
                    abx = line.split()[5].split(":")[1]
                    file_out = open(output_fp + "/" + library + "-" +  ".fna", 'a')
                file_out.write(line)
    elif args.within_lib:
        if args.protein_fp:
            for line in open(args.protein_fp, 'r'):
                if line.startswith('>'):
                    library = line.split()[1].split(":")[1]
                    file_out = open(output_fp + "/" + library + ".faa", 'a')
                file_out.write(line)
        else:
            for line in open(args.nucleotide_fp, 'r'):
                if line.startswith('>'):
                    library = line.split()[1].split(":")[1]
                    file_out = open(output_fp + "/" + library + ".fna", 'a')
                file_out.write(line)
    elif args.within_abx:
        if args.protein_fp:
            for line in open(args.protein_fp, 'r'):
                if line.startswith('>'):
                    abx = line.split()[5].split(":")[1]
                    file_out = open(output_fp + "/" + abx + ".faa", 'a')
                file_out.write(line)
        else:
            for line in open(args.nucleotide_fp, 'r'):
                if line.startswith('>'):
                    abx = line.split()[5].split(":")[1]
                    file_out = open(output_fp + "/" + abx + ".fna", 'a')
                file_out.write(line)
    else:
        if args.protein_fp:
            file_out = open(output_fp + "/all_proteins.faa", 'a')
            for line in open(args.protein_fp, 'r'):
                file_out.write(line)
        else:
            file_out = open(output_fp + "/all_nucleotides.fna", 'a')
            for line in open(args.nucleotide_fp, 'r'):
                file_out.write(line)
    file_out.close()

def cluster_fasta(output_fp, args):
    for fasta in os.listdir(output_fp):
        basename = os.path.splitext(os.path.basename(fasta))[0]
        if args.protein_fp:
            cluster_fp =  output_fp + '/' + basename + '_unique.faa'
            command = 'cd-hit -i ' + output_fp + '/' + fasta + ' -o ' + cluster_fp + ' -c ' + args.percent_id + ' -aS 1.0 -g 1 -d 0'
        elif args.nucleotide_fp:
            cluster_fp = output_fp + '/' + basename + '_unique.fna'
            command = 'cd-hit-est -i ' + output_fp + '/' + fasta + ' -o ' + cluster_fp + ' -c ' + args.percent_id + ' -aS 1.0 -g 1 -d 0 -r 1'
        h.run_command(command)

    if args.within_lib and args.within_abx:
        if args.nucleotide_fp:
            command = 'cat ' + output_fp + '/*_unique.fna > ' + output_fp + '/unique_within_lib-abx_' + args.percent_id + '.fna'
        else:
            command = 'cat ' + output_fp + '/*_unique.faa > ' + output_fp + '/unique_within_lib-abx_'+ args.percent_id + '.faa'
            h.run_command(command)
    elif args.within_lib:
        if args.nucleotide_fp:
            command = 'cat ' + output_fp + '/*_unique.fna > ' + output_fp + '/unique_within_lib_' + args.percent_id + '.fna'
        else:
            command = 'cat ' + output_fp + '/*_unique.faa > ' + output_fp + '/unique_within_lib_'+ args.percent_id + '.faa'
        h.run_command(command)
    elif  args.within_abx:
        if args.nucleotide_fp:
            command = 'cat ' + output_fp + '/*_unique.fna > ' + output_fp + '/unique_within_abx_'+ args.percent_id + '.fna'
        else:
            command = 'cat ' + output_fp + '/*_unique.faa > ' + output_fp + '/unique_within_abx_'+ args.percent_id + '.faa'
        h.run_command(command)
            

def create_resistance_profile(output_fp, args):
    cluster_to_abx = {}
    gene_to_cluster = {}

    if args.protein_fp:
        cluster_fp = output_fp + '/all_proteins_unique.faa'
        all_orfs_fp = output_fp + '/all_proteins.faa'
    else:
        cluster_fp = output_fp + '/all_nucleotides_unique.fna'
        all_orfs_fp = output_fp + '/all_nucleotides.fna'

    for gene in open(cluster_fp + ".bak.clstr", 'r'):
        if gene.rstrip().endswith('%'):
            [cluster_num, length, gene_id,trash, percent] = gene.rstrip().split()
        else:
            [cluster_num, length, gene_id, percent] = gene.rstrip().split()

        gene_to_cluster[gene_id.rstrip('.').strip('>')] = cluster_num
        abx = parse_mapping.main(args.mapping_fp).abx[gene_id.split('.')[0].strip('>')]
        if not cluster_num in cluster_to_abx.keys():
            cluster_to_abx[cluster_num] = []
        if not abx in cluster_to_abx[cluster_num]:
            cluster_to_abx[cluster_num].append(abx)

    res_profile = {}
    updated_records = []
    for record in SeqIO.parse(all_orfs_fp, 'fasta'):
        
        profile = ",".join(cluster_to_abx[gene_to_cluster[record.id]])
        res_profile[record.id] = profile

        record.description = record.description + " abx_profile:" + profile
        updated_records.append(record)

    if args.protein_fp:
        output_seqs = output_fp + "/all_proteins_abx_profile.fna"
        SeqIO.write(updated_records, output_seqs, 'fasta')
    else:
        output_seqs = output_fp + "/all_nucleotides_abx_profile.fna"
        SeqIO.write(updated_records, output_seqs, 'fasta')

    output_map = open(output_fp + "/abx_profile_map.txt", 'w')
    output_map.write("orf_name\tres_profile\n")
    for gene in res_profile.keys():
        output_map.write(gene + "\t" + res_profile[gene] + "\n")

def summarize_output(output_fp, args):
    output_file = open(output_fp + "/gene_summary.txt", 'w')
    for fasta_file in os.listdir(output_fp):
        if fasta_file.endswith("_unique.faa"):
            sequences = list(SeqIO.parse(output_fp + "/" + fasta_file, "fasta"))
            
            output_file.write(fasta_file + "\t" + str(len(sequences)) + "\n")



if __name__ == "__main__":
    main()



