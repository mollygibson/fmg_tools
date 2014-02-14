#!/usr/bin/env python                                                                                                                                                          

__author__ = "Molly Gibson"
__email__ = "molly.gibson@wustl.edu"


# Python imports
import argparse, subprocess, os, itertools, operator
# Other imports
import parse_config, parse_mapping

def main():
    parser = argparse.ArgumentParser(description="This script collapses all redundant proteins \
                        or ORFs by percent identity, can be done within either library or antibiotic.")
    # arguments
    parser.add_argument('-prot', dest="protein_fp", help="Path to updated proteins from annotation output")
    parser.add_argument('-nucl', dest="nucleotide_fp", help="Path to updated nucleotides from annotation output")
    parser.add_argument('-id', dest="percent_id", help="Percent identity to collapse genes on", default="1.0")
    parser.add_argument('--lib', dest="within_lib", help="Collapse within library", action="store_true", default=False)
    parser.add_argument('--abx', dest="within_abx", help="Collapse within antibiotic", action="store_true", default=False)
    parser.add_argument('-o', dest="output_fp", help="Path to output directory")
    parser.add_argument('-f', dest="override", help="Force override otuput directory", action="store_true", default=False)
    args = parser.parse_args()
    
    # check to see that either nucleotides or proteins are provided
    if not args.nucleotide_fp and not args.protein_fp:
        parser.exit(status=0, message="You must provide either nucleotides or proteins to collapse. \n Check usage with 'collapse_redundant_genes.py -h'.\n\n")

    if args.within_lib and args.within_abx:
        parser.exit(status=0, message="You can only collapse on either library OR antibiotic selection. If you would like to collapse the entire \
                          functional selection, don't provide either flat. \n Check usage with 'collapse_redundant_genes.py -h'.\n\n")

    # Figure out the output directory if it isn't given
    if not args.output_fp:
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

    # Make output directory
    if os.path.isdir(output_fp):
        if args.override:
            subprocess.call('rm -r ' + output_fp, shell=True)
            subprocess.call('mkdir ' + output_fp, shell=True)
        else:
            parser.exit(status=0, message="If you want to override the output directory use the -f flag. \n Check usage with \
                          'collapse_redundant_genes.py -h'. \n\n")
    else:
        subprocess.call('mkdir ' + output_fp, shell=True)

    # Split the files
    if args.within_lib:
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
            file_out = open(output_fp + "/protiens_updated.faa", 'a')
            for line in open(args.protein_fp, 'r'):
                file_out.write(line)
        else:
            file_out = open(output_fp + "/proteins_updated.fna", 'a')
            for line in open(args.nucleotide_fp, 'r'):
                file_out.write(line)

    # Run clustering
    for fasta in os.listdir(output_fp):
        if args.protein_fp:
            command = 'cd-hit -i ' + output_fp + '/' + fasta + ' -o ' + output_fp + '/' + fasta.split('.')[0] + '_unique.faa -c ' + args.percent_id + ' -aS 1.0 -g 1 -d 0'
        elif args.nucleotide_fp:
            command = 'cd-hit-est -i ' + output_fp + '/' + fasta + ' -o ' + output_fp + '/' + fasta.split('.')[0] + '_unique.fna -c ' + args.percent_id + ' -aS 1.0 -g 1 -d 0 -r 1'
        subprocess.call(command, shell=True)

    # Concatenate and get annotation files
    if args.within_lib:
        if args.nucleotide_fp:
            command = 'cat ' + output_fp + '/*_unique.fna > ' + output_fp + '/unique_within_lib_' + args.percent_id + '.fna'
        else:
            command = 'cat ' + output_fp + '/*_unique.faa > ' + output_fp + '/unique_within_lib_'+ args.percent_id + '.faa'
    elif  args.within_abx:
        if args.nucleotide_fp:
            command = 'cat ' + output_fp + '/*_unique.fna > ' + output_fp + '/unique_within_abx_'+ args.percent_id + '.fna'
        else:
            command = 'cat ' + output_fp + '/*_unique.faa > ' + output_fp + '/unique_within_abx_'+ args.percent_id + '.faa'
    else:
        if args.nucleotide_fp:
            command = 'cat ' + output_fp + '/*_unique.fna > ' + output_fp + '/all_unique_'+ args.percent_id + '.fna'
        else:
            command = 'cat ' + output_fp + '/*_unique.faa > ' + output_fp + '/all_unique_'+ args.percent_id + '.faa'

    subprocess.call(command, shell=True)

if __name__ == "__main__":
    main()



