#!/usr/bin/env python                                                                                                                                                          

__author__ = "Molly Gibson"
__email__ = "molly.gibson@wustl.edu"


# Python imports
import argparse, subprocess, os, itertools, operator
# Other imports
import parse_config, parse_mapping

def main():
    parser = argparse.ArgumentParser(description="This script collapses all redundant contigs \
                        by percent identity, can be done within either library or antibiotic.")
    # arguments
    parser.add_argument('-contigs', dest="contig_fp", help="Path to contigs from PARFuMS")
    parser.add_argument('-m', dest="mapping_fp", help="Path to mapping file")
    parser.add_argument('-id', dest="percent_id", help="Percent identity to collapse genes on", default="1.0")
    parser.add_argument('--lib', dest="within_lib", help="Collapse within library", action="store_true", default=False)
    parser.add_argument('--abx', dest="within_abx", help="Collapse within antibiotic", action="store_true", default=False)
    parser.add_argument('-o', dest="output_fp", help="Path to output directory")
    parser.add_argument('-f', dest="override", help="Force override otuput directory", action="store_true", default=False)
    parser.add_argument('-min', dest="min_len", help="Minimum length of contig to retain")
    parser.add_argument('-annotations', dest="annotations_fp", help="Filter out annotations that are redundant as well")
    args = parser.parse_args()
    
    # check to see that either nucleotides or proteins are provided
    if not args.contig_fp:
        parser.exit(status=0, message="You must provide a contig file to collapse. \n Check usage with 'collapse_redundant_contigs.py -h \n\n")

    if (args.within_lib or args.within_abx) and not args.mapping_fp:
        parser.exit(status=0, message="You must provide a mapping file to collapse within library or antibiotic. \n Check usage with 'collapse_redundant_contigs.py -h \n\n")

    if args.within_lib and args.within_abx:
        parser.exit(status=0, message="You can only collapse on either library OR antibiotic selection. If you would like to collapse the entire \
                          functional selection, don't provide either flat. \n Check usage with 'collapse_redundant_contigs.py -h'.\n\n")

    # Figure out the output directory if it isn't given
    if not args.output_fp:
        if args.within_lib:
            output_fp = 'collapse_contig_lib_' + str(args.percent_id)
        elif  args.within_abx:
            output_fp = 'collapse_contig_abx_' + str(args.percent_id)
        else:
            output_fp = 'collapse_contig_' + str(args.percent_id)
    else:
        output_fp = args.output_fp

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

    # set min_len
    if args.min_len:
        min_len = args.min_len
    else:
        min_len = 0

    # Split the files
    if args.within_lib:
        for line in open(args.contig_fp, 'r'):
            if line.startswith('>'):
                library = parse_mapping.main(args.mapping_fp).id[line.split("_Contig_")[0].strip(">")]
                file_out = open(output_fp + "/" + library + ".fna", 'a')

                sample_name = line.split("_Contig_")[0].strip(">").rstrip()
                contig_num = line.split("_Contig_")[1].split("_")[0].rstrip()
                if len(line.split("_Mean:")) > 1:
                    mean = line.split("_Mean:")[1].split("_")[0].rstrip()
                else:
                    mean = "NA"
                if len(line.split("_Len:")) > 1:
                    length = line.split("_Len:")[1].rstrip()
                else:
                    length = "NA"
                newHeader = ">" + sample_name + "_" + contig_num +" ID:" + parse_mapping.main(args.mapping_fp).id[sample_name] + " Contig:" + contig_num + " Mean:" + mean + " Len:" + length  + " abx:" + parse_mapping.main(args.mapping_fp).abx[sample_name] 
            else:
                if len(line) >= int(min_len):
                    file_out.write(newHeader + "\n")
                    file_out.write(line)
    elif args.within_abx:
        for line in open(args.contig_fp, 'r'):
            if line.startswith('>'):
                abx = parse_mapping.main(args.mapping_fp).abx[line.split("_Contig_")[0].strip(">")]
                file_out = open(output_fp + "/" + abx + ".fna", 'a')
                
                sample_name = line.split("_Contig_")[0].strip(">").rstrip()
                contig_num = line.split("_Contig_")[1].split("_")[0].rstrip()

                if len(line.split("_Mean:")) > 1:
                        mean = line.split("_Mean:")[1].split("_")[0].rstrip()
                else:
                    mean = "NA"
                if len(line.split("_Len:")) > 1:
                    length = line.split("_Len:")[1].rstrip()
                else:
                    length = "NA"
                newHeader = ">" + sample_name + "+" + contig_num + " ID:" + parse_mapping.main(args.mapping_fp).id[sample_name] + " Contig:" + contig_num + " Mean:" + mean + " Len:" + length  + " abx:" + parse_mapping.main(args.mapping_fp).abx[sample_name]

            else:
                if len(line) >= int(min_len):
                    file_out.write(newHeader + "\n")
                    file_out.write(line)
    else:
        file_out = open(output_fp + "/contigs_updated.fna", 'a')
        save_header=""
        for line in open(args.contig_fp, 'r'):
            if line.startswith('>'):
                save_header = line.rstrip()
            else:
                if len(line) >= int(min_len):
                    file_out.write(save_header + "\n")
                    file_out.write(line)
    file_out.close()

    # Run clustering
    for fasta in os.listdir(output_fp):
        command = 'cd-hit-est -i ' + output_fp + '/' + fasta + ' -o ' + output_fp + '/' + fasta.split('.')[0] + '_unique.fna -c ' + args.percent_id + ' -aS 1.0 -g 1 -d 0 -r 1'
        subprocess.call(command, shell=True)

    # Concatenate and get annotation files
    if args.within_lib:
        output_contigs = output_fp + '/unique_within_lib_' + args.percent_id + '.fna'
        command = 'cat ' + output_fp + '/*_unique.fna > ' + output_fp + '/unique_within_lib_' + args.percent_id + '.fna'
    elif  args.within_abx:
        output_contigs = output_fp + '/unique_within_abx_'+ args.percent_id + '.fna'
        command = 'cat ' + output_fp + '/*_unique.fna > ' + output_fp + '/unique_within_abx_'+ args.percent_id + '.fna'
    else:
        output_contigs = output_fp + '/all_unique_'+ args.percent_id + '.fna'
        command = 'cat ' + output_fp + '/*_unique.fna > ' + output_fp + '/all_unique_'+ args.percent_id + '.fna'

    subprocess.call(command, shell=True)

    # Filter annotations
    if args.annotations_fp:
        save_contig = []
        for line in open(output_contigs):
            if line.startswith('>'):
                save_contig.append(line.rstrip())
        print save_contig

        output_anno = open(output_fp + '/unique_annotations.tab', 'w')
        first = True
        for line in open(args.annotations_fp):
            if first:
                output_anno.write(line)
                first = False
            else:
                if line.rstrip() in save_contig:
                    output_anno.write(line)
            

if __name__ == "__main__":
    main()



