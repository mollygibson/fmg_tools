#!/usr/bin/env python                                                                                                                                                          

__author__ = "Molly Gibson"
__email__ = "molly.gibson@wustl.edu"


# Python imports
import argparse, subprocess, os, itertools, operator
import parse_config

def main():
    parser = argparse.ArgumentParser(description="This script analyzes output from \
              antibiotic functional metagenomic selections assembled using PARFuMS.")

    # arguments
    parser.add_argument('-contigs', dest="contig_fp", help="Path to contig file from PARFuMS output")
    parser.add_argument('-proteins', dest="protein_fp", help="Path to protein file from PARFums output")
    parser.add_argument('-o', dest="output_fp", help="Path to output directory for annotation")
    args = parser.parse_args()
    
    # check to see that either contigs or proteins are provided
    if not args.contig_fp and not args.protein_fp:
        parser.exit(status=0, message="You must provide either contigs or proteins to annotate\nCheck usage with 'annotate_functional_selections.py -h'.\n\n")

    # set the output path and make the directory
    prefix = os.path.splitext(args.contig_fp)[0]
    if not args.output_fp:
        output_fp = os.path.splitext(args.contig_fp)[0] + "_annotation"
    else:
        output_fp = args.output_fp
    if not os.path.isdir(output_fp):
        subprocess.call('mkdir ' + output_fp, shell=True)

    # call orfs (if contigs) and annotate proteins
    if not args.protein_fp:
        protein_fp = call_orfs(args.contig_fp, output_fp, prefix)
    else:
        protein_fp = args.protein_fp
    annotation_fp = annotate_proteins(protein_fp, output_fp)
    
    # merge hmm files
    merge_hmm_files(output_fp, prefix)


def call_orfs(contigs, output_fp, prefix):    

    if not os.path.isfile(output_fp + '/sequence.gff'):
        command = 'gmhmmp -a -d -f G -m ' + parse_config.main().mgm_model_fp + ' -o ' + output_fp + '/sequence.gff ' + contigs
        subprocess.call(command, shell=True)

    if not os.path.isfile(output_fp + '/nucleotides.fasta'):
        command = 'nt_from_gff.pl < ' + output_fp + '/sequence.gff > ' + output_fp + '/nucleotides.fasta'
        subprocess.call(command, shell=True)

    if not os.path.isfile(output_fp + '/proteins.fasta'):
        command = 'aa_from_gff.pl < ' + output_fp + '/sequence.gff > ' + output_fp + '/proteins.fasta'
        subprocess.call(command, shell=True)

    nucleotideFile = open(output_fp + '/nucleotides.fasta', 'r')
    newNucFile = open(output_fp + '/nucleotides_updated.fasta', 'w')
    proteinFile = open(output_fp + '/proteins.fasta', 'r')
    sequenceFile = open(output_fp + '/sequence.gff', 'r')
    newProteinFilePath = output_fp + '/proteins_updated.fasta'
    newProteinFile = open(output_fp + '/proteins_updated.fasta', 'w')
    contigNames = open(output_fp + '/contig_names.txt', 'w')
    
    for line in nucleotideFile:
        if line.startswith('>'):
            # Create new header                                                                                                                                                                                                                   
            header = line.strip('>').split('_')
            searchS = header[0] + '_' + header[1] + " " +  header[2]
            # Find the line in sequenceFile that contains searchS                                                                                                                                                                                 
            for headLine in sequenceFile:
                if not headLine.startswith('#'):
                    if searchS in headLine:
                        sequenceHeader = headLine
                        break
            headerPieces = sequenceHeader.split('\t')
            newHeader = "> ID:" + headerPieces[0] + " abx:" + headerPieces[3] + " start:" + headerPieces[4] + " stop:" + headerPieces[6]
            newNucFile.write(newHeader + "\n")
            contigNames.write(newHeader + "\n")
        else:
            newNucFile.write(line)

    newNucFile.close()
    contigNames.seek(0)
    sequenceFile.seek(0)
    
    for line in proteinFile:
        if line.startswith('>'):
            # Create new header
            header = line.strip('>').split('_')
            searchS = header[0] + '_' + header[1] + " " +  header[2]
            # Find the line in sequenceFile that contains searchS
            for headLine in sequenceFile:
                if not headLine.startswith('#'):
                    if searchS in headLine:
                        sequenceHeader = headLine
                        break
            headerPieces = sequenceHeader.split('\t')
            newHeader = ">" + headerPieces[0] + "_" + headerPieces[3] + "_" + headerPieces[4] + "_" + headerPieces[6]
            newProteinFile.write(newHeader + "\n")
            contigNames.write(newHeader + "\n")
        else:
            newProteinFile.write(line)
        
    newProteinFile.close()
    contigNames.close()

    return(newProteinFilePath)

# annotate proteins
def annotate_proteins(proteins,output_fp):
    if not os.path.isfile(output_fp + '/Pfam-targets.txt'):
        command = 'hmmscan -o /dev/null --cut_ga --tblout ' + output_fp + '/Pfam-targets.txt ' + parse_config.main().pfam_database_fp + ' ' +  proteins
        subprocess.call(command, shell=True)
    if not os.path.isfile(output_fp + '/Tigrfam-targets.txt'):
        command = 'hmmscan -o /dev/null --cut_ga --tblout ' + output_fp + '/Tigrfam-targets.txt ' + parse_config.main().tigrfam_database_fp + ' ' + proteins
        subprocess.call(command, shell=True)
    if not os.path.isfile(output_fp + '/ResFams-targets.txt'):
        command = 'hmmscan -o /dev/null --cut_ga --tblout ' + output_fp + '/ResFams-targets.txt ' + parse_config.main().resfam_database_fp + ' ' + proteins
        subprocess.call(command, shell=True)

    # merge hmm annotations 

def merge_hmm_files(output_fp, prefix):
    pfam = open(output_fp + '/Pfam-targets.txt', 'r')
    tigrfam = open(output_fp + '/Tigrfam-targets.txt', 'r')
    resfams = open(output_fp + '/ResFams-targets.txt', 'r')
    output = open(output_fp + '/' + prefix + '.HMMAnnotation.txt', 'w')
    contigs = open(output_fp + '/contig_names.txt', 'r') 

    for line in contigs:
        pfam.seek(0)
        tigrfam.seek(0)
        resfams.seek(0)

        hits = []
        contig = line.strip('>').rstrip()
            # Process Pfam File                                                                                                                                                                                                            
        for item in pfam:
            if contig in item:
                pieces = item.rstrip().split()
                description = pieces[18]
                if len(pieces) > 19:
                    description = " ".join(pieces[18:len(pieces)])
                array = ["Pfam", pieces[1], pieces[2], pieces[3], float(pieces[4]), float(pieces[5]), float(pieces[6]), float(pieces[7]), float(pieces[8]), float(pieces[9]), description, pieces[0]]
                hits.append(array)

            # Process TIGRFam File                                                                                                                                                                                                         
        for item in tigrfam:
            if contig in item:
                pieces = item.rstrip().split()
                array = ["TIGRFam", pieces[1], pieces[2], pieces[3], float(pieces[4]), float(pieces[5]), float(pieces[6]), float(pieces[7]), float(pieces[8]), float(pieces[9]), " ".join(pieces[18:len(pieces)]), pieces[0] ]
                hits.append(array)

        # Process ResFams Files                                                                                                                                                                                                            
        for item in resfams:
            if contig in item:
                pieces = item.rstrip().split()
                if pieces[1].startswith('P'):
                    array = ["ResFam", pieces[1], pieces[2], pieces[3], float(pieces[4]), float(pieces[5]), float(pieces[6]), float(pieces[7]), float(pieces[8]), float(pieces[9]), " ".join(pieces[18:len(pieces)])]
                elif pieces[1].startswith('T'):
                    array = ["ResFam", pieces[1], pieces[2], pieces[3], float(pieces[4]), float(pieces[5]), float(pieces[6]), float(pieces[7]), float(pieces[8]), float(pieces[9]), " ".join(pieces[18:len(pieces)])]
                else:
                    array = ["ResFam", pieces[1], pieces[2], pieces[3], float(pieces[4]), float(pieces[5]), float(pieces[6]), float(pieces[7]), float(pieces[8]), float(pieces[9]), pieces[0]]
                hits.append(array)

        output.write(">" +  contig.rstrip() + "\n")

        final_list = []
        for item in hits:
            pieces = item[2].split('_')
            start = pieces[len(pieces)-3]
            stop = pieces[len(pieces)-2]
            final_list.append([start, stop, item[0], item[1], item[10], item[4], item[5]])

        pieces = contig.rstrip().split('_')
        print pieces
        output.write(pieces[len(pieces)-3] + "\t" + pieces[len(pieces)-2] + "\t" + pieces[len(pieces)-1] + "\n")

        sorted_hits = sorted(final_list, reverse=True,  key=lambda hit: hit[6])

        for key,group in itertools.groupby(sorted_hits,operator.itemgetter(6)):
            sorted_group = sorted(list(group), key=lambda hit: hit[5])

            for item in sorted_group:
                if item:
                    output.write("\t" + str(item[2]) + "\t" + str(item[3]) + "\t" + str(item[4]) + "\t" + str(item[5]) + "\t" + str(item[6]) + "\n")


if __name__ == "__main__":
    main()
