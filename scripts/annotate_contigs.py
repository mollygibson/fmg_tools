#!/usr/bin/env python

__author__ = "Molly Gibson"
__email__ = "molly.gibson@wustl.edu"

import subprocess
import os
import parse_config

# Call orfs using MetaGeneMark
def call_orfs(contigs, prefix):    

    if not os.path.isdir(prefix + '_annotation_temp_files'):
        subprocess.call('mkdir ' + prefix + '_annotation_temp_files', shell=True)

    if not os.path.isfile(prefix + '_annotation_temp_files/sequence.gff'):
        command = 'gmhmmp -a -d -f G -m ' + parse_config.main().mgm_model_fp + ' -o ' + prefix + '_annotation_temp_files/sequence.gff ' + contigs
        subprocess.call(command, shell=True)

    if not os.path.isfile(prefix + '_annotation_temp_files/nucleotides.fasta'):
        command = 'nt_from_gff.pl < ' + prefix + '_annotation_temp_files/sequence.gff > ' + prefix + '_annotation_temp_files/nucleotides.fasta'
        subprocess.call(command, shell=True)

    if not os.path.isfile(prefix + '_annotation_temp_files/proteins.fasta'):
        command = 'aa_from_gff.pl < ' + prefix + '_annotation_temp_files/sequence.gff > ' + prefix + '_annotation_temp_files/proteins.fasta'
        subprocess.call(command, shell=True)

    nucleotideFile = open(prefix + '_annotation_temp_files/nucleotides.fasta', 'r')
    newNucFile = open(prefix + '_annotation_temp_files/nucleotides_updated.fasta', 'w')
    proteinFile = open(prefix + '_annotation_temp_files/proteins.fasta', 'r')
    sequenceFile = open(prefix + '_annotation_temp_files/sequence.gff', 'r')
    newProteinFilePath = prefix + '_annotation_temp_files/proteins_updated.fasta'
    newProteinFile = open(prefix + '_annotation_temp_files/proteins_updated.fasta', 'w')
    contigNames = open(prefix + '_annotation_temp_files/contig_names.txt', 'w')
    
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
            newHeader = ">" + headerPieces[0] + "_" + headerPieces[3] + "_" + headerPieces[4] + "_" + headerPieces[6]
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
def annotate_proteins(proteins, prefix):
    if not os.path.isfile(prefix + '_annotation_temp_files/Pfam-targets.txt'):
        command = 'hmmscan -o /dev/null --cut_ga --tblout ' + prefix + '_annotation_temp_files/Pfam-targets.txt ' + parse_config.main().pfam_database_fp + ' ' +  proteins
        subprocess.call(command, shell=True)
    if not os.path.isfile(prefix + '_annotation_temp_files/Tigrfam-targets.txt'):
        command = 'hmmscan -o /dev/null --cut_ga --tblout ' + prefix + '_annotation_temp_files/Tigrfam-targets.txt ' + parse_config.main().tigrfam_database_fp + ' ' + proteins
        subprocess.call(command, shell=True)
    if not os.path.isfile(prefix + '_annotation_temp_files/ResFams-targets.txt'):
        command = 'hmmscan -o /dev/null --cut_ga --tblout ' + prefix + '_annotation_temp_files/ResFams-targets.txt ' + parse_config.main().resfam_database_fp + ' ' + proteins
        subprocess.call(command, shell=True)

    # merge hmm annotations 


