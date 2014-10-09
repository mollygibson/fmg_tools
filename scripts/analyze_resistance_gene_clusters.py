#!/usr/bin/env python                                                                                                                                                                                                                        
__author__ = "Molly Gibson"
__email__ = "molly.gibson@wustl.edu"

# Python imports                                                                                                                                                                                                                             
import argparse, subprocess, os, itertools, operator
# Other imports                                                                                                                                                                                                                              
from Bio import SeqIO
import parse_config, parse_mapping
import helper as h
import pandas

def main():
    parser = argparse.ArgumentParser(description="Analyze clusters for gene sharing and co-resistance")

    # arguments                                                     
    parser.add_argument('-clusters', dest="cluster_fp", help="Path to .clstr data")
    parser.add_argument('-fasta', dest="fasta_fp", help="Path to fasta of genes or proteins we care about")
    parser.add_argument('-o', dest="output_fp", help="Path to output directory")
    parser.add_argument('-m', dest="mapping_fp", help="Mapping file")
    parser.add_argument('-anno_tab', dest="anno_tab", help="Annotation tab file")
    parser.add_argument('-f', dest="override", help="Force override otuput directory", action="store_true", default=False)
    args = parser.parse_args()


    if not os.path.exists(args.output_fp):
        h.run_command("mkdir " + args.output_fp)

    # only look at the antibiotic resistance genes
    res_genes = []
    for line in open(args.fasta_fp, 'r'):
        if line.startswith(">"):
            res_genes.append(line.split()[0].strip(">"))

    # analyze the cluster
    clusters = {}
    first = True
    save_cluster = False
    for line in open(args.cluster_fp, 'r'):
        if line.startswith(">"):
            if not first and save_cluster:
                clusters[cluster_num] = cluster_data
            # Reset data
            cluster_data = []
            cluster_num = line.rstrip().strip(">")
            # Reset flags
            first = False
            save_cluster = False
        else:
            gene_member = line.split(">")[1].split("...")[0]
            if gene_member in res_genes:
                save_cluster = True
            cluster_data.append(gene_member)
    if save_cluster:
        clusters[cluster_num] = cluster_data

    sum_output = open(args.output_fp + "/summary.txt", 'w')
    lib_output = open(args.output_fp + "/libraries.txt", 'w')
    abx_output = open(args.output_fp + "/abxs.txt", 'w')
    sample_to_sample = open(args.output_fp + "/sample_to_sample.txt", 'w')
    coselection_out = open(args.output_fp + "/coselection.txt", 'w')
    shared_out = open(args.output_fp + "/shared.txt", 'w')
    genes_out = open(args.output_fp + "/gene_resistance.txt", 'w')

    sum_output.write("Cluster\tClusterSize\tAbxSize\tLibSize\tAbxs\tLibs\tAnnotations\n")
    genes_out.write("ORF\tCluster\tClusterSize\tAbxSize\tLibSize\tAbxs\tLibs\tAnnotations\n")
    coselection = {}
    shared = {}
    total_libs = {}
    total_abx = {}
    for cluster_num in clusters:
        # analyze individuals and antibiotics
        antibiotics = []
        libraries = []
        annotations = []
        for gene in clusters[cluster_num]:
            sample = gene.split(".")[0]
            antibiotic = parse_mapping.main(args.mapping_fp).abx[sample]
            library = parse_mapping.main(args.mapping_fp).id[sample]
            annotation = retrieve_annotations(gene, args.anno_tab)

            if not annotation.split("\t")[0] in annotations and not annotation == "NONE":
                annotations.append(annotation.split("\t")[0])
            
            if not antibiotic in antibiotics:
                antibiotics.append(antibiotic)

            if not library in libraries:
                libraries.append(library)
            
        sum_output.write(cluster_num + "\t" + str(len(clusters[cluster_num])) + "\t" + str(len(antibiotics)) + "\t" + str(len(libraries)) + "\t" + ";".join(antibiotics) + "\t" + ";".join(libraries) + "\t" + ";".join(annotations) + "\n")

        for gene in clusters[cluster_num]:
            genes_out.write(gene + "\t" + cluster_num + "\t" + str(len(clusters[cluster_num])) + "\t" + str(len(antibiotics)) + "\t" + str(len(libraries)) + "\t" + ";".join(antibiotics) + "\t" + ";".join(libraries) + "\t" + ";".join(annotations) + "\n")

        # Libraries
        for library in libraries:
            lib_output.write(cluster_num + "\t" + library + "\t" + ";".join(annotations) + "\n")

            if library in total_libs:
                total_libs[library] = total_libs[library] + 1
            else:
                total_libs[library] = 1
    

            for lib2 in libraries:
                if library > lib2:
                    if library in shared.keys():
                        if lib2 in shared[library].keys():
                            shared[library][lib2] = shared[library][lib2] + 1
                        else:
                            shared[library][lib2] = 1
                    else:
                        shared[library] = {}
                        shared[library][lib2] = 1

        for n in range(0,len(libraries)):
            for m in range(n+1,len(libraries)):
                sample_to_sample.write(cluster_num + "\t" + libraries[n] + "\t" + libraries[m] + "\n")
        # Abx
        for abx in antibiotics:
            abx_output.write(cluster_num + "\t" + abx + "\t" + ";".join(annotations) + "\n")

            if abx in total_abx:
                total_abx[abx] = total_abx[abx] + 1
            else:
                total_abx[abx] = 1

            for abx2 in antibiotics:
                if abx in coselection.keys():
                    if abx2 in coselection[abx].keys():
                        coselection[abx][abx2] = coselection[abx][abx2] + 1
                    else:
                        coselection[abx][abx2] = 1
                else:
                    coselection[abx] = {}
                    coselection[abx][abx2] = 1
            
    for abx in coselection:
        for abx2 in coselection[abx]:
            if abx < abx2:
                coselection_out.write(abx + "\t" + abx2 + "\t" + str(coselection[abx][abx2]) + "\t" + str(total_abx[abx]) + "\t" + str(total_abx[abx2]) +  "\n")

    for lib in shared:
        for lib2 in shared[lib]:
            shared_out.write(lib + "\t" + lib2 + "\t" + str(shared[lib][lib2]) + "\t" + str(total_libs[lib]) + "\t" + str(total_libs[lib2]) + "\n")


def retrieve_annotations(gene, anno_fp):
    gene_id = gene.split('.')[0].strip('>')
    contig = int(gene.split('.')[1].split(':')[0])
    gene_start = int(gene.split('.')[1].split(':')[1].split('-')[0])
    gene_stop = int(gene.split('.')[1].split(':')[1].split('-')[1])

    annotations = pandas.io.parsers.read_csv(anno_fp, sep="\t")

    gene_annotation = annotations[(annotations.sample==gene_id)&(annotations.contig_num==int(contig))&(annotations.start==int(gene_start))&(annotations.stop==int(gene_stop))&(annotations.database=='ResFam')]

    if not gene_annotation.empty:
        return gene_annotation[:1][["anno_id","description"]].values[0][0] + "\t" + gene_annotation[:1][["anno_id","description"]].values[0][1]
    else:
        return "NONE"


if __name__ == "__main__":
    main()
