#!/usr/bin/env python

''' Distance clustering from phylogenetic tree in newick format.

Input:
Takes newick file as input and returns leaf names of clade clusters.

Optional Input:
Number of clusters to use for Kmeans and similar algorithms.
Uniprot Fasta file to parse species and gene annotations from deflines.

Business:
Converts newick file to a distance matrix and clusters using HDBSCAN.

Output:
Distance matrix as tab separated file (TSV).
Leaf node names with cluster label as TSV file with species and gene
annotations if optional fasta file is provided.

Required Python Packages:
conda install -c conda-forge biopython
conda install -c conda-forge hdbscan
conda install -c intel scikit-learn
* scikit-learn installs with hdbscan

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: February, 2022
License :: GNU GPLv3
Copyright 2022 Roth Conrad
All rights reserved
-------------------------------------------
'''

import argparse, itertools
import pandas as pd
from Bio import Phylo
from collections import defaultdict
import hdbscan
import sklearn.cluster as cluster
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt


def newick_to_distmat(newick, outpre):
    """Reads in a newick file and returns a distance matrix as pd df"""

    # Read in and parse the newick file with BioPython Phylo package
    nwk = Phylo.read(newick, 'newick')

    # initialize dict of dicts to store distance matrix
    dstmt = defaultdict(lambda: defaultdict(list))
    # loop over all points in the tree and calculate distances
    for x, y in itertools.combinations(nwk.get_terminals(), 2):
        dst = float(nwk.distance(x, y))
        dstmt[x.name][y.name] = dst
        dstmt[y.name][x.name] = dst
    for x in nwk.get_terminals():
        dstmt[x.name][x.name] = 0
    
    # Create dataframe from dictionary
    df = pd.DataFrame(dstmt)
    # sort columns and index for symmetric matrix
    df = df.sort_index(axis=1)
    df = df.sort_index(axis=0)

    # Write matrix to file
    df.to_csv(f'{outpre}_distmat.tsv', sep='\t')

    return df


def distance_matrix_cluster(df, n_clusters, outpre):
    """Takes a distance matrix and returns clusters"""

    # HDBSCAN
    print('\n\t\t\t\tRunning HDBSCAN algorithm ...')
    hdb = hdbscan.HDBSCAN(metric='precomputed')
    hdb.fit_predict(df)

    # USE HDBSCAN clusters for n_clusters if not user provided
    if not n_clusters: n_clusters = hdb.labels_.max()

    # DBSCAN
    print('\n\t\t\t\tRunning DBSCAN algorithm ...')
    db = cluster.DBSCAN()
    db.fit_predict(df)

    '''
    # AgglomerativeClustering
    print('\n\t\t\t\tRunning Agglomerative algorithm ...')
    agl = cluster.AgglomerativeClustering(
                                        linkage='single',
                                        n_clusters=n_clusters,
                                        affinity='precomputed'
                                        )
    agl.fit_predict(df.to_numpy())

    # SpectralClustering
    print('\n\t\t\t\tRunning Spectral algorithm ...')
    spc = cluster.SpectralClustering(
                                    n_clusters=n_clusters,
                                    affinity='precomputed'
                                    )
    spc.fit_predict(df)

    # MeanShift
    print('\n\t\t\t\tRunning Mean Shift algorithm ...')
    ms = cluster.MeanShift(cluster_all=False)
    ms.fit_predict(df)

    # AffinityPropagation
    print('\n\t\t\t\tRunning Affinity Propagation algorithm ...')
    afp = cluster.AffinityPropagation(affinity='precomputed')
    afp.fit_predict(df)
    '''

    # KMeans
    print('\n\t\t\t\tRunning KMeans algorithm ...')
    km = cluster.KMeans(n_clusters=n_clusters)
    km.fit_predict(df)

    # Add the labels for each algorithm to the df.
    print('\n\t\t\t\tUpdating DataFrame and writing to file ...')
    lbld_data = {
                "Gene_Name": df.index,
                "HDBSCAN": hdb.labels_,
                "DBSCAN": db.labels_,
                #"Agglomerative": agl.labels_,
                #"Spectral": spc.labels_,
                #"MeanShift": ms.labels_,
                #"Affinity": afp.labels_,
                "K-Means": km.labels_,
                }
    df_lbld = pd.DataFrame(lbld_data)

    df_lbld.to_csv(f'{outpre}_labeled.tsv', sep='\t', index=False)

    return df


def read_fasta(fp):
    name, seq = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name: yield (name, ''.join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name: yield (name, ''.join(seq))


def add_gene_info(fasta, outpre):
    """Reads in UniProt Fasta and adds species and gene annotations"""

    # initialize dict for gene info parsed from fasta file
    gene_info = {}
    # initialize dict for fasta file output of each cluster
    byCluster = defaultdict(list)

    # parse gene info out of fasta file deflines
    # uniprot fasta has annotation and species info in the deflines
    with open(fasta, 'r') as file:
        for name, seq in read_fasta(file):
            X = name.rstrip().split('//')
            if len(X) > 1:
                uid, status, func, spec, tax = X[0][1:], X[1], X[2], X[3], X[4]
            else:
                R = 'REFERENCE SEQUENCE'
                uid, status, func, spec, tax = X[0][1:], R, R, R, R

            gene_info[uid] = [status, func, spec, tax]

    # Read the labeled dataframe file and add gene info
    data = f'{outpre}_labeled.tsv'
    outfile = f'{outpre}_annotated.tsv'
    with open(data, 'r') as inf, open(outfile, 'w') as outf:
        header = inf.readline().rstrip()
        outf.write(
            f'{header}\tStatus\tGene_function\tSpecies_name\tTaxonomy\n'
            )

        for line in inf:
            X = line.rstrip().split('\t')
            name = X[0]
            clstr = X[3] # select the kmeans row
            info = gene_info[name]
            X.extend(info)
            outf.write('\t'.join(X) + '\n')
            # add fasta sequence for cluster to byCluster dict.
            fa = f'>{info[2]}\n{info[3]}\n'
            byCluster[clstr].append(fa)

    # write fasta file for each cluster
    for cluster, seqs in byCluster.items():
        with open(f'{outpre}_{int(cluster):03}.fa', 'w') as file:
            for entry in seqs:
                file.write(entry)

    return True


def unique_lists_for_colors(outpre):
    """reads back the annotated file and writes out unique names lists"""

    file = f'{outpre}_annotated.tsv'

    clusters = {}
    genes = {}
    genus = {}
    phylum = {}
    tclass = {}

    with open(file, 'r') as f:
        header = f.readline()
        for line in f:
            X = line.rstrip().split('\t')
            clusters[X[3]] = ''
            genes[X[4]] = ''
            genus[X[5].split(' ')[0]] = ''
            try: tclass[name] = X[7].split(';')[2]
            except: tclass[name] = "n/a"
            tclass[X[7].split(';')[2]] = ''

    with open(f'{outpre}_colors_clusters.tsv', 'w') as outfile:
        for name in clusters.keys():
            outfile.write(f'{name}\n')

    with open(f'{outpre}_colors_genes.tsv', 'w') as outfile:
        for name in genes.keys():
            outfile.write(f'{name}\n')

    with open(f'{outpre}_colors_genus.tsv', 'w') as outfile:
        for name in genus.keys():
            outfile.write(f'{name}\n')

    with open(f'{outpre}_colors_phylum.tsv', 'w') as outfile:
        for name in phylum.keys():
            outfile.write(f'{name}\n')

    with open(f'{outpre}_colors_class.tsv', 'w') as outfile:
        for name in tclass.keys():
            outfile.write(f'{name}\n')

    return True


def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-i', '--input_newick_file',
        help='Please specify the input newick file!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-o', '--output_prefix',
        help='How do you want to name the output files?',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-n', '--number_of_clusters',
        help='(Optional) Number of clusters to return.',
        metavar='',
        type=int,
        required=False
        )
    parser.add_argument(
        '-f', '--uniprot_fasta_file',
        help='(Optional) Fasta file from Uniprot.',
        metavar='',
        type=str,
        required=False
        )
    args=vars(parser.parse_args())

    # Do what you came here to do:
    print('\n\nRunning Script ...')

    # define input params
    newick = args['input_newick_file']
    outpre = args['output_prefix']
    n_clusters = args['number_of_clusters']
    fasta = args['uniprot_fasta_file']

    # generate distance matrix from newick file
    print('\n\t\tParsing Newick File and Converting to Distance Matrix ...')
    df = newick_to_distmat(newick, outpre)

    # cluster the distance matrix
    print('\n\t\tClustering the Distance Matrix ...')
    df_clstrd = distance_matrix_cluster(df, n_clusters, outpre)

    # Add species and gene annotations to dataframe from Uniprot fasta file
    if fasta:
        print('\n\t\tAdding species and gene annotations to DataFrame ...')
        _ = add_gene_info(fasta, outpre)

    # write out unique name sets for tree color annotation
    _ = unique_lists_for_colors(outpre)

    print(f'\n\nComplete success space cadet!! Finished without errors.\n\n')


if __name__ == "__main__":
    main()
