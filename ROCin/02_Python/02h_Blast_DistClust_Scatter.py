#!/usr/bin/env python

'''Plot Scatter plot of Blast result colored by distance clustering.

This tool takes the following input parameters:

    * Sequences.blast - tabular Blast File
    * newick distance matrix cluster result

This script returns the following files:

    * 3 scatter Plot .pdf format

This script requires the following packages:

    * argparse
    * matplotlib
    * seaborn

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: March 2022
License :: GNU GPLv3
Copyright 2022 Roth Conrad
All rights reserved
-------------------------------------------
'''


import argparse
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd


def parse_blast(file):
    """ parse blast file for pidents returns dict {query: [pident, pml]} """

    data = {}

    with open(file, 'r') as f:
        header = f.readline()
        for l in f:
            X = l.rstrip().split('\t')
            name = X[0]
            pident = float(X[1])
            alen = int(X[2])
            qlen = int(X[3])
            pml = alen / qlen

            data[name] = [pident, pml]

    return data


def parse_distclust(file, blast):
    """ parse distance clustered tsv file returns longform df """

    data = {'Name': [], 'Cluster': [], 'Label': [], 'Pident': [], 'PML': []}

    with open(file, 'r') as f:
        header = f.readline()
        for l in f:
            X = l.rstrip().split('\t')
            name = X[0]
            cluster = X[1]
            label = X[2]

            if name in blast:
                blst = blast[name]       

                data['Name'].append(name)
                data['Cluster'].append(cluster)
                data['Label'].append(label)
                data['Pident'].append(blst[0])
                data['PML'].append(blst[1])

    df = pd.DataFrame(data).set_index('Name')
    df['Cluster'] = df['Cluster'].astype(str)
    print('\n\n', df)

    return(df)


def plot_scatter(df, out):

    sns.set_style("ticks")
    sns.set_context("talk")
    sns.color_palette()

    clust_order = [
                'LnuF', 'LnuB', 'LnuG', '0', '1', '2', '3', '4', '5', '6', '7',
                '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18',
                '19'
                ]

    g = sns.relplot(
        x="Pident", y="PML", kind="scatter", data=df,
        #style='Label', legend=False,
        col='Cluster', col_order=clust_order, col_wrap=6,
        #style_order=['NCBI', 'TrEMBL'],
        markers=["s", "o"],
        s=200,
        alpha=0.3,
        #hue='Cluster', 
        #hue_order=clust_order,
        )
    #g.ax.legend(ncol=2)

    # adjust layout, save, and close
    #sns.despine(offset=10, trim=True)
    plt.savefig(out)
    plt.close()


def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-b', '--blast_input_file',
        help='Please specify the input file name!',
        metavar=':',
        type=str,
        required=True
        )
    parser.add_argument(
        '-d', '--distclust_input_file',
        help='Please specify the input file name!',
        metavar=':',
        type=str,
        required=True
        )
    parser.add_argument(
        '-o', '--output_prefix',
        help='Please specify the prefix to use for output files!',
        metavar=':',
        type=str,
        required=True
        )
    args=vars(parser.parse_args())

    # Do what you came here to do:
    print('\nRunning Script...\n')
    
    blast = args['blast_input_file']
    distclust = args['distclust_input_file']
    out = args['output_prefix']

    print('\n\nParsing Blast input file ...')
    data = parse_blast(blast)
    print('\n\nParsing DistClust input file ...')
    df = parse_distclust(distclust, data)
    print('\n\nPlotting ... and scheming ... maybe ...')
    _ = plot_scatter(df, out)

    print('\n\nComplete success space cadet! Hold on to your boots.\n\n')


if __name__ == "__main__":
    main()
