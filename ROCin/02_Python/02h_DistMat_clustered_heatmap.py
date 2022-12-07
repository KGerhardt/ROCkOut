#!/usr/bin/env python

'''Plots clustered heatmap of percent identity matrix (P.I.M.).

This tool takes the following input parameters:

    * Sequences.dst - percent identity or distance matrix in tsv format.
      - Should not have column headers
      - should be a square matrix - all verse all
      - columns should be separated with tabs (tsv format)
      - each line is row (or rows are separated by new lines)
      - 1st column should be sequence name and is used as row names
      - the row names will be used as the column names
      - the column order should be the same as the row order

This script returns the following file:

    * clustered heatmap in .pdf format

This script requires the following packages:

    * argparse
    * pandas
    * matplotlib
    * seaborn

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: April 2022
License :: GNU GPLv3
Copyright 2022 Roth Conrad
All rights reserved
-------------------------------------------
'''


import argparse
import pandas as pd
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
import seaborn as sns


def plot_clustred_heatmap(df, outfile):

    # build the plot
    g = sns.clustermap(
                    df, figsize=(16,9),
                    xticklabels=True,
                    yticklabels=True
                    )

    # adjust layout, save, and close
    #plt.gca().invert_yaxis()
    #fig.set_tight_layout(True)
    g.savefig(outfile)
    plt.close()


def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-i', '--input_file',
        help='Please specify the input file name!',
        metavar=':',
        type=str,
        required=True
        )
    parser.add_argument(
        '-o', '--output_file',
        help='Please specify the output file name!',
        metavar=':',
        type=str,
        required=True
        )
    args=vars(parser.parse_args())

    # Do what you came here to do:
    print('\n\nRunning Script...')
    
    # define parameters
    infile = args['input_file']
    outfile = args['output_file']

    # read in the tsv file with pandas.
    df = pd.read_csv(infile, sep='\t', header=0, index_col=0)
    #df.columns = df.index

    print('\n\n', df)

    # create the plot
    _ = plot_clustred_heatmap(df, outfile)

    print('\n\nComplete success space cadet! Hold on to your boots.\n\n')


if __name__ == "__main__":
    main()
