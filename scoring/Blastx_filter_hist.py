#!/usr/bin/env python

'''Plots histograms from best hit filtered tabular blastx file

Plots histograms for:

            'pid': 'Sequence Identity Histogram (%)',
            'alen': 'Alignment Length Histogram (AA)',
            'pml': 'Sequence Alignment Ratio (%)',
            'qlen': 'Query Length Histogram (AA)'

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

import argparse, random
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt


def parse_blast(file):
    """ parse blast file for pidents returns list of floats """

    data = {'pid': [], 'alen': [], 'pml': [], 'qlen': [], 'bitscore': []}

    with open(file, 'r') as f:
        for l in f:
            X = l.rstrip().split('\t')
            pident = float(X[2])
            alen = int(X[3])
            bitscore = float(X[11])
            qlen = int(X[12]) / 3 # blastx amino acid correction
            pml = alen / qlen

            data['pid'].append(pident)
            data['alen'].append(alen)
            data['pml'].append(pml)
            data['qlen'].append(qlen)
            data['bitscore'].append(bitscore)

    return data


def plot_hist(data, outpre, key):

    # Define the titles
    plot_titles = {
                    'pid': 'Sequence Identity Histogram (%)',
                    'alen': 'Alignment Length Histogram (AA)',
                    'pml': 'Sequence Alignment Ratio (%)',
                    'qlen': 'Query Length Histogram (AA)',
                    'bitscore': 'Bitscore Histogram'
                    }

    print(f'\t\tPlotting {plot_titles[key]}...')

    # Set the colors
    bar_color = '#2171b5'
    gridM = '#bdbdbd'
    gridm = '#d9d9d9'
    alpha = 0.6

    # Build the plot
    fig, ax = plt.subplots(figsize=(10, 7))

    # Plot titles
    ax.set_title(
        plot_titles[key],
        fontsize=20, y=1.02
        )

    # Plot labels
    ax.set_ylabel('Count', fontsize=14)
    ax.set_xlabel('Value', fontsize=14)

    # Set plot/grid style
    ax.minorticks_on()
    ax.tick_params(
        which='minor', axis='both', left=False, bottom=False
        )
    ax.tick_params(
                which='major', axis='both',
                left=True, bottom=True,
                size=6, width=2, tickdir='inout',
                labelsize=12, zorder=10
                )
    ax.yaxis.grid(
        which="minor", color=gridm, linestyle='--',
        linewidth=1, alpha=0.6, zorder=1
        )
    ax.yaxis.grid(
        which="major", color=gridM, linestyle='--',
        linewidth=1.5, alpha=0.4, zorder=1
        )
    ax.set_axisbelow(True)
    for spine in ax.spines.values(): spine.set_linewidth(2)

    # Plot the data
    ax.hist(
        data[key],
        bins=30,
        #orientation='horizontal',
        rwidth=0.9,
        color=bar_color,
        alpha=alpha,
        )

    # Set plot axis ranges
    #ax.set_xlim(left=0, right=int((max(d['xs'])+min(d['xs']))))

    # adjust layout, save, and close
    #plt.gca().invert_yaxis()
    fig.set_tight_layout(True)
    plt.savefig(f'{outpre}_{key}_histogram.pdf')
    plt.close()


def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-i', '--in_file',
        help='Please specify the tabular magic blast file!',
        metavar='',
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
    print('\n\nRunning Script...\n')


    data = parse_blast(args['in_file'])

    for key in ['pid', 'alen', 'pml', 'qlen', 'bitscore']:
        _ = plot_hist(data, args['output_prefix'], key)

    print('\n\nComplete success space cadet! Hold on to your boots.\n\n')

if __name__ == "__main__":
    main()
