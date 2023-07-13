#!/usr/bin/env python

'''Score the rocker model against static filters

Takes inputs:

1) original mock metagenome fasta
2) best hit filtered hmmsearch result (besthit_filter_hmm.py)
3) best hit filtered blastx alignment (besthit_filter_blast.py)
4) best hit filtered passing rocker alignments
5) best hit filtered failing rocker alignments

* expects inputs to be best hit filtered.
* newest rocker should have built in best hit filter.

Outputs:

1) tsv data table
2) pdf bar plot

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: April 2023
License :: GNU GPLv3
Copyright 2023 Roth Conrad
All rights reserved
-------------------------------------------
'''

import argparse
import pandas as pd
import matplotlib.pyplot as plt

def read_fasta(fp):
    ''' parses a fasta file into name, seq '''
    name, seq = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name: yield (name, ''.join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name: yield (name, ''.join(seq))


def parse_mm(mm):
    ''' Read fasta and counts positive labeled reads. returns count '''

    data = {}
    pos_count = 0

    with open(mm, 'r') as file:
        for name, seq in read_fasta(file):
            name = name[1:]
            label = name.split(';')[-1]
            data[label] = ''
            if label == 'Positive':
                pos_count += 1

    print('Labels in this dataset:', list(data.keys()))

    return pos_count


def score_results(TP, FP, TN, FN):

    P = TP + FP
    N = TN + FN
    Total = P + N
    FNR = FN/(TP + FN)
    FPR = FP/P
    Sensitivity = 1 - FNR
    Specificity = 1 - FPR
    Accuracy = (TP + TN) / (P + N)
    Precision = TP / (TP + FP)
    Recall = TP / (TP + FN)
    F1 = 2 * (Precision * Recall) / (Precision + Recall)

    results = [
                Total, P, N, TP, FP, TN, FN, FNR, FPR,
                Sensitivity, Specificity, Accuracy,
                Precision, Recall, F1
                ]

    return results


def score_blastx(bx, all_pos_reads):

    customA = {'TP': 0, 'FP': 0, 'TN': 0, 'FN': 0}
    customB = {'TP': 0, 'FP': 0, 'TN': 0, 'FN': 0}
    e30 = {'TP': 0, 'FP': 0, 'TN': 0, 'FN': 0}
    e20 = {'TP': 0, 'FP': 0, 'TN': 0, 'FN': 0}
    e10 = {'TP': 0, 'FP': 0, 'TN': 0, 'FN': 0}
    e03 = {'TP': 0, 'FP': 0, 'TN': 0, 'FN': 0}

    pos_count = 0

    with open(bx, 'r') as file:
        for line in file:
            X = line.rstrip().split('\t')
            query = X[0] # read identifier
            label = query.split(';')[-1]
            evalue = float(X[10]) # evalue
            pid = float(X[2]) # percent sequence identity
            aLen = int(X[3]) # read alignment length
            qLen = int(X[12]) / 3 # full length of read divided by 3 amino acids
            pMatch = aLen / qLen # percent match length of read length

            # count positive reads in blastX alignment
            if label == 'Positive':
                pos_count += 1

            # count customA filter
            if pMatch >= 0.7 and pid >= 95:
                if label == 'Positive':
                    customA['TP'] += 1
                else: customA['FP'] += 1
            else:
                if label == 'Positive':
                    customA['FN'] += 1
                else: customA['TN'] += 1

            # count customB filter
            if pMatch >= 0.5 and pid >= 90:
                if label == 'Positive':
                    customB['TP'] += 1
                else: customB['FP'] += 1
            else:
                if label == 'Positive':
                    customB['FN'] += 1
                else: customB['TN'] += 1

            # count e30 filter
            if evalue <= 1e-30:
                if label == 'Positive':
                    e30['TP'] += 1
                else: e30['FP'] += 1
            else:
                if label == 'Positive':
                    e30['FN'] += 1
                else: e30['TN'] += 1

            # count e20 filter
            if evalue <= 1e-20:
                if label == 'Positive':
                    e20['TP'] += 1
                else: e20['FP'] += 1
            else:
                if label == 'Positive':
                    e20['FN'] += 1
                else: e20['TN'] += 1

            # count e10 filter
            if evalue <= 1e-10:
                if label == 'Positive':
                    e10['TP'] += 1
                else: e10['FP'] += 1
            else:
                if label == 'Positive':
                    e10['FN'] += 1
                else: e10['TN'] += 1

            # count e03 filter
            if evalue <= 1e-03:
                if label == 'Positive':
                    e03['TP'] += 1
                else: e03['FP'] += 1
            else:
                if label == 'Positive':
                    e03['FN'] += 1
                else: e03['TN'] += 1

    # group the data for ease
    data = [customA, customB, e30, e20, e10, e03]

    # check all positive reads are counted
    print('\n\nPositive read count from metagenome:', all_pos_reads)
    print('\nPositive read count from BlastX alignment:', pos_count)
    pos_diff = all_pos_reads - pos_count
    print('Positive reads not aligned by BlastX:', pos_diff)

    return data


def score_rocker(rp, rf, all_pos_reads):

    roc = {'TP': 0, 'FP': 0, 'TN': 0, 'FN': 0}

    passing_reads = {}
    pos_count = 0

    # rocker passing TP, FP
    with open(rp, 'r') as file:
        for line in file:
            X = line.rstrip().split('\t')
            query = X[0] # read identifier
            label = query.split(';')[-1]
            evalue = float(X[10]) # evalue
            passing_reads[query] = ''

            if label == 'Positive':
                roc['TP'] += 1
                pos_count += 1
            else: roc['FP'] += 1

    # rocker failing FN, TN
    with open(rf, 'r') as file:
        for line in file:
            X = line.rstrip().split('\t')
            query = X[0] # read identifier
            label = query.split(';')[-1]
            evalue = float(X[10]) # evalue

            if query in passing_reads: continue

            if label == 'Positive':
                roc['FN'] += 1
                pos_count += 1
            else: roc['TN'] += 1

    # check all positive reads are counted
    print('\nPositive read count from ROCkOut:', pos_count)
    pos_diff = all_pos_reads - pos_count
    print('Positive reads not reported by ROCkOut:', pos_diff)

    return roc


def score_hmm(hm, all_pos_reads):

    hDf = {'TP': 0, 'FP': 0, 'TN': 0, 'FN': 0}
    h30 = {'TP': 0, 'FP': 0, 'TN': 0, 'FN': 0}
    h20 = {'TP': 0, 'FP': 0, 'TN': 0, 'FN': 0}
    h10 = {'TP': 0, 'FP': 0, 'TN': 0, 'FN': 0}
    h03 = {'TP': 0, 'FP': 0, 'TN': 0, 'FN': 0}

    pos_count = 0

    with open(hm, 'r') as file:
        for line in file:
            X = line.rstrip().split('\t')
            query = X[0] # read identifier
            label = query.split(';')[-1]
            evalue = float(X[1]) # evalue

            # count hmm default
            if evalue <= 0.01:
                if label == 'Positive':
                    hDf['TP'] += 1
                else: hDf['FP'] += 1
            else:
                if label == 'Positive':
                    hDf['FN'] += 1
                else: hDf['TN'] += 1

            # count h30 filter
            if evalue <= 1e-30:
                if label == 'Positive':
                    h30['TP'] += 1
                else: h30['FP'] += 1
            else:
                if label == 'Positive':
                    h30['FN'] += 1
                else: h30['TN'] += 1

            # count h20 filter
            if evalue <= 1e-20:
                if label == 'Positive':
                    h20['TP'] += 1
                else: h20['FP'] += 1
            else:
                if label == 'Positive':
                    h20['FN'] += 1
                else: h20['TN'] += 1

            # count h10 filter
            if evalue <= 1e-10:
                if label == 'Positive':
                    h10['TP'] += 1
                else: h10['FP'] += 1
            else:
                if label == 'Positive':
                    h10['FN'] += 1
                else: h10['TN'] += 1

            # count h03 filter
            if evalue <= 1e-03:
                if label == 'Positive':
                    h03['TP'] += 1
                else: h03['FP'] += 1
            else:
                if label == 'Positive':
                    h03['FN'] += 1
                else: h03['TN'] += 1

            # Count positive reads mapped by hmmsearch
            if label == 'Positive':
                pos_count += 1

    # group the data for ease
    data = [h30, h20, h10, h03, hDf]

    # check all positive reads are counted
    print('\nPositive read count from hmmsearch:', pos_count)
    pos_diff = all_pos_reads - pos_count
    print('Positive reads not aligned by hmmsearch:', pos_diff)

    return data


def build_bar_plots(df, out):

    fs = 12 # set font size
    metrics = ['F1', 'FPR', 'FNR']
    labels = df.columns.to_list()

    fig, axes = plt.subplots(len(metrics), 1, figsize=(4.25, 10), sharex=True)

    for i, met in enumerate(metrics):
        ax = axes[i]
        data = df.T[met].to_list()
        _ = ax.bar(labels, data, color='#636363', width=0.5, alpha=0.75)
        ax.set_ylabel(met, fontsize=fs)
        ax.set_ylim([0, 1])
        ax.vlines(x=0.5, ymin=0, ymax=1, color='#000000', ls='--', lw=1.5,)
        ax.vlines(x=5.5, ymin=0, ymax=1, color='#000000', ls='--', lw=1.5,)
        ax.text(3, 0.92, 'BLASTx', fontweight='black', ha='center')
        ax.hlines(y=0.91, xmin=1.9, xmax=4.1, color='#000000', ls='-', lw=1.5)
        ax.text(8.2, 0.92, 'HMM', fontweight='black', ha='center')
        ax.hlines(y=0.91, xmin=7.45, xmax=8.95, color='#000000', ls='-', lw=1.5)
        ax.yaxis.grid(
            which="major", color='#d9d9d9', linestyle='--',
            linewidth=1, zorder=1
            )
        ax.set_axisbelow(True)

    plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
         rotation_mode="anchor")

    fig.set_tight_layout(True)
    plt.savefig(f'{out}_bar_plot.pdf')
    plt.close()

    return True


def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-mm', '--mock_metagenome_fasta',
        help='Original mock metagenome fasta file.',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-bx', '--blastx_alignments_file',
        help='Best hit filtered blastx alignment file.',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-hm', '--filtered_hmmer_file',
        help='best hit filtered hmmsearch results.',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-rp', '--rocker_passing_alignments',
        help='BlastX alignments passing the rocker filter.',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-rf', '--rocker_failing_alignments',
        help='BlastX alignments failing the rocker filter.',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-o', '--output_prefix',
        help='Prefix to use for output files.',
        metavar='',
        type=str,
        required=True
        )
    args=vars(parser.parse_args())

    # define input params
    mm = args['mock_metagenome_fasta'] # original mock metagenome fasta
    bx = args['blastx_alignments_file'] # best hit filtered blastx file
    hm = args['filtered_hmmer_file'] # best hit filtered hmmsearch file
    rp = args['rocker_passing_alignments'] # rocker passing alignments
    rf = args['rocker_failing_alignments'] # rocker failing alignments
    out = args['output_prefix']  # output file prefix

    # Do what you came here to do:
    print('\n\nRunning Script...\n\n')

    # Get all reads labeled as positive
    all_pos_reads = parse_mm(mm)
    # Count TP, FP, TN, FN for blastx
    bx_data = score_blastx(bx, all_pos_reads)
    bx_results = []
    for d in bx_data:
        r = score_results(d['TP'], d['FP'], d['TN'], d['FN'])
        bx_results.append(r)

    # Count TP, FP, TN, FN for rocker
    roc_data = score_rocker(rp, rf, all_pos_reads)
    roc_results = score_results(
                            roc_data['TP'], roc_data['FP'],
                            roc_data['TN'], roc_data['FN']
                            )

    # Count TP, FP, TN, FN for hmm
    hm_data = score_hmm(hm, all_pos_reads)
    hm_results = []
    for d in hm_data:
        r = score_results(d['TP'], d['FP'], d['TN'], d['FN'])
        hm_results.append(r)

    # output table and pdf
    scores = {
                'ROCker': roc_results, 'Custom-A': bx_results[0],
                'Custom-B': bx_results[1], 'evalue 1e-30': bx_results[2],
                'evalue 1e-20': bx_results[3], 'evalue 1e-10': bx_results[4],
                'evalue 1e-3': bx_results[5], 'HMM 1e-30': hm_results[0],
                'HMM 1e-20': hm_results[1], 'HMM 1e-10': hm_results[2],
                'HMM 1e-3': hm_results[3], 'HMM default': hm_results[4],
                }
    rows = [
            'Total', 'P', 'N', 'TP', 'FP', 'TN', 'FN', 'FNR', 'FPR',
            'Sensitivity', 'Specificity', 'Accuracy',
            'Precision', 'Recall', 'F1'
            ]

    df = pd.DataFrame(scores, index=rows).round(2)
    df.to_csv(f'{out}_score_table.tsv', sep='\t')

    print('\n\n', df)

    _ = build_bar_plots(df, out)

    print('\n\nComplete success space cadet! Hold on to your boots.\n\n')

if __name__ == "__main__":
    main()
