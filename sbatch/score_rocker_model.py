#!/usr/bin/env python

'''Score the rocker model against static filters

writes out a tsv table and a pdf bar plot of results.

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

import argparse, random
import pandas as pd
import matplotlib.pyplot as plt


def best_hits(query, bitscore, d, line):
    """ Filters the besthit based on bitscore """

    if query in d:
        old_bitscore = float(d[query][0].split('\t')[11])

        if bitscore > old_bitscore:
            d[query] = [line]

        elif bitscore == old_bitscore:
            d[query].append(line)

    else:
        d[query] = [line]

    return d


def best_hits_filter(infile):

    d = {} # initialize dictionary for bitscore besthits
    total = 0

    with open(infile, 'r') as file:

        for line in file:
            total += 1
            X = line.rstrip().split('\t')
            query = X[0] # read identifier
            bitscore = float(X[11]) # bitscore

            d = best_hits(query, bitscore, d, line)

    # store best hits
    bhs = []
    for k,v in d.items():
        bhs.append(random.choice(v))

    print(f'{infile}: {total}, {len(bhs)}')

    return bhs


def dedup_pos_neg(pos, neg):

    pos_reads = {}

    P, N, D = 0, 0, 0

    for line in pos:
        label = ';'.join(line.split(';')[:-1])
        pos_reads[label] = ''
        P += 1

    for line in neg:
        label = ';'.join(line.split(';')[:-1])
        if label in pos_reads:
            D += 1
        else:
            N += 1

    return P, N, D


def score_results(TP, FN, FP, TN):

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


def get_read_lengths(prl, nrl):

    read_lens = {} # readname: readlength

    with open(prl, 'r') as file:
        for name, seq in read_fasta(file):
            read_lens[name[1:]] = len(seq)

    with open(nrl, 'r') as file:
        for name, seq in read_fasta(file):
            read_lens[name[1:]] = len(seq)

    return read_lens


def static_filter_blast(result, read_lens):

    strpass, strfail = 0, 0
    e30pass, e30fail = 0, 0
    e20pass, e20fail = 0, 0
    e10pass, e10fail = 0, 0
    e03pass, e03fail = 0, 0

    for line in result:
        X = line.rstrip().split('\t')
        query = X[0] # read identifier
        qX = query.split('_')
        label = qX[0]
        if len(qX) == 5: target = '_'.join(qX[-2:])
        else: target = qX[-1]
        evalue = float(X[10]) # evalue
        pid = float(X[2]) # percent sequence identity
        aLen = int(X[3]) # read alignment length
        qLen = read_lens[query] / 3 # full length of read divided by 3 amino acids
        pMatch = aLen / qLen # percent match length of read length

        if pMatch >= 0.5 and pid >= 90 and evalue <= 1e-3:
            strpass += 1
        else:
            strfail += 1

        if evalue <= 1e-30:
            e30pass += 1
        else:
            e30fail += 1

        if evalue <= 1e-20:
            e20pass += 1
        else:
            e20fail += 1

        if evalue <= 1e-10:
            e10pass += 1
        else:
            e10fail += 1

        if evalue <= 1e-03:
            e03pass += 1
        else:
            e03fail += 1


    data = [
            strpass, strfail, e30pass, e30fail, e20pass,
            e20fail, e10pass, e10fail, e03pass, e03fail
            ]

    return data


def score_static(pabh, nabh, read_lens):

    pos = static_filter_blast(pabh, read_lens)
    neg = static_filter_blast(nabh, read_lens)

    staramer = score_results(pos[0], pos[1], neg[0], neg[1])
    e30 = score_results(pos[2], pos[3], neg[2], neg[3])
    e20 = score_results(pos[4], pos[5], neg[4], neg[5])
    e10 = score_results(pos[6], pos[7], neg[6], neg[7])
    e03 = score_results(pos[8], pos[9], neg[8], neg[9])

    return staramer, e30, e20, e10, e03


def build_bar_plots(df, out):

    fs = 12 # set font size
    metrics = ['F1', 'FPR', 'FNR']
    labels = df.columns.to_list()

    fig, axes = plt.subplots(len(metrics), 1, figsize=(4.25, 10), sharex=True)

    for i, met in enumerate(metrics):
        ax = axes[i]
        data = df.T[met].to_list()
        _ = ax.bar(labels, data, color='#bdbdbd', width=0.8)
        ax.set_ylabel(met, fontsize=fs)
        ax.xaxis.set_tick_params(rotation=45, labelsize=fs)
        ax.set_ylim([0, 1])

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
        '-pp', '--pos_pass_file',
        help='Positive read passing rocker filter (TP).',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-pf', '--pos_fail_file',
        help='Positive read failing rocker filter (FN).',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-np', '--neg_pass_file',
        help='Negative read passing rocker filter (FP).',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-nf', '--neg_fail_file',
        help='Negative read failing rocker filter (TN).',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-pa', '--pos_alignments_file',
        help='Positive read unfiltered alignment file.',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-na', '--neg_alignments_file',
        help='Negative read unfiltered alignment file.',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-prl', '--pos_reads_fasta',
        help='Positive reads fasta file (aligned reads).',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-nrl', '--neg_reads_fasta',
        help='Negative reads fasta file (aligned reads).',
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
    pp = args['pos_pass_file'] # rocker filter TP
    pf = args['pos_fail_file'] # rocker filter FN
    np = args['neg_pass_file'] # rocker filter FP
    nf = args['neg_fail_file'] # rocker filter TN
    # these are to test the static filters
    pa = args['pos_alignments_file'] # unfiltered positive alignments
    na = args['neg_alignments_file'] # unfilter negative alignments
    prl = args['pos_reads_fasta'] # pos reads fasta
    nrl = args['neg_reads_fasta'] # neg reads fasta
    out = args['output_prefix']  # output file prefix

    # Do what you came here to do:
    print('\n\nRunning Script...\n\n')

    print('\t\tBest hit filtering: before, after\n')

    # score ROCker model
    # best hit filter
    TPx = best_hits_filter(pp)
    FNx = best_hits_filter(pf)
    TP, FN, D = dedup_pos_neg(TPx, FNx)
    print(f'\tPassing read duplicates removed from failing reads: {D}\n')
    TPx, FNx = None, None

    FPx = best_hits_filter(np)
    TNx = best_hits_filter(nf)
    FP, TN, D = dedup_pos_neg(FPx, TNx)
    print(f'\tPassing read duplicates removed from failing reads: {D}\n')
    FPx, TNx = None, None

    # then score the results
    rocker = score_results(TP, FN, FP, TN)

    # score static blast filters
    # best hit filter first
    pabh = best_hits_filter(pa)
    nabh = best_hits_filter(na)
    # get read lengths for filtering
    read_lens = get_read_lengths(prl, nrl)


    # score the results
    staramer, e30, e20, e10, e03 = score_static(pabh, nabh, read_lens)

    # output table and pdf
    scores = {
                'ROCker': rocker, 'STARAMR': staramer,
                'evalue 1e-30': e30, 'evalue 1e-20': e20,
                'evalue 1e-10': e10, 'evalue 1e-3': e03,
                }
    rows = [
            'Total', 'P', 'N', 'TP', 'FP', 'TN', 'FN', 'FNR', 'FPR',
            'Sensitivity', 'Specificity', 'Accuracy',
            'Precision', 'Recall', 'F1'
            ]

    df = pd.DataFrame(scores, index=rows)
    df.to_csv(f'{out}_score_table.tsv', sep='\t')

    print('\n\n', df)

    _ = build_bar_plots(df, out)


    print('\n\nComplete success space cadet! Hold on to your boots.\n\n')

if __name__ == "__main__":
    main()
