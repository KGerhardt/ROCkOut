#!/usr/bin/env python

'''Filter and score tabular blast file

This script is intended to work with reads that have been labeled by
"04a_relabel_reads_fasta.py" before the blast search.

This script will use a static filter on the blast results and report
the F1, accuracy, recall, precision, TP, TN, FP, FN scores in a tsv file.

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: January 2023
License :: GNU GPLv3
Copyright 2023 Roth Conrad
All rights reserved
-------------------------------------------
'''

import argparse

# not using best hits filter currently because as far as i can tell
# siyu/brittany did not and staramr does not
def best_hits(query, bitscore, d, line, dups):
    """ Filters the besthit based on bitscore """

    if query in d:
        dups += 1
        old_bitscore = float(d[query].split('\t')[11])

        if bitscore > old_bitscore:
            d[query] = line

    else:
        d[query] = line

    return d, dups


def score_blast(infile, pml, psi, evl):
    # filter and score the blast results

    TP, FP, TN, FN = 0, 0, 0, 0

    with open(infile, 'r') as f:
        for l in f:
            X = l.rstrip().split('\t')
            query = X[0] # read identifier
            qX = query.split('_')
            label = qX[0]
            if len(qX) == 5: target = '_'.join(qX[-2:])
            else: target = qX[-1]
            evalue = float(X[10]) # evalue
            pid = float(X[2]) # percent sequence identity
            aLen = int(X[3]) # read alignment length
            qLen = int(X[12]) / 3 # full length of read divided by 3 amino acids
            pMatch = aLen / qLen # percent match length of read length

            if pMatch >= pml and pid >= psi and evalue <= evl:
                if label == "positive" and target == "Target": TP += 1
                elif label == "positive" and target == "Non_Target": FP += 1
                elif label == "negative": FP += 1
                else: print("positive result error:\n", l)

            else:
                if label == "negative": TN += 1
                elif label == "positive" and target == "Non_Target": TN += 1
                elif label == "positive" and target == "Target": FN += 1
                else: print("negative result error:\n", l)

                
    precision = TP / (TP + FP)
    recall = TP / (TP + FN)
    F1 = 2 * (precision * recall) / (precision + recall)

    print(f"TP: {TP}")
    print(f"FP: {FP}")
    print(f"TN: {TN}")
    print(f"FN: {FN}")
    print(f"FPR: {FP / (TP + FP)}")
    print(f"FNR: {FN / (TP + TN)}")
    print(f"precision: {precision}")
    print(f"recall: {recall}")
    print(f"F1: {F1}")

    print(f"Pos Total: {TP + FP}")
    print(f"Neg Total: {TN + FN}")
    print(f"Total: {TP + FP + TN + FN}")
    print(f"Line Total: {total_lines}")

def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-i', '--in_file',
        help='Please specify the tabular blast input file name!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-pml', '--percent_match_length',
        help='(OPTIONAL) Percent match length to filter for (Default=0.5).',
        metavar='',
        type=float,
        required=False,
        default=0.5
        )
    parser.add_argument(
        '-psi', '--percent_sequence_identity',
        help='(OPTIONAL) Percent sequence identity to filter for (Default=90).',
        metavar='',
        type=float,
        required=False,
        default=90.0
        )
    parser.add_argument(
        '-evl', '--e_value',
        help='(OPTIONAL) e-value to filter for (Default=0.001).',
        metavar='',
        type=float,
        required=False,
        default=0.001
        )
    args=vars(parser.parse_args())

    # Do what you came here to do:
    print('Running Script...')

    score_blast(
        args['in_file'],
        args['percent_match_length'],
        args['percent_sequence_identity'],
        args['e_value']
        )




if __name__ == "__main__":
    main()
