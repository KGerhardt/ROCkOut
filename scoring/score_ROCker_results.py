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

def score_rocker(passing, failing):
    # filter and score the blast results

    TP, FP, TN, FN = 0, 0, 0, 0
    total_lines = 0

    with open(passing, 'r') as f:
        for l in f:
            total_lines += 1
            X = l.rstrip().split('\t')
            query = X[0] # read identifier
            qX = query.split('_')
            label = qX[0]
            if len(qX) == 5: target = '_'.join(qX[-2:])
            else: target = qX[-1]

            if label == "positive" and target == "Target": TP += 1
            elif label == "positive" and target == "Non_Target": FP += 1
            elif label == "negative": FP += 1
            else: print("Passing error:", l)

    with open(failing, 'r') as f:
        for l in f:
            total_lines += 1
            X = l.rstrip().split('\t')
            query = X[0] # read identifier
            qX = query.split('_')
            label = qX[0]
            if len(qX) == 5: target = '_'.join(qX[-2:])
            else: target = qX[-1]

            if label == "positive" and target == "Target": FN += 1
            elif label == "positive" and target == "Non_Target": TN += 1
            elif label == "negative": TN += 1
            else: print("Failing error:", l)

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
        '-p', '--passing_file',
        help='Please specify the rocker filter passing file!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-f', '--failing_file',
        help='Please specify the rocker filter failing file!',
        metavar='',
        type=str,
        required=True
        )
    args=vars(parser.parse_args())

    # Do what you came here to do:
    print('Running Script...')

    score_rocker(
        args['passing_file'],
        args['failing_file']
        )

if __name__ == "__main__":
    main()
