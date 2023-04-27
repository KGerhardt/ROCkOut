#!/usr/bin/env python

''' Sort simulated reads from ROCkOut from pos and neg IDs

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

def sort_pos_neg(pos, neg):

    # store the read name from the positive target sequences
    # do not write neg sequences that match positive target read names
    # will also check read names to not write duplicates
    pos_dict = {}
    neg_dict = {}

    # open output files
    pos_out = open('01_simulated_raw_reads_pos.fasta', 'w')
    neg_out = open('02_simulated_raw_reads_neg.fasta', 'w')

    # read the positive file. Reads labeled as "Target" are positives
    # reads labeled as "Non_Target" are negatives
    # reads labeled as "Foreign_Target" belong to another protein ID
    # and they are ignored for our purposes here.
    with open(pos, 'r') as file:
        for name, seq in read_fasta(file):
            X = name.split(';')
            label = X[-1]
            newName = ';'.join(X[:-1])

            # if we already have the read, do not write it to file
            if newName in pos_dict or newName in neg_dict: continue

            # if we don't have the read, let's sort it.
            if label == "Target":
                pos_dict[newName] = ''
                pos_out.write(f'{name}\n{seq}\n')
            elif label == "Non_Target":
                neg_dict[newName] = ''
                neg_out.write(f'{name}\n{seq}\n')

    # check for optional neg file to continue
    if neg:
        # process neg file
        with open(neg, 'r') as file:
            for name, seq in read_fasta(file):
                X = name.split(';')
                label = X[-1]
                newName = ';'.join(X[:-1])

                # if we already have the read, do not write it to file
                if newName in pos_dict or newName in neg_dict: continue

                # if we don't have the read, let's sort it.
                # all reads from neg not labeled as 'Foreign_Target' get 
                # written to neg_out
                if label != "Foreign_Target":
                    neg_dict[newName] = ''
                    neg_out.write(f'{name}\n{seq}\n')

    return True

def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-pos', '--reads_from_positive',
        help='Reads fasta from positive IDs!',
        metavar=':',
        type=str,
        required=True
        )
    parser.add_argument(
        '-neg', '--reads_from_negative',
        help='Reads fasta from negative IDs!',
        metavar=':',
        type=str,
        required=False
        )
    args=vars(parser.parse_args())

    # Define inputs
    pos = args['reads_from_positive']
    neg = args['reads_from_negative']

    # sort the buggers
    _ = sort_pos_neg(pos, neg)

    
if __name__ == "__main__":
    main()

