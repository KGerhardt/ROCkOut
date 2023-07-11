#!/usr/bin/env python

'''Filters Blast+ Tabular Output for best hit.

Blast frequently returns multiple matches per query read.
This script uses bitscore to select the best hit.

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: Wednesday, August 28th, 2019
License :: GNU GPLv3
Copyright 2019 Roth Conrad
All rights reserved
-------------------------------------------
'''

import argparse, random


def best_hits(query, bitscore, d, entry, dups):
    """ Filters the besthit based on bitscore """

    if query in d:
        dups += 1
        old_bitscore = float(d[query][0][1])
        if bitscore > old_bitscore:
            d[query] = [entry]

        elif bitscore == old_bitscore:
            d[query].append(entry)

    else:
        d[query] = [entry]

    return d, dups


def tabular_BlastPlus_filter(infile):

    d = {} # initialize dictionary for bitscore besthits
    dups = 0 # counter for number of duplicate matches for reads
    total = 0 # counter for total blast entries in file

    with open(infile, 'r') as f:
        # skip the hmm header
        for _ in range(14):
            next(f)
        # read throught the results
        for l in f:
            # remove whitespace and split fields
            X = l.lstrip().rstrip().split()
            # break the loop before the domain annotation output
            if X[1] == 'inclusion': break
            # sort out results
            total += 1 # add total
            # define fields
            query = '_'.join(X[8].split('_')[1:])
            bitscore = float(X[1])
            # track duplicate hits 
            d, dups = best_hits(query, bitscore, d, X, dups)

    print('Total number of entries in blast file:', total)
    print('Number of duplicate blast matches:', dups)

    return d


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
        '-o', '--out_file',
        help='Please specify the output file name!',
        metavar='',
        type=str,
        required=True
        )
    args=vars(parser.parse_args())

    # Do what you came here to do:
    # define input param
    infile = args['in_file']
    outfile = args['out_file']
    print('\n\nRunning Script...\n')
    filtered_best_hits = tabular_BlastPlus_filter(infile)

    # Write output file
    with open(outfile, 'w') as o:
        for k,v in filtered_best_hits.items():
            X = random.choice(v)
            # Fields are query, evaule, score, bias
            lineout = f'{k}\t{X[0]}\t{X[1]}\t{X[2]}\n'
            o.write(lineout)
        print(
            'Number of best hit entries written to new file:',
            len(filtered_best_hits), '\n\n'
            )

if __name__ == "__main__":
    main()
