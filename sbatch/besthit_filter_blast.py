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

import argparse, random, subprocess


def best_hits(query, bitscore, d, line, dups):
    """ Filters the besthit based on bitscore """

    if query in d:
        dups += 1
        old_bitscore = float(d[query][0].split('\t')[11])

        if bitscore > old_bitscore:
            d[query] = [line]

        elif bitscore == old_bitscore:
            d[query].append(line)

    else:
        d[query] = [line]

    return d, dups


def tabular_BlastPlus_filter(infile):

    d = {} # initialize dictionary for bitscore besthits
    dups = 0 # counter for number of duplicate matches for reads
    total = 0 # counter for total blast entries in file

    with open(infile, 'r') as f:

        for l in f:
            total += 1
            X = l.rstrip().split('\t')
            query = X[0] # read identifier
            bitscore = float(X[11]) # bitscore

            d, dups = best_hits(query, bitscore, d, l, dups)

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
    args=vars(parser.parse_args())

    # Do what you came here to do:
    # define input param
    infile = args['in_file']
    print('\n\nRunning Script...\n')
    filtered_best_hits = tabular_BlastPlus_filter(infile)

    # Write output file
    outfile = args['in_file'].split('.')[0] + '.fltrdBstHts.blst'
    with open(outfile, 'w') as o:
        for k,v in filtered_best_hits.items(): o.write(random.choice(v))
        print(
            'Number of best hit entries written to new file:',
            len(filtered_best_hits), '\n\n'
            )

    _ = subprocess.run(['mv', outfile, infile])

if __name__ == "__main__":
    main()
