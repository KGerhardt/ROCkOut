#!/usr/bin/env python

''' Translate nucleotide sequence to its 6 frames of amino acid sequence

This script returns an amino acid sequence fasta with 6 entries per
nucleotide fasta input. It uses the biopython seq module to translate
all six reading frames of a nucleotide sequence.

This script was written and used for the ROCkOut project to translate
simulated illumina reads to amino acid sequence similar to blastx. The
output is used with hmmsearch to score the models.

Requires BioPython:
conda install -c conda-forge biopython
or
pip install biopython

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: June 2023
License :: GNU GPLv3
Copyright 2023 Roth Conrad
All rights reserved
-------------------------------------------
'''

import argparse
from Bio.Seq import Seq

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

def translate_reading_frames(infile, outfile):

    with open(infile, 'r') as inf, open(outfile, 'w') as outf:
        for name, seq in read_fasta(inf):

            f1 = str(Seq(seq[0:]).translate())
            f2 = str(Seq(seq[1:]).translate())
            f3 = str(Seq(seq[2:]).translate())
            revseq = str(Seq(seq).reverse_complement())
            f4 = str(Seq(revseq[0:]).translate())
            f5 = str(Seq(revseq[1:]).translate())
            f6 = str(Seq(revseq[2:]).translate())

            n = name[1:]
            outstr = (
                        f'>f1_{n}\n{f1}\n'
                        f'>f2_{n}\n{f2}\n'
                        f'>f3_{n}\n{f3}\n'
                        f'>f4_{n}\n{f4}\n'
                        f'>f5_{n}\n{f5}\n'
                        f'>f6_{n}\n{f6}\n'
                        )

            outf.write(outstr)


def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-i', '--input_file',
        help='Please specify the nucleotide fasta file!',
        metavar=':',
        type=str,
        required=True
        )
    parser.add_argument(
        '-o', '--output_file',
        help='Please specify output amino acid fasta file!',
        metavar=':',
        type=str,
        required=True
        )
    args=vars(parser.parse_args())

    # define input params
    infile = args['input_file']
    outfile = args['output_file']

    # Do what you came here to do:
    print('\n\nRunning Script...\n\n')

    _ = translate_reading_frames(infile, outfile)

    print('\n\nComplete success space cadet! Hold on to your boots.\n\n')
    
if __name__ == "__main__":
    main()

