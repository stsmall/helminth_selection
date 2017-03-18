#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 14 16:32:37 2017

@author: scott
"""

import numpy as np
import pickle
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('INfasta', metavar="INfasta", type=str,
                    help='path to fasta file')
parser.add_argument('-s', '--seqs', type=int, required=True,
                    help="number of sequences")
parser.add_argument('-l', "--seqlength", type=int, required=True,
                    help="length of sequence")
parser.add_argument('-p', '--popsizelist', nargs='+', type=int, required=True,
                    help="population list")
args = parser.parse_args()


class MSgenotype(object):
    def __init__(self, haplotype=[], positions=[], names=[], popsizelist=[]):

        self.haps = haplotype
        self.pos = positions
        self.names = names
        self.pops = popsizelist


def fasta2pylibseq(fastaFile, numberseqs, lensequence, popsizelist):
    '''take a fasta alignmnet with out group as first sequence and convert it
    to ms-type format.

    Parameters
    ----------
    fasta : file
        alignment in fasta format
    numberseq : int
        number of sequences to expect for array size not including reference
    lensequence : int
        length of sequence
    Returns
    -------
    gtbinary : object
        genotype object of genotypes in strings
    positions : array
        list of mutation positions 0 - 1
    '''
    seqnames = []
    positions = np.arange(0.0, 1.0, 1.0/lensequence)
    genotypes = np.zeros(shape=[numberseqs, len(positions)])
    tmpfasta = []
    with open("{}".format(fastaFile), 'r') as fasta_aln:
        for line in fasta_aln:
            if line.startswith(">"):
                seqnames.append(line[1:].strip())
            else:
                tmpfasta.append(line.strip())

    refseq = tmpfasta[0]
    for ind, seq in enumerate(tmpfasta[1:]):
        try:
            muts = [i for i in range(len(seq)) if refseq[i] != seq[i]]
        except IndexError:
            pass
        genotypes[ind][muts] = 1

    fastaname = MSgenotype(genotypes, positions, seqnames[1:], popsizelist)

    with open('{}.pkl'.format('fastaout'), 'wb') as output:
        pickler = pickle.Pickler(output, -1)
        fastaout = MSgenotype(fastaname.haps, fastaname.pos, fastaname.names,
                              fastaname.pops)
        pickler.dump(fastaout)
    del fastaout

    return(None)

if __name__ == '__main__':
    fasta2pylibseq(args.INfasta, args.seqs, args.seqlength, args.popsizelist)
