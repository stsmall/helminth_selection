#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue May  2 13:00:27 2017

@author: scott
"""

import numpy as np
import argparse
from libsequence.polytable import simData
from libsequence.fst import fst

parser = argparse.ArgumentParser()
parser.add_argument('INmsfile', metavar="INmsfile", type=str,
                    help='path to ms-formatted file')
parser.add_argument('-p', '--pops', nargs='+', type=int, required=True,
                    help="population list")
parser.add_argument('-m', "--perms", type=int, required=True,
                    help="number of permutations")
args = parser.parse_args()


def permuteFst(INmsfile, pops, perms):
    """does a permutation test for an ms model, just swaps labels then
    recalculates FST via pylibseq
    """
    nhap = sum(pops)
    fst_t = []
    fst_obs = []
    ix = 0
    popiix = []
    # read inmsfile
    with open(INmsfile, 'r') as msfile:
        for line in msfile:
            line = line.decode('utf-8')
            if line.startswith("positions"):
                pos = np.array(line.strip().split()[1:], dtype=np.float64)
                gt_array = np.zeros((nhap, pos.shape[0]), dtype=np.uint8)
                cix = 0
                while cix < nhap:
                    line = next(msfile)
                    line = line.decode('utf-8')
                    line = list(line.strip())
                    gt_array[cix, :] = np.array(line, dtype=np.uint8)
                    cix += 1
                genonp = gt_array
    # Observed FST
    for p in pops:
        popiix.append(range(ix, ix + p))
        ix += p
    for i, pix in enumerate(popiix):
        for j, jix in enumerate(popiix):
            if i > j:
                popX = genonp[pix]
                popY = genonp[jix]
                sdfst = simData()
                geno_fst = np.vstack([popX, popY])
                gtpop_fst = [''.join(str(n) for n in y) for y in geno_fst]
                sdfst.assign_sep(pos, gtpop_fst)
                size = [popX.shape[0], popY.shape[0]]
                f1 = fst(sdfst, size)
                fst_obs.append(f1.slatkin())
    # FST random permutations
    for p in range(perms):
        popX = genonp[np.random.randint(0, nhap, pops[0])]
        popY = genonp[np.random.randint(0, nhap, pops[0])]
        sdfst = simData()
        geno_fst = np.vstack([popX, popY])
        gtpop_fst = [''.join(str(n) for n in y) for y in geno_fst]
        sdfst.assign_sep(pos, gtpop_fst)
        size = [popX.shape[0], popY.shape[0]]
        f1 = fst(sdfst, size)
        fst_t.append(f1.slatkin())
    # mark significant FST
    fst_tnp = np.array(fst_t)
    Fstdist = [len(np.where(f > fst_tnp)[0]) for f in fst_obs]
    return(fst_obs, [f/float(perms) for f in Fstdist])

if __name__ == '__main__':
    fst_obs, fstpvalue = permuteFst(args.INmsfile, args.pops, args.perms)
    print(fst_obs)
    print(fstpvalue)
