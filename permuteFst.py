#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue May  2 13:00:27 2017

@author: scott
"""
from __future__ import print_function

import numpy as np
import argparse
import random
import subprocess
import pandas as pd
from libsequence.polytable import simData
from libsequence.fst import fst

parser = argparse.ArgumentParser()
parser.add_argument('-m', "--msfile", default=None, type=str,
                    help='path to ms-formatted file')
parser.add_argument('-p', '--pops', nargs='+', type=int, required=True,
                    help="population list")
parser.add_argument('-np', "--n_perm", type=int, default=1000,
                    help="number of permutations")
parser.add_argument('-N', '--effectivesize', type=int, default=5E5,
                    help='effective population size')
parser.add_argument('-r', "--reps", type=int, default=1000,
                    help="number of reps")
parser.add_argument("--perm", action="store_true",
                    help="run permutation test")
args = parser.parse_args()


def readms(msfile, pops):
    """does a permutation test for an ms model, just swaps labels then
    recalculates FST via pylibseq
    """
    gtdict = {}
    posdict = {}
    nhap = sum(pops)
    rep = -1
    # read inmsfile
    with open(msfile, 'r') as msfile:
        for line in msfile:
            line = line.decode('utf-8')
            if line.startswith('//'):
                rep += 1
            elif line.startswith("positions"):
                pos = np.array(line.strip().split()[1:], dtype=np.float64)
                gt_array = np.zeros((nhap, pos.shape[0]), dtype=np.uint8)
                cix = 0
                while cix < nhap:
                    line = next(msfile)
                    line = line.decode('utf-8')
                    line = list(line.strip())
                    gt_array[cix, :] = np.array(line, dtype=np.uint8)
                    cix += 1
                gtdict[str(rep)] = gt_array
                posdict[str(rep)] = pos
    return(posdict, gtdict)


def readms2(msout, pops):
    """does a permutation test for an ms model, just swaps labels then
    recalculates FST via pylibseq
    """
    gtdict = {}
    posdict = {}
    nhap = sum(pops)
    rep = -1
    # read ms from stout
    msfile = iter(msout.stdout.readline, '')
    for line in msfile:
        line = line.decode('utf-8')
        if line.startswith('//'):
            rep += 1
        elif line.startswith("positions"):
            pos = np.array(line.strip().split()[1:], dtype=np.float64)
            gt_array = np.zeros((nhap, pos.shape[0]), dtype=np.uint8)
            cix = 0
            while cix < nhap:
                line = next(msfile)
                line = line.decode('utf-8')
                line = list(line.strip())
                gt_array[cix, :] = np.array(line, dtype=np.uint8)
                cix += 1
            gtdict[str(rep)] = gt_array
            posdict[str(rep)] = pos
    return(posdict, gtdict)


def calcfst(pops, posdict, gtdict):
    """
    """
    fst_obs = []
    ix = 0
    popiix = []
    pw = len(pops)
    fstarray = np.zeros([len(posdict.keys()), (pw*(pw-1))/2])
    # Observed FST
    for p in pops:
        popiix.append(range(ix, ix + p))
        ix += p
    for r in gtdict.keys():
        fst_obs = []
        for i, pix in enumerate(popiix):
            for j, jix in enumerate(popiix):
                if i > j:
                    popX = gtdict[r][pix]
                    popY = gtdict[r][jix]
                    sdfst = simData()
                    geno_fst = np.vstack([popX, popY])
                    gtpop_fst = [''.join(str(n) for n in y) for y in geno_fst]
                    sdfst.assign_sep(posdict[r], gtpop_fst)
                    size = [popX.shape[0], popY.shape[0]]
                    f1 = fst(sdfst, size)
                    fst_obs.append(f1.slatkin())
        fstarray[int(r), :] = fst_obs
    return(fstarray)


def permtest(gtdict, posdict, pops, n_perm, fstarray):
    """
    """
    nhap = sum(pops)
    fst_t = []
    r = random.choice(gtdict.keys())
    # FST random permutations
    for p in range(n_perm):
        popX = gtdict[r][np.random.randint(0, nhap, pops[0])]
        popY = gtdict[r][np.random.randint(0, nhap, pops[0])]
        sdfst = simData()
        geno_fst = np.vstack([popX, popY])
        gtpop_fst = [''.join(str(n) for n in y) for y in geno_fst]
        sdfst.assign_sep(posdict[r], gtpop_fst)
        size = [popX.shape[0], popY.shape[0]]
        f1 = fst(sdfst, size)
        fst_t.append(f1.slatkin())
    # mark significant FST
    fst_tnp = np.array(fst_t)
    Fstdist = [len(np.where(f > fst_tnp)[0]) for f in fstarray]
    return([1-(f/float(n_perm)) for f in Fstdist])


def runmssims(ms, Ne, migp, pops, reps, theta, gens, time):
    """
    """
    nhap = sum(pops)
    demes = len(pops)
    fsts = []
    for m in migp:
        for t in time:
            ms_params = {
                        'ms': ms,
                        'nhaps': nhap,
                        'nreps': reps,
                        'theta': theta,
                        'demes': "{} {}".format(demes,
                                                " ".join(map(str, pops))),
                        'Nm': m * 4 * Ne,
                        'time': (gens * t) / (4.0 * Ne)}
            msms_base = ("{ms} {nhaps} {nreps} -t {theta} "
                         "-I {demes} {Nm} -ej "
                         "{time} 2 1 -ej {time} 3 1 "
                         "-ej {time} 4 1 -ej {time} 5 1 ")
            mscmd = msms_base.format(**ms_params)
            print(mscmd)
            msout = subprocess.Popen(mscmd, shell=True, stdout=subprocess.PIPE)
            # parse
            posdict, gtdict = readms2(msout, pops)
            # calc FST
            fstarray = calcfst(pops, posdict, gtdict)
            fsts.append(np.mean(fstarray, axis=1))
#    import ipdb;ipdb.set_trace()
    dfFig1a = pd.DataFrame({'mig': np.repeat(migp, len(time) * reps),
                            'time': list(np.repeat(time, reps)) * len(migp),
                            'fst': np.concatenate(fsts).ravel()
                            })
    dfFig1a = dfFig1a.loc[:, ['time', 'mig', 'fst']]
    dfFig1a.to_csv("Fig1A_helminth.csv")
    return(None)


if __name__ == '__main__':
    pops = args.pops  # list
    reps = args.reps  # default 1000
    Ne = args.effectivesize  # default 1E6
    ms = '/usr/bin/ms'  # default path
    gens = 12  # gens per year
    time = [20, 50, 100]  # time in years
    migp = [0, 0.0001, 0.001, 0.01]  # migration proportion
    theta = 8.28  # UK 8.28, India 6.80, France 3.6, China 4.76, theta
    if args.msfile is not None:
        posdict, gtdict = readms(args.msfile, pops)
        fstarray = calcfst(pops, posdict, gtdict)
        a = np.round(np.mean(fstarray, axis=0), 2)
        print(a)
        if args.perm:
            fstpvalue = permtest(gtdict, posdict, pops, args.n_perm, a)
            print("[%s]" % ", ".join(map(str, fstpvalue)))
    else:
        runmssims(ms, Ne, migp, pops, reps, theta, gens, time)
