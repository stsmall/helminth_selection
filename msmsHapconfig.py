#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 19 17:31:47 2017

@author: scott
"""
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('INmsfile', metavar="INmsfile", type=str,
                    help='path to ms-formatted file')
parser.add_argument('-p', '--populations', nargs='+', type=int, required=True,
                    help="population list")
parser.add_argument('-r', "--reps", type=int, required=True,
                    help="number of reps")
args = parser.parse_args()


def hapfreqfx(config):
    """
    """
    uniqhaps = np.array([np.array(x) for x in set(tuple(x) for x in config)])
    hapfreq = np.array([len(config[np.all(config == x, axis=1)])
                        for x in uniqhaps], dtype=int)
    return(uniqhaps, hapfreq)


def haploconfig_fx(gtarray, positions, pos):
    """construct haplotype configuration array

    Parameters
    ---------


    Returns
    -------

    """
    # haplotype configuration from gtarray
    gtarrayR =
    Rsite = positions.index(pos)
    for indv in gtarray:
        if indv[Rsite] != 0:
            #resistant
        else:
            #susceptible

    uniqhaps, hapfreq = hapfreqfx(gtarray)
    n = sum(hapfreq)
    C_freq, C_count = np.unique(hapfreq, return_counts=True)
    C = np.zeros(n)
    C[C_freq - 1] = C_count  # full hap configuration
    # haplotype diversity
#    H = 1 - sum([(((i+1)/float(n))**2) * c for i,c in enumerate(C_v)])
#    M = max(C_v)
#    K = sum(C_v)

    # from gtlist
    # K = set(gtlist)

    return(C)


def hapstats(positions, gtdict, popsizelist):
    """calculate stats from ms-type file

    Parameters
    ----------
    poistions : array, float
        matrix of mutational positions
    gtdict : default dict
        dictionary of genotypes per pop; gtdict['pop1'].append((gt_string))

    Returns
    -------
    ms_otherstats : pandas df

    """
    hapconfig_commonR = []
    hapR = []
    hapT = []
    popiix = []
    ix = 0
    pos = 0.5  # need -Sp 0.5 and -Smark
    for pop in popsizelist:
        popiix.append(range(ix, ix + pop))
        ix += pop
    for rep in range(len(gtdict.keys())):
        hapt = []
        hapr = []
        for iix in popiix:
            gtarray = gtdict[str(rep + 1)][iix]
            # hapconfig
            r, t = haploconfig_fx(gtarray, positions[rep], pos)
            hapr.append(r)
            hapt.append(t)
        hapT.append(hapt)
        hapR.append(hapr)
    hap_summT = [np.mean(pop, axis=0) for i, pop in enumerate(zip(*hapT))]
    hap_summR = [np.mean(pop, axis=0) for i, pop in enumerate(zip(*hapR))]
    haparrayT = [np.vstack(i) for i in zip(*hapT)]
    haparrayR = [np.vstack(i) for i in zip(*hapR)]
    for config in haparrayR:
        uniq, hapfreq = hapfreqfx(config)
        hapconfig_commonR.append(zip(uniq, hapfreq))
    return(hap_summR, hapconfig_commonR)


def parse_msfile(msfiles_in, popsizelist, number_reps):
    """
    """
    nind = sum(popsizelist)
    gtdict = {}
    rep = 0
    positions = []
    with open(msfiles_in, 'r') as msstyle:
        # each ms-type file could have 1000s of reps
        for line in msstyle:
            line = line.decode('utf-8')
            if line.startswith("positions"):
                pos = np.array(line.strip().split()[1:], dtype=np.float64)
                positions.append(pos)
                gt_array = np.zeros((nind, pos.shape[0]), dtype=np.uint8)
                rep += 1
                cix = 0
                line = next(msstyle)
                try:
                    while line.startswith("0") or line.startswith("1"):
                        # by size of populations, add to gtdict
                        line = line.decode('utf-8')
                        hap = np.array(list(line.strip()), dtype=np.uint8)
                        gt_array[cix, :] = hap
                        cix += 1
                        line = next(msstyle)
                except StopIteration:
                    pass
                gtdict[str(rep)] = gt_array
            else:
                pass
    hapconfig = hapstats(positions, gtdict, popsizelist)
    return(None)

if __name__ == '__main__':
    parse_msfile(args.INmsfile, args.populations, args.reps)
