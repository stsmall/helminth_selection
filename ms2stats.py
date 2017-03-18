#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 13 16:43:08 2017

@author: scott t small
"""
import pickle
import numpy as np
import pandas as pd
from itertools import groupby, combinations
import argparse
from libsequence.polytable import simData
from libsequence.summstats import garudStats
# from libsequence.summstats import polySIM
# from libsequence.summstats import ld
# from libsequence.summstats import lhaf
# from libsequence.windows import Windows

parser = argparse.ArgumentParser()
parser.add_argument('INmsfile', metavar="INmsfile", type=str,
                    help='path to ms-formatted file')
parser.add_argument('-p', '--populations', nargs='+', type=int, required=True,
                    help="population list")
parser.add_argument('-r', "--reps", type=int, required=True,
                    help="number of reps")
parser.add_argument('-u', "--priors", type=str, required=True,
                    help="file with priors")
args = parser.parse_args()


def site_freqspec_fx(gtarray, pos):
    """calculates the site frequency spectrum. If >1 pop also does jSFS

    Parameters
    ---------

    Returns
    -------

    """
    # calculate SFS
    freqsum = np.sum(gtarray, axis=0)
    nind = gtarray.shape[0]
    # sfs
    sfsfreq, sfscount = np.unique(np.sort(freqsum), return_counts=True)
    sfs = np.zeros(gtarray.shape[0])
    sfs[sfsfreq - 1] = sfscount

    # derived allele freq
    derived = freqsum[np.int(len(freqsum)/(1.0/pos))] / nind

    return(sfs, derived)


def hapconfigfx(config):
    """
    """
    uniqhaps = np.array([np.array(x) for x in set(tuple(x) for x in config)])
    hapfreq = np.array([len(config[np.all(config == x, axis=1)])
                        for x in uniqhaps], dtype=int)
    return(uniqhaps, hapfreq)


def jointsfs_fx(gtdict, popiix):
    """makes the joint site frequency spectrum

    """
    for rep in gtdict.keys():
        freqsumlist = [np.sum(gtdict[rep][pop]) for pop in popiix]

        joint = [(freqsumlist[i][site], freqsumlist[j][site])
                 for i, j in combinations(range(len(freqsumlist)), 2)
                 for site in range(len(freqsumlist[0]))]
        jointsfs = [(key, len(list(group))) for key, group in groupby(joint)]
    return(jointsfs)


def haploconfig_fx(gtarray):
    """construct haplotype configuration array

    Parameters
    ---------


    Returns
    -------

    """
    # haplotype configuration from gtarray
    uniqhaps, hapfreq = hapconfigfx(gtarray)
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


def calc_otherstats(positions, gtdict, popsizelist):
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
    hapconfig_common = []
    gH12 = []
    gH1 = []
    gH21 = []
    sfs = []
    hap = []
    dfreq = []
    popiix = []
    ix = 0
    pos = 0.5  # need -Sp 0.5 and -Smark
    for pop in popsizelist:
        popiix.append(range(ix, ix + pop))
        ix += pop
    for rep in range(len(gtdict.keys())):
        dfreqt = []
        sfst = []
        hapt = []
        gH12t = []
        gH1t = []
        gH21t = []
        for iix in popiix:
            gtarray = gtdict[str(rep+1)][iix]
            sdpop1 = simData()
            gtpop1 = [''.join(str(n) for n in y) for y in gtarray]
            sdpop1.assign_sep(positions[rep], gtpop1)
            # pspop1 = polySIM(sdpop1)
            # stats
            garudStats_t = (garudStats(sdpop1))  # garud 2015
            gH12t.append(garudStats_t['H12'])
            gH1t.append(garudStats_t['H1'])
            gH21t.append(garudStats_t['H2H1'])
            # lhaf_t = lhaf(sdpop1,1) #1-HAF is most common; Ronen 2016
            # sfs
            sfsarray, derived = site_freqspec_fx(gtarray, pos)
            sfst.append(sfsarray)
            # hapconfig
            hapt.append(haploconfig_fx(gtarray))
            # derived
            dfreqt.append(derived)
        gH12.append(gH12t)
        gH1.append(gH1t)
        gH21.append(gH21t)
        sfs.append(sfst)
        hap.append(hapt)
        dfreq.append(dfreqt)
    garudH12 = [np.mean(pop) for pop in zip(*gH12)]
    garudH1 = [np.mean(pop) for pop in zip(*gH1)]
    garudH21 = [np.mean(pop) for pop in zip(*gH21)]
    sfs_summ = [np.mean(pop, axis=0)/(popsizelist[i])
                for i, pop in enumerate(zip(*sfs))]
    hap_summ = [np.mean(pop, axis=0) for i, pop in enumerate(zip(*hap))]
    haparray = [np.vstack(i) for i in zip(*hap)]
    for config in haparray:
        uniq, hapfreq = hapconfigfx(config)
        hapconfig_common.append(zip(uniq, hapfreq))
    derivedfreq = [np.mean(pop) for pop in zip(*dfreq)]
    # jointsfs_fx(gtdict, popiix)
    return(garudH12, garudH1, garudH21, sfs_summ, hap_summ, derivedfreq)


def parse_msfile(msfiles_in, popsizelist, number_reps, priors):
    """calculates stats from ms files

    Parameters
    ----------
    msfiles_in : int
        number of msfiles to loop, prefix of the msfiles, assumes
        "${}_{}.msout".format(prefix, int)
    popsizelist : list, int
        population sizes for each subpop, if
    number_reps : int
        number of reps in each file
    params : int
        total number of files, combinations

    Returns
    ------
    positions : array
        array of positions
    genotype : array
        array of genotypess
    """

    nind = sum(popsizelist)
    garudH12 = []
    garudH1 = []
    garudH21 = []
    hapconfig_summary = []
    sfs = []
    origin_count = []
    derived = []
    positions = []
    for params in range(priors):  # for each file
        gtdict = {}
        rep = 0
        origin = []
        # read in msfile
        with open("{}_{}.msout".format("out.sm", 1), 'r') as msstyle:
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
                elif line.startswith("Origin"):  # if msms
                    origin.append(int(line.strip().split(":")[1]))
                else:
                    pass

        # run file through stats
        gH12t, gH1t, gH21t, sfst, hapt, dert = calc_otherstats(positions,
                                                               gtdict,
                                                               popsizelist)
        # collect pop mean for ABC
        garudH12.append(gH12t)  # extend?
        garudH1.append(gH1t)
        garudH21.append(gH21t)
        sfs.append(sfst)
        hapconfig_summary.append(hapt)
        # unused except for msms
        origin_count.append(np.mean(origin))
        derived.append(dert)

    return(garudH12, garudH1, garudH21, sfs, hapconfig_summary, derived,
           origin_count)


def parse_msstats(msstats_in, priors, demes):
    """parses file from msstats into pops and means values

    Parameters
    ----------
    msfiles_in : int
        number of msfiles to loop, prefix of the msfiles, assumes
        "${}_{}.msout".format(prefix, int)
    popsizelist : list, int
        population sizes for each subpop
    number_reps : int
        number of reps in each file

    Returns
    -------

    """
    abcstats = np.zeros(priors * demes, 14)  # some number
    for f in range(priors):
        msstat = pd.read_table("{}_{}".format(msstats_in, f), header=0,
                               sep='\t')
        meanstat = msstat.groupby('pop').mean().values[2:]
        # add to abcstats
        abcstats.append(meanstat)
    cols = msstat.columns[2:]

    return(abcstats, cols)


def format4ABC(priors_in, msfileprefix, popsizelist, number_reps):
    """creates file for R containing stats of ms or msms runs

    Parameters
    ----------
    popsizelist : int, list
        [10, 12, 15, 8, 10, 12] list of population sizes
    priors_in : file
        priors file
    msfileprefix : str
        files prefix for ms and msstats
    number_reps : int
        number of reps in each ms file

    Returns
    -------
    ABCstatTable : pandas df object
    """
    priors = pd.read_table("{}".format(priors_in), header=0, sep='\t')
    lenpriors = priors.shape[0]
    demes = len(popsizelist)
    gH12, gH1, gH21, sfs, hap, der, orig = parse_msfile(msfileprefix,
                                                        popsizelist,
                                                        number_reps,
                                                        lenpriors)
    abcstats, cols = parse_msstats(msfileprefix, lenpriors, demes)
    priors = ["theta", "selection"]
    names = priors + "pop" + cols + ["garudH12", "garudH1", "garudH21", "sfs",
                                     "hapconfig", "derived", "origins"]
    pop = range(1, demes + 1) * lenpriors
    origin = np.repeat(orig, demes)

    ABCstatTable = pd.DataFrame({np.repeat(priors.theta.values, demes),
                                 pop,
                                 abcstats,
                                 gH12,
                                 gH1,
                                 gH21,
                                 sfs,
                                 hap,
                                 der,
                                 origin
                                 })
    ABCstatTable.to_csv("ABCstatTable.csv")
    return(None)

if __name__ == '__main__':
    format4ABC(args.priors, args.INmsfile, args.populations, args.reps)
