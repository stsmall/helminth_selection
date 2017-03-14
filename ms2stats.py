#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 13 16:43:08 2017

@author: scott
"""
import numpy as np
import pandas as pd
from collection import defaultdict
from itertools import groupby
import argparse
from libsequence.polytable import simData
from libsequence.summstats import polySIM
from libsequence.summstats import garudStats
from libsequence.summstats import ld
from libsequence.summstats import lhaf
from libsequence.windows import Windows

parser = argparse.ArgumentParser()
parser.add_argument('INmsfile', metavar="INmsfile",type=str,help='path to ms-formatted file')
parser.add_argument('-p', '--populations', nargs='+', type=int, required=True, help="population list")
parser.add_argument('-r',"--reps", type=int, required=True, help="number of reps")
args = parser.parse_args()

def site_freqspec_fx(gtarray):
    '''calculates the site frequency spectrum. If >1 pop also does jSFS

    Parameters
    ---------

    Returns
    -------

    '''
    #calculate SFS
#    sfs_t = np.zeros(len(gtlist),len(gtlist[0]))
    #gtarray = np.array(map(lambda x: map(int, list(x)), gtlist))
    freqsum = np.sum(gtarray, axis=0)
    sfsfreq, sfscount = np.unique(np.sort(freqsum), return_counts = True)
    sfs = np.zeros(gtarray.shape[0])
    sfs[sfsfreq - 1] = sfscount
    return(sfs)

def haploconfig_fx(gtarray):
    '''construct haplotype configuration array

    Parameters
    ---------


    Returns
    -------

    '''
    #haplotype configuration from gtarray
    uniqhaps = np.array([np.array(x) for x in set(tuple(x) for x in gtarray)])
    hapfreq = np.array([len(gtarray[np.all(gtarray==x, axis=1)]) for x in uniqhaps],dtype=int)
    n = sum(hapfreq)
    C_freq, C_count = np.unique(hapfreq, return_counts=True)
    C = np.zeros(n)
    C[C_freq - 1] = C_count #full hap configuration
    #haplotype diversity
#    H = 1 - sum([(((i+1)/float(n))**2) * c for i,c in enumerate(C_v)])
#    M = max(C_v)
#    K = sum(C_v)

    #from gtlist
    #K = set(gtlist)

    return(C)

def calc_otherstats(positions, gtdict):
    '''calculate stats from ms-type file

    Parameters
    ----------
    poistions : array, float
        matrix of mutational positions
    gtdict : default dict
        dictionary of genotypes per pop; gtdict['pop1'].append((gt_string))

    Returns
    -------
    ms_otherstats : pandas df

    '''
    #numpops = len(gtdict.keys())
    reps = positions.shape[0]
    garudH12 = []
    garudH1 = []
    garudH2 = []
    sfs = []
    hapconfig = []
    for pop in gtdict.keys():
        garudStats_t = []
        sfs_t = np.array([])
        hapconfig_t = []
        for r in range(reps):
            sdpop1 = simData()
            gtarray = gtdict[pop][r]
            gtpop1 = [''.join(str(n) for n in y) for y in gtarray]
            sdpop1.assign_sep(positions[r], gtpop1)
            #pspop1 = polySIM(sdpop1)
            #stats
            garudStats_t.append(garudStats(sdpop1)) #garud 2015
            #lhaf_t = lhaf(sdpop1,1) #1-HAF is most common; Ronen 2016
            #sfs
            sfsarray = site_freqspec_fx(gtarray)
            np.concatenate(sfs_t, sfsarray)
            #hapconfig
            hapconfig_t.append(haploconfig_fx(gtarray))

        garudH12.extend([sum(i[0])/float(r) for i in garudStats_t])
        garudH1.extend([sum(i[1])/float(r) for i in garudStats_t])
        garudH2.extend([sum(i[2])/float(r) for i in garudStats_t])
        sfs.append(np.sum(sfs_t)/sfs_t.shape[0])
        #collapse all configs to common
        hapconfig.append([(key,len(list(group))) for key, group in groupby(hapconfig_t)])

    return(garudH12, garudH1, garudH2, sfs, hapconfig)

def parse_msfile(msfiles_in, popsizelist, number_reps, priors):
    '''calculates stats from ms files

    Parameters
    ----------
    msfiles_in : int
        number of msfiles to loop, prefix of the msfiles, assumes
        "${}_{}.msout".format(prefix, int)
    popsizelist : list, int
        population sizes for each subpop
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
    '''
    demes = len(popsizelist)
    garudH12 = []
    garudH1 = []
    garudH2 = []
    hapconfig_summary = []
    sfs = []
    for params in range(priors): #for each file
        positions = np.array([]).reshape(number_reps, 0)
        gtdict = defaultdict(list)
        rep = 1
        #gtdict['pop1'].append((strings)
        #read in msfile
        with open("{}.{}.msout".format(msfiles_in, params), 'r') as msstyle:
            #each ms-type file could have 1000s of reps
            for line in msstyle:
                if line.startswith("positions"):
                    np.concatenate(positions, map(float,line.strip()))
                    rep += 1
                elif line.startswith("0") or line.startswith("1"): #by size of populations, add to gtdict
                    gtdict[str(rep)].append(np.array(list(line.strip()), dtype=np.uint8))
        #run file through stats
        garudH12t, garudH1t, garudH2t, sfst, hapconfigt = calc_otherstats(positions, gtdict)

        #collect pop mean for ABC
        garudH12.append(garudH12t)
        garudH1.append(garudH1t)
        garudH2.append(garudH2t)
        sfs.append(sfst)
        hapconfig_summary.append(hapconfigt)

    haploconfig = pd.DataFrame({"priors" : np.repeat([1, priors], demes),
                             "pop" : range(1, demes + 1) * priors,
                             "garudH12" : garudH12,
                             "garudH1" : garudH1,
                             "garudH2" : garudH2,
                             "sfs" : sfs,
                             "hapconfig" : hapconfig_summary})

    haploconfig = haploconfig.loc[:, ["priors", "pop", "garudH12","garudH1", "garudH2", "sfs", "hapconfig"]]

    return(haploconfig)


def parse_msstats_out(msstats_in, popsizelist, number_reps):
    '''parses file from msstats into pops and means values

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

    '''




def format4ABC(priors_in, msstats_out, popsizelist):
    '''
    popsizelist = [10, 12, 15, 8, 10, 12]
    '''

    #read in prior file


    #read in msstats


    #summarry of msstats per pop
    #read in a file to pd.DF
    #mean of stat by population


    '''
    1.2 0.05 1 -2.5 [([5, 0, 1, 0], 10), ([5, 0, 2, 1, 1], 5)] [0, 0.23, 0.35]
    1.2 0.05 2 -1.3 [([5, 0, 1, 0], 10), ([5, 0, 2, 1, 1], 5)] [0.23, 0, 0.35]
    1.2 0.05 3 -1.2 [([5, 0, 1, 0], 10), ([5, 0, 2, 1, 1], 5)] [0.23, 0.35, 0]
    1.2 0.05 4
    1.2 0.05 5
    1.2 0.05 6
    '''

    ABCstatTable = pd.DataFrame({"prior1" : [],
                                 "prior2" : [],
                                 "pop1" : [],
                                 "pop2" : [],
                                 "stat1" : [],
                                 "stat2" : [],
                                 "hapconfig" : [],
                                 "FST_p" : []
                                 })
    ABCstatTable = ABCstatTable.loc[:, ['prior1', 'prior2', 'pop1', 'pop2', 'stat1', 'stat2', 'hapconfig', 'FST_p']]
    ABCstatTable.to_csv("ABCstatTable.csv")

    return(None)