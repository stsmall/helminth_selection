#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 13 16:43:08 2017

@author: scott t small
"""
import numpy as np
import pandas as pd
from itertools import groupby, combinations
import argparse
from libsequence.polytable import simData
from libsequence.summstats import garudStats
#from libsequence.summstats import polySIM
#from libsequence.summstats import ld
#from libsequence.summstats import lhaf
#from libsequence.windows import Windows

parser = argparse.ArgumentParser()
parser.add_argument('INmsfile', metavar="INmsfile",type=str,help='path to ms-formatted file')
parser.add_argument('-p', '--populations', nargs='+', type=int, required=True, help="population list")
parser.add_argument('-r',"--reps", type=int, required=True, help="number of reps")
parser.add_argument('-u',"--priors", type=str, required=True, help="file with priors")
args = parser.parse_args()

def site_freqspec_fx(gtarray, pos):
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
    nind = gtarray.shape[0]
    #sfs
    sfsfreq, sfscount = np.unique(np.sort(freqsum), return_counts = True)
    sfs = np.zeros(gtarray.shape[0])
    sfs[sfsfreq - 1] = sfscount

    #derived allele freq
    derived = freqsum[len(freqsum)/2.0] / nind
    
    return(sfs,derived)

def jointsfs_fx(gtdict, popiix):
    '''makes the joint site frequency spectrum
    
    '''    
    for rep in gtdict.keys():
        freqsumlist = [np.sum(gtdict[rep][pop]) for pop in popiix]        
        
        joint = [(freqsumlist[i][site],freqsumlist[j][site])
                 for i, j in combinations(range(len(freqsumlist)), 2)
                 for site in range(len(freqsumlist[0]))]         
        jointsfs = [(key,len(list(group))) for key, group in groupby(joint)]
    return(jointsfs)                
                    
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

def calc_otherstats(positions, gtdict, popsizelist):
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
    garudH12 = []
    garudH1 = []
    garudH2 = []
    sfs = []
    hapconfig = []
    derivedfreq = []
    popiix = []
    ix = 0
    pos = .5 #need -Sp 0.5 and -Smark
    for pop in popsizelist:
        popiix.append(range(ix, ix + pop))
        ix += pop
        
    for rep in gtdict.keys():
        derivedfreq_t = []
        sfs_t = np.array([])
        hapconfig_t = []
        garudH12_t = []
        garudH1_t = []
        garudH2_t = []
        for iix in popiix:
            gtarray = gtdict[rep][iix]
            sdpop1 = simData()
            gtpop1 = [''.join(str(n) for n in y) for y in gtarray]
            sdpop1.assign_sep(positions[rep], gtpop1)
            #pspop1 = polySIM(sdpop1)
            #stats
            garudStats_t = (garudStats(sdpop1)) #garud 2015
            garudH12_t.append(garudStats_t[0])
            garudH1_t.append(garudStats_t[1])
            garudH2_t.append(garudStats_t[2])
            #lhaf_t = lhaf(sdpop1,1) #1-HAF is most common; Ronen 2016
            #sfs
            sfsarray, derived = site_freqspec_fx(gtarray, pos)
            np.concatenate(sfs_t, sfsarray)
            #hapconfig
            hapconfig_t.append(haploconfig_fx(gtarray))
            #derived
            derivedfreq_t.append(derived)
        garudH12.append([np.mean(pop) for pop in zip(*garudH12_t)])
        garudH1.append([np.mean(pop) for pop in zip(*garudH1_t)])
        garudH2.append([np.mean(pop) for pop in zip(*garudH2_t)])
        sfs.append(np.sum(sfs_t)/sfs_t.shape[0])
        #collapse all configs to common
        hapconfig.append([(key,len(list(group))) for key, group in groupby(hapconfig_t)])
        #msms
        derivedfreq.append([np.mean(pop) for pop in zip(*derivedfreq_t)])
        #jsfs
        #jointsfs_fx(gtdict, popiix)
    return(garudH12, garudH1, garudH2, sfs, hapconfig, derivedfreq)

def parse_msfile(msfiles_in, popsizelist, number_reps, priors):
    '''calculates stats from ms files

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
    '''
    
    nind = sum(popsizelist)
    garudH12 = []
    garudH1 = []
    garudH2 = []
    hapconfig_summary = []
    sfs = []
    origin_count = []
    derived = []
    for params in range(priors): #for each file
        positions = np.array([]).reshape(number_reps, 0)
        gtdict = {}
        rep = 0
        origin = []
        #read in msfile
        with open("{}.{}.msout".format(msfiles_in, params), 'r') as msstyle:
            #each ms-type file could have 1000s of reps
            for line in msstyle:
                line = line.decode('utf-8')                   
                if line.startswith("positions"):
                    pos = np.array(line.strip().split()[1:],dtype=np.float64)
                    np.concatenate(positions, pos)
                    gt_array = np.zeros((nind, pos.shape[0]), dtype=np.uint8)
                    rep += 1
                    cix = 0
                elif line.startswith("0") or line.startswith("1"): #by size of populations, add to gtdict
                    line = line.decode('utf-8')
                    hap = np.array(list(line.strip()), dtype=np.uint8)
                    gt_array[cix, :] = hap
                    cix += 1
                elif line.startswith("Origin"): #if msms
                    origin.append(int(line.strip().split(":")[1]))
                else: 
                    pass
                gtdict[str(rep)] = gt_array
        
        #run file through stats
        garudH12t, garudH1t, garudH2t, sfst, hapconfigt, derivedt = calc_otherstats(positions, gtdict, popsizelist)
        #collect pop mean for ABC
        garudH12.append(garudH12t)
        garudH1.append(garudH1t)
        garudH2.append(garudH2t)
        sfs.append(sfst)
        hapconfig_summary.append(hapconfigt)
        #unused except for msms
        origin_count.append(np.mean(origin))
        derived.append(derivedt)

    return(garudH12, garudH1, garudH2, sfs, hapconfig_summary, derived, origin_count)

def parse_msstats(msstats_in, priors, demes):
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
    abcstats = np.zeros(priors * demes, 14) #some number
    for f in range(priors):
        msstat = pd.read_table("{}_{}".format(msstats_in, f), header=0, sep='\t')
        meanstat = msstat.groupby('pop').mean().values[2:]
        #add to abcstats
        abcstats.append(meanstat)
    cols = msstat.columns[2:]
    
    return(abcstats, cols)

def format4ABC(priors_in, msfileprefix, popsizelist, number_reps):
    '''creates file for R containing stats of ms or msms runs
    
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
    '''
    priors = pd.read_table("{}".format(priors_in),header=0, sep='\t') 
    lenpriors = priors.shape[0]
    demes = len(popsizelist)
    garudH12, garudH1, garudH2, sfs, hapconfig_summary, derived, origincount = parse_msfile(msfileprefix, 
                                                                                             popsizelist, 
                                                                                             number_reps, 
                                                                                             lenpriors)
    abcstats, cols = parse_msstats(msfileprefix, lenpriors, demes)
    priors = ["theta", "selection"]
    names = priors + "pop" + cols + ["garudH12", "garudH1", "garudH2", "sfs", "hapconfig", "derived", "origins"]
    pop = range(1,demes + 1) * lenpriors
    origin = np.repeat(origincount, demes)
      
    ABCstatTable = pd.DataFrame({np.repeat(priors.theta.values, demes),
                                 pop,
                                 abcstats,
                                 garudH12,
                                 garudH1,
                                 garudH2,
                                 sfs,
                                 hapconfig_summary,
                                 derived,
                                 origincount
                                 })
    ABCstatTable.to_csv("ABCstatTable.csv")
    return(None)
    
if __name__ == '__main__':
    format4ABC(args.priors, args.INmsfile, args.populations, args.reps)    