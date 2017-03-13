#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 13 16:43:08 2017

@author: scott
"""
import pandas as pd
import numpy as np


def msout2stats(msfiles_in, popsizelist):
    '''calculates stats from ms files

    Parameters
    ----------
    msfiles_in : int
        number of msfiles to loop
    popsizelist : list, int
        population sizes for each subpop

    Returns
    ------

    '''
    #read in msfile
    with open("ms.{}".format(rep), 'r') as msstyle:


    #calculate stats in pylibseq




    #calculate haplotype config

    return(hapconfig, libseqstats)

def formatms4ABC(priors_in, msstats_out, popsizelist):
    '''
    popsizelist = [10, 12, 15, 8, 10, 12]
    '''

    #read in prior file


    #read in msstats


    #summarry of msstats per pop
    #read in a file to pd.DF
    #mean of stat by population


    hapconfig, libseqstats = msout2stats(msfiles, popsizelist)

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
                                 "hapconfig" : [].
                                 "FST_p" : []
                                 })
    ABCstatTable = ABCstatTable.loc[:, ['prior1', 'prior2', 'pop1', 'pop2', 'stat1', 'stat2', 'hapconfig', 'FST_p']]
    ABCstatTable.to_csv("ABCstatTable.csv")

    return(None)