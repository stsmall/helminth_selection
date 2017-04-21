#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 19 17:31:47 2017

@author: scott
"""
import numpy as np
import argparse
import subprocess
import pandas as pd
import bisect

parser = argparse.ArgumentParser()
parser.add_argument('-N', '--effectivesize', type=int,
                    help='effective population size', required=True)
parser.add_argument('-p', '--populations', nargs='+', type=int, required=True,
                    help="population list")
parser.add_argument('-r', "--reps", type=int, required=True,
                    help="number of reps")
args = parser.parse_args()


def parse_msfile(msout, nhap, num_reps):
    """
    """
    gtdict = {}
    posdict = {}
    origcount = np.zeros(num_reps)
    freqtrace = {}
    rep = -1
    msmsfile = iter(msout.stdout.readline, '')
    for line in msmsfile:
        line = line.decode('utf-8')
        if line.startswith("//"):
            rep += 1
            trace = []
        elif line.startwith("Frequency"):
            line = next(msmsfile)
            while not line.startswith("segsites"):
                trace.append(list(line.strip()))
                line = next(msmsfile)
            freqtrace[str(rep)] = np.vstack((trace))
        elif line.startswith("positions"):
            pos = np.array(line.strip().split()[1:], dtype=np.float64)
            gt_array = np.zeros((nhap, pos.shape[0]), dtype=np.uint8)
            cix = 0
            while cix < nhap:
                line = next(msmsfile)
                line = line.decode('utf-8')
                line = list(line.strip())
                try:
                    gt_array[cix, :] = np.array(line, dtype=np.uint8)
                except IndexError:
                    break
                cix += 1
            gtdict[str(rep)] = gt_array
            posdict[str(rep)] = pos
        elif line.startswith("OriginCount"):

            origcount[rep] = line.strip().split(":")[1]
        else:
            pass
    return(gtdict, posdict, origcount, freqtrace)


def fig1a_stats(freqtrace, time, Ne, gens):
    """
    freq of resistant allele at different time points in the past with
    varying amounts of selection
    """
    gens_past = [((gens * i) / 4 * Ne) for i in time]
    freq = []
    for rep in freqtrace.keys():
        time_trace = freqtrace[rep][:, 0]
        freq.append([freqtrace[rep][bisect(time_trace, g), 2]
                    for g in gens_past])
    freq_mean = [np.mean(ftime) for ftime in zip(*freq)]
    return(freq_mean)


def fig1a(msms, Ne, pops, reps, s, rho, theta, sp, smu, sAa, saa, sft, sff,
          gens):
    """
    Recreates data for Figure 1A
    """
    selp = np.arange(0, 0.1, 0.001)  # selection coefficient
    time = [20, 40, 60, 80]  # time in years
    freq = []
    inds = pops[0]
    nhap = inds * 2
    sAa = 0  # 1 if codominant
    for selco in selp:
        ms_params = {
                    'msms': msms,
                    'Ne': Ne,
                    'nhaps': nhap,
                    'nreps': reps,
                    'seg': s,
                    'rho': rho,
                    'selpos': sp,
                    'smu': smu,
                    'sAA': selco * 2 * Ne,
                    'sAa': selco * Ne * sAa,
                    'saa': saa,
                    'sft': sft,
                    'sff': sff}
        msms_base = ("{msms} -N {Ne} -ms {nhaps} {nreps} -s {seg} -r {rho} "
                     "-Sp {selpos} -Smu {smu} -SAA {sAA} -SAa {sAa} -Saa {saa}"
                     " -SF {sft} {sff} -oOC -Smark -oTrace")
        mscmd = msms_base.format(**ms_params)
        print(mscmd)
        msout = subprocess.Popen(mscmd, shell=True, stdout=subprocess.PIPE)
        gtdict, posdict, origcount, freqtrace = parse_msfile(msout, nhap, reps)
        # calc stats
        freq.extend(fig1a_stats(freqtrace, time, Ne, gens))
    dfFig1a = pd.DataFrame({'SAA': np.repeat(selp * 2 * Ne, len(time)),
                            'SAa': np.repeat(selp * Ne * sAa, len(time)),
                            'time': time * len(selco),
                            'freq': freq
                            })
    return(dfFig1a)


def fig1b_stats(gtdict, posdict, demesizelist, sp):
    """calculates the haplotype diversite of haplotypes carrying the resistant
       allele. Also: number of origins, total haplotype diversity, resistant
       haplotype congfig
    """
    # TO DO: calc below for all pops
    nR = []  # number of resistant haps
    nRmax = []  # max freq of resistant hap
    hd = []  # haplotype diversity; Depaulis and Veuille
    hapconfig = []  # haplotype configuration
    rfreq = []  # frequency of resistant allele
    popiix = []
    ix = 0
    for pop in demesizelist:
        popiix.append(range(ix, ix + pop))
        ix += pop
    for rep in gtdict.keys():
        smark = posdict[rep](sp)
        piix = gtdict[popiix[0]]
        riix = np.where(piix[:, smark] > 0)
        # riix = np.where(gtdict[:, smark] > 0)
        hapr = gtdict[rep][riix]
        uniqhaps = np.array([np.array(x) for x in set(tuple(x) for x in hapr)])
        hapfreq = np.array([len(hapr[np.all(hapr == x, axis=1)])
                            for x in uniqhaps], dtype=int)
        # full hap config
        n = sum(hapfreq)
        C_freq, C_count = np.unique(hapfreq, return_counts=True)
        C = np.zeros(n)
        C[C_freq - 1] = C_count

        # haplotype diversity
        Hd = 1 - sum([(((i+1)/float(n))**2) * c for i, c in enumerate(C)])
        M = max(C)  # greatest freq
        K = sum(C)  # number of haps
#        pdist = np.zeros((hapr.shape[0], hapr.shape[0]))
#        [np.count_nonzero(a!=b) for i, a in enumerate(hapr)
#         for j, b in enumerate(hapr) if j >= i]
        rfreq.append(riix.shape[0])
        hapconfig.append(C)
        nR.append(K)
        nRmax.append(M)
        hd.append(Hd)
#    hapsummR = [np.mean(pop, axis=0) for i, pop in enumerate(zip(*hapconfig))]
#    haparrayR = [np.vstack(i) for i in zip(*hapconfig)]
#    for conf in haparrayR:
#        uniq = np.array([np.array(x) for x in set(tuple(x) for x in conf)])
#        hapfreq = np.array([len(conf[np.all(conf == x, axis=1)])
#                            for x in uniqhaps], dtype=int)
#        hapconfig_commonR.append(zip(uniq, hapfreq))

    return(np.mean(nR), np.mean(hd), np.mean(nRmax), np.mean(rfreq))


def fig1b(msms, Ne, pops, reps, s, rho, theta, sp, smu, sAA, sAa, saa, sit):
    """Recreates data for Figure 1B
    """
    freqp = np.arange(0, 0.3, 0.001)  # freq start
    migp = np.arange(0, 0.3, 0.001)  # migration proportion
    inds = sum(pops)
    nhap = inds * 2
    demes = len(pops)
    for fq, m in zip(freqp, migp):
        orig = []
        ms_params = {
                    'msms': msms,
                    'Ne': Ne,
                    'nhaps': nhap,
                    'nreps': reps,
                    'seg': s,
                    'rho': rho,
                    'demes': demes,
                    'selpos': sp,
                    'smu': smu,
                    'sAA': sAA * 2 * Ne,
                    'sAa': sAa * Ne,  # 5000
                    'saa': saa,
                    'sit': sit,
                    'sif': fq,
                    'Nm': m * 4 * Ne}
        msms_base = ("{msms} -N {Ne} -ms {nhaps} {nreps} -s {seg} -r {rho} "
                     "-I {demes} 20 20 20 20 20 {Nm} "
                     "-Sp {selpos} -Smu {smu} -SAA {sAA} -SAa {sAa} -Saa {saa}"
                     " -SI {sit} {demes} {sif} {sif} {sif} {sif} {sif}"
                     " -oOC -Smark")
        mscmd = msms_base.format(**ms_params)
        print(mscmd)
        msout = subprocess.Popen(mscmd, shell=True, stdout=subprocess.PIPE)
        gtdict, posdict, origcount, freqtrace = parse_msfile(msout, nhap, reps,
                                                             sp)
        # calc stats
        nR, hd, nRmax, rfreq = fig1b_stats(gtdict, posdict, pops)
        orig.extend(sum(i > 1 for i in origcount) / float(reps))
    dfFig1b = pd.DataFrame({'freq': np.repeat(freqp, demes),
                            'mig': np.repeat(migp, demes),
                            'nRhap': nR,
                            'Rmaxf': nRmax,
                            'Rhapdiv': hd,
                            'Rfreq': rfreq,
                            'origins': np.repeat(orig, demes)
                            })
    return(dfFig1b)

if __name__ == '__main__':
    Ne = args.effectivesize
    pops = args.populations
    reps = args.reps
    msms = '/home/scott/programs_that_work/msms/bin/msms'
    s = 50  # number of seg sites
    rho = 15  # recombination rate, rho
    theta = 0  # theta, population mutation rate
    sp = 0.5  # position of the selected locus
    smu = 0.01  # mutation rate from wildtype to derived
    sAA = 0.01  # selection coeff homo; fig1B
    sAa = 0  # selection coeff het; fig1B
    saa = 0  # selection coeff wildtype homo; fig1B
    sit = 0.00024  # time of selection starting in 4 *Ne * gens; fig1B
    sft = 0  # time of selection stopping; fig1A
    sff = 0.30  # final freq of allele; fig1A
    gens = 12  # gens per year
    fig1a(msms, Ne, pops, reps, s, rho, theta, sp, smu, sAa, saa, sft, sff,
          gens)
    fig1b(msms, Ne, pops, reps, s, rho, theta, sp, smu, sAA, sAa, saa, sit)
