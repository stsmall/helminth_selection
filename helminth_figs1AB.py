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
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument('-N', '--effectivesize', type=int,
                    help='effective population size', required=True)
parser.add_argument('-p', '--populations', nargs='+', type=int, required=True,
                    help="population list")
parser.add_argument('-r', "--reps", type=int, required=True,
                    help="number of reps")
parser.add_argument('-a', "--popall", action="store_true",
                    help="dont split by pops")
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
        elif line.startswith("Frequency"):
            line = next(msmsfile)
            while not line.strip():
                line = next(msmsfile)
            while not line.startswith("segsites"):
                trace.append(line.strip().split())
                line = next(msmsfile)
            trace = [map(float, x) for x in trace]
            freqtrace[str(rep)] = np.vstack(trace)
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
                except ValueError:
                    #print(line)
                    sdig = [i for i, s in enumerate(line) if not s.isdigit()]
                    x = 1
                    for s in sdig:
                        try:
                            line[s] = int(line[s], 36)
                        except ValueError:
                            line[s] = 35 + x
                            x += 1
                    gt_array[cix, :] = np.array(line, dtype=np.uint8)
                cix += 1
            gtdict[str(rep)] = gt_array
            posdict[str(rep)] = pos
        elif line.startswith("OriginCount"):
            origcount[rep] = line.strip().split(":")[1]
            print("Orig:{}".format(line.strip().split(":")[1]))
        else:
            pass
    return(gtdict, posdict, origcount, freqtrace)


def fig1a_stats(freqtrace, time, Ne, gens):
    """
    freq of resistant allele at different time points in the past with
    varying amounts of selection
    """
    gens_past = [(gens * i) / (4.0 * Ne) for i in time]
    freq = []
    for rep in freqtrace.keys():
        time_trace = freqtrace[rep][:, 0]
        freq.append([freqtrace[rep][np.argmin(time_trace > g), 2]
                    for g in gens_past])
    freq_mean = [np.mean(ftime) for ftime in zip(*freq)]
    return(freq_mean)


def fig1a(msms, Ne, pops, reps, s, rho, theta, sp, smu, sAa, saa, sft, sff,
          gens, selp, time):
    """
    Recreates data for Figure 1A
    """
    freq = []
    nhap = pops[0]
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
                            'time': time * len(selp),
                            'freq': freq
                            })
    dfFig1a.to_csv("Fig1A_helminth.csv")
    return(dfFig1a)


def fig1b_stats(gtdict, posdict, demesizelist, sp, popall):
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
    try:
        if popall:
            popiix.append(range(0, sum(demesizelist)))
        else:
            for pop in demesizelist:
                popiix.append(range(ix, ix + pop))
                ix += pop
        for rep in gtdict.keys():
            smark = np.where(posdict[rep] == sp)[0]
            for pop in popiix:
                piix = gtdict[rep][pop]
                riix = np.where(piix[:, smark] > 0)[0]
                if any(riix):
                    hapr = gtdict[rep][riix]
                    uniqhaps = np.array([np.array(x) for x in set(tuple(x)
                                         for x in hapr)])
                    hapfreq = np.array([len(hapr[np.all(hapr == x, axis=1)])
                                        for x in uniqhaps], dtype=int)
                    # full hap config
                    n = sum(hapfreq)
                    C_freq, C_count = np.unique(hapfreq, return_counts=True)
                    C = np.zeros(n)
                    C[C_freq - 1] = C_count
                    # haplotype diversity
                    Hd = 1 - sum([(((i+1)/float(n))**2) * c
                                  for i, c in enumerate(C)])
                    M = max(np.nonzero(C)[0]) + 1  # greatest non-zero position
                    K = sum(C)  # number of haps
            #        pdist = np.zeros((hapr.shape[0], hapr.shape[0]))
            #        [np.count_nonzero(a!=b) for i, a in enumerate(hapr)
            #         for j, b in enumerate(hapr) if j >= i]
                else:
                    C = np.zeros(piix.shape[0])
                    K = 0
                    M = 0
                    Hd = 0
                rfreq.append(riix.shape[0])
                hapconfig.append(C)
                nR.append(K)
                nRmax.append(M)
                hd.append(Hd)
    except:
        import ipdb; ipdb.set_trace()
#    hapsummR = [np.mean(pop, axis=0) for i, pop in enumerate(zip(*hapconfig))]
#    haparrayR = [np.vstack(i) for i in zip(*hapconfig)]
#    for conf in haparrayR:
#        uniq = np.array([np.array(x) for x in set(tuple(x) for x in conf)])
#        hapfreq = np.array([len(conf[np.all(conf == x, axis=1)])
#                            for x in uniqhaps], dtype=int)
#        hapconfig_commonR.append(zip(uniq, hapfreq))

    return(np.mean(nR), np.mean(hd), np.mean(nRmax), np.mean(rfreq))


def fig1b(msms, Ne, pops, reps, s, rho, theta, sp, smu, sAA, sAa, saa, sit,
          popall, freqp, migp):
    """Recreates data for Figure 1B
    """
    nhap = sum(pops)
    demes = len(pops)
    orig = []
    for fq in freqp:
        for m in migp:
            ms_params = {
                        'msms': msms,
                        'Ne': Ne,
                        'nhaps': nhap,
                        'nreps': reps,
                        'seg': s,
                        'rho': rho,
                        'demes': "{} {}".format(demes, (str(pops[0]) + ' ')*demes),
                        'selpos': sp,
                        'smu': smu,
                        'sAA': sAA * 2 * Ne,
                        'sAa': sAa * Ne,  # 5000
                        'saa': saa,
                        'sit': sit,
                        'sif': "{} {}".format(demes, (str(fq) + ' ') * demes),
                        'Nm': m * 4 * Ne}
            msms_base = ("{msms} -N {Ne} -ms {nhaps} {nreps} -s {seg} -r {rho} "
                         "-I {demes} {Nm} "
                         "-Sp {selpos} -Smu {smu} -SAA {sAA} -SAa {sAa} -Saa {saa}"
                         " -SI {sit} {sif}"
                         " -oOC -Smark")
            mscmd = msms_base.format(**ms_params)
            print(mscmd)
            msout = subprocess.Popen(mscmd, shell=True, stdout=subprocess.PIPE)
            gtdict, posdict, origcount, freqtrace = parse_msfile(msout, nhap, reps)
            # calc stats
            nR, hd, nRmax, rfreq = fig1b_stats(gtdict, posdict, pops, sp, popall)
            orig.append(sum([i > 1 for i in origcount]) / float(reps))
    if popall:
        arraylen = len(list(migp) * len(freqp))
        dfFig1b = pd.DataFrame({'freq': np.repeat(freqp, len(migp)),
                                'mig': list(migp) * len(freqp),
                                'nRhap': nR,
                                'Rmaxf': nRmax,
                                'Rhapdiv': hd,
                                'Rfreq': rfreq,
                                'origins': np.repeat(orig, arraylen)
                                })
        dfFig1b.to_csv("Fig1B_helminth.csv")
    else:
        dfFig1b = pd.DataFrame({'freq': np.repeat(freqp, demes),
                                'mig': np.repeat(migp, demes),
                                'nRhap': nR,
                                'Rmaxf': nRmax,
                                'Rhapdiv': hd,
                                'Rfreq': rfreq,
                                'origins': np.repeat(orig, arraylen*demes)
                                })
        dfFig1b.to_csv("Fig1B_helminth.csv")
    return(dfFig1b)

if __name__ == '__main__':
    Ne = args.effectivesize
    pops = args.populations
    reps = args.reps
    msms = '/home/scott/programs_that_work/msms/bin/msms'
    popall = args.popall
    # Fig1A
    selp = np.arange(0, 0.1, 0.001)  # selection coefficient
    time = [20, 40, 60, 80]  # time in years
    sft = 0  # time of selection stopping; fig1A
    sff = 0.30  # final freq of allele; fig1A
    # Fig1B
    freqp = np.arange(0, 0.2, 0.01)  # freq start
    migp = np.arange(0.000001, 0.01, 0.001)  # migration proportion
    sit = 0.00024  # time of selection starting in 4 *Ne * gens; fig1B
    sAA = 0.01  # selection coeff homo; fig1B
    sAa = 0  # selection coeff het; fig1B
    saa = 0  # selection coeff wildtype homo; fig1B
    # general
    s = 50  # number of seg sites
    rho = 15  # recombination rate, rho
    theta = 0  # theta, population mutation rate
    sp = 0.5  # position of the selected locus
    smu = 0.01  # mutation rate from wildtype to derived
    gens = 12  # gens per year
    # fig1a(msms, Ne, pops, reps, s, rho, theta, sp, smu, sAa, saa, sft, sff,
    #       gens, selp, time)
    fig1b(msms, Ne, pops, reps, s, rho, theta, sp, smu, sAA, sAa, saa, sit,
          popall, freqp, migp)
