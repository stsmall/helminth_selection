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
# import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument('-N', '--effectivesize', type=int,
                    help='effective population size', required=True)
parser.add_argument('-p', '--populations', nargs='+', type=int, required=True,
                    help="population list")
parser.add_argument('-r', "--reps", type=int, required=True,
                    help="number of reps")
parser.add_argument('-fa', "--figA", action="store_true",
                    help="runfigA")
parser.add_argument('-fb', "--figB", action="store_true",
                    help="runfigB")
parser.add_argument('-D', "--dominant", action="store_true",
                    help="effect is dominant")
parser.add_argument('-A', "--additive", action="store_true",
                    help="effect is additive")
parser.add_argument('-t', "--threads", type=int, required=False,
                    default=1, help="threads")
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
            # print("Orig:{}".format(line.strip().split(":")[1]))
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


def fig1a(msms, Ne, pops, reps, s, rho, theta, sp, smu, sAAc, sAac, saac, sft,
          sff, gens, selp, time, threads):
    """
    Recreates data for Figure 1A
    """
    freq = []
    nhap = pops[0]
    orig = []
    for selco in selp:
        ms_params = {
                    'msms': msms,
                    'Ne': Ne,
                    'nhaps': nhap,
                    'nreps': reps,
                    'seg': s,
                    'theta': theta,
                    'rho': rho,
                    'selpos': sp,
                    'smu': smu,
                    'sAA': selco * sAAc * Ne,
                    'sAa': selco * sAac * Ne,
                    'saa': saa,
                    'sft': sft,
                    'sff': sff,
                    't': threads}
        msms_base = ("{msms} -N {Ne} -ms {nhaps} {nreps} -s {seg} -r {rho} "
                     "-Sp {selpos} -Smu {smu} -SAA {sAA} -SAa {sAa} -Saa {saa}"
                     " -SF {sft} {sff} -oOC -Smark -oTrace -threads {t}")
        mscmd = msms_base.format(**ms_params)
        print(mscmd)
        msout = subprocess.Popen(mscmd, shell=True, stdout=subprocess.PIPE)
        gtdict, posdict, origcount, freqtrace = parse_msfile(msout, nhap, reps)
        # calc stats
        freq.extend(fig1a_stats(freqtrace, time, Ne, gens))
        orig.append(np.repeat(sum([i > 1 for i in origcount])/float(reps),
                              len(time)))
    dfFig1a = pd.DataFrame({'SAA': np.repeat(selp * sAAc, len(time)),
                            'SAa': np.repeat(selp * sAac, len(time)),
                            'time': time * len(selp),
                            'freq': freq,
                            'orig': np.concatenate(orig).ravel()
                            })
    dfFig1a = dfFig1a.loc[:, ['SAA', 'SAa', 'time', 'freq', 'orig']]
    dfFig1a.to_csv("Fig1A_helminth.csv")
    return(dfFig1a)


def fig1b_stats(gtdict, posdict, demesizelist, sp):
    """calculates the haplotype diversite of haplotypes carrying the resistant
       allele. Also: number of origins, total haplotype diversity, resistant
       haplotype congfig
    """
    # TO DO: calc below for all pops
#    nR = []  # number of resistant haps
#    nRmax = []  # max freq of resistant hap
#    hd = []  # haplotype diversity; Depaulis and Veuille
#    hapconfig = []  # haplotype configuration
    rfreq = []  # frequency of resistant allele
#    ev = []
#    pdist = []
#    hapconfig_commonR = []
    Rplot = []
    Splot = []
    for rep in gtdict.keys():
        smark = np.where(posdict[rep] == sp)[0]
        piix = gtdict[rep][0:demesizelist[0]]
        riix = np.where(piix[:, smark] > 0)[0]
        if riix.any():
            hapr = gtdict[rep][riix]
            uniqhaps = np.array([np.array(x) for x in set(tuple(x)
                                 for x in hapr)])
            hapfreq = np.array([len(hapr[np.all(hapr == x, axis=1)])
                                for x in uniqhaps], dtype=int)
#                # full hap config
#                n = sum(hapfreq)
#                C_freq, C_count = np.unique(hapfreq, return_counts=True)
#                C = np.zeros(n)
#                C[C_freq - 1] = C_count
#                # haplotype diversity
#                Hd = 1 - sum([(((i+1)/float(n))**2) * c
#                              for i, c in enumerate(C)])
#                M = max(np.nonzero(C)[0]) + 1  # greatest non-zero position
#                K = sum(C)  # number of haps
#                # eveness
#                Ds = 1.0/sum([(float(hf)/(hapfreq))**2 for hf in hapfreq])
#                Ev = Ds/uniqhaps.shape[0]
            # pdist = np.zeros((hapr.shape[0], hapr.shape[0]))
            # fill pdist with below
            # [np.count_nonzero(a != b) for i, a in enumerate(hapr)
            # for j, b in enumerate(hapr) if j > i]
        else:
            hapfreq = np.array([0])
#                C = np.zeros(piix.shape[0])
#                K = 0
#                M = 0
#                Hd = 0
#                Ev = 0
#                pdist = 0
        Rplot.append(np.sort(hapfreq))
        Splot.append(piix.shape[0] - riix.shape[0])
        rfreq.append(riix.shape[0] / float(piix.shape[0]))
#            hapconfig.append(C)
#            print(C)
#            nR.append(K)
#            nRmax.append(M)
#            hd.append(Hd)
#            ev.append(Ev)
#            pdist.append()
    Rplot_m = np.zeros(max([len(i) for i in Rplot]))
    for r in Rplot:
        for i, hap in enumerate(r):
            Rplot_m[i] += hap
    Piplot = np.append(Rplot_m, sum(Splot))
    Rfreq = np.repeat(np.mean(rfreq), len(Piplot))
#    hapsummR = [np.mean(pop, axis=0) for i, pop in enumerate(zip(*hapconfig))]
#    haparrayR = [np.vstack(i) for i in zip(*hapconfig)]
#    for conf in haparrayR:
#        uniq = np.array([np.array(x) for x in set(tuple(x) for x in conf)])
#        hapf = np.array([len(conf[np.all(conf == x, axis=1)])
#                         for x in uniq], dtype=int)
#    hapconfig_commonR.append(zip(uniq, hapf))
#    np.mean(nR)
#    np.mean(hd)
#    np.mean(nRmax)
#    np.mean(rfreq)
    return(Piplot, Rfreq)


def fig1b(msms, Ne, pops, reps, s, rho, theta, sp, smu, sAAc, sAac, saa, sit,
          sif_l, selp, migp, threads, join_times, gens):
    """Recreates data for Figure 1B
    """
    nhap = sum(pops)
    demes = len(pops)
    orig = []
    origsd = []
    origmax = []
    piplot = []
    selpdf = []
    mdf = []
    strdf = []
    freqr = []
    sif = sif_l * demes
    for selco in selp:
        for m in migp:
            ms_params = {
                        'msms': msms,
                        'Ne': Ne,
                        'nhaps': nhap,
                        'nreps': reps,
                        'seg': s,
                        'theta': theta,
                        'rho': rho,
                        'demes': "{} {}".format(demes,
                                                " ".join(map(str, pops))),
                        'selpos': sp,
                        'smu': smu,
                        'sAA': selco * sAAc * Ne,
                        'sAa': selco * sAac * Ne,
                        'saa': saa,
                        'sit': (gens * sit) / (4.0 * Ne),
                        'sif': "{} {}".format(demes, " ".join(map(str, sif))),
                        'Nm': m * 4 * Ne,
                        't': threads,
                        'time': (gens * join_times) / (4.0 * Ne)}
            msms_base = ("{msms} -N {Ne} -ms {nhaps} {nreps} -s {seg} "
                         "-r {rho} -I {demes} {Nm} -Sp {selpos} -Smu {smu} "
                         "-SAA {sAA} -SAa {sAa} -Saa {saa}"
                         " -SI {sit} {sif} -ej {time} 2 1 -ej {time} 3 1 "
                         "-ej {time} 4 1 -ej {time} 5 1 "
                         " -oOC -Smark -SFC -threads {t}")
            mscmd = msms_base.format(**ms_params)
            print(mscmd)
            msout = subprocess.Popen(mscmd, shell=True, stdout=subprocess.PIPE)
            gtdict, posdict, origcount, freqtrace = parse_msfile(msout, nhap,
                                                                 reps)
            # calc stats
            Piplot, freqR = fig1b_stats(gtdict, posdict, pops, sp)
            piplot.extend(Piplot)
            freqr.extend(freqR)
#            orig.append(np.repeat(sum([i > 1 for i in origcount])/float(reps),
#                        len(time)))
            orig.append(np.repeat(np.mean(origcount), len(Piplot)))
            origsd.append(np.repeat(np.std(origcount), len(Piplot)))
            origmax.append(np.repeat(np.max(origcount),len(Piplot)))
            selpdf.extend([selco]*len(Piplot))
            mdf.extend([m]*len(Piplot))
            strtemp = ["R"] * (len(Piplot) - 1)
            strtemp.extend("S")
            strdf.extend(strtemp)
    dfFig1b = pd.DataFrame({'sel': selpdf,
                            'mig': mdf,
                            'idHap': strdf,
                            'piHap': piplot,
                            'freqR': freqr,
                            'orig': np.concatenate(orig).ravel(),
                            'origsd': np.concatenate(origsd).ravel(),
                            'origmax': np.concatenate(origmax).ravel()
                            })
    dfFig1b = dfFig1b.loc[:, ['sel', 'mig', 'idHap', 'piHap', 'freqR', 'orig',
                              'origsd', 'origmax']]
    dfFig1b.to_csv("Fig1B_helminth.csv")
    return(None)

if __name__ == '__main__':
    # universal
    Ne = args.effectivesize
    pops = args.populations
    reps = args.reps
    msms = '/home/ssmall2/programs_that_work/msms/bin/msms'
    threads = args.threads
    s = 34  # number of seg sites
    rho = 8.0  # recombination rate, rho
    theta = 8.28  # UK 8.28, India 6.80, France 3.6, China 4.76, theta
    sp = 0.5  # position of the selected locus
    smu = 0.01  # mutation rate from wildtype to derived
    gens = 12  # gens per year
    saa = 0  # selection coeff wildtype homo; fig1B
    # selection coefficient
    selpA = np.arange(0.00, 0.1, 0.001)  # selection coefficient

    # dominance
    if args.dominant:
        # dominant
        sAAc = 2
        sAac = 2
    elif args.additive:
        # recessive
        sAAc = 2
        sAac = 1
    else:
        sAAc = 2
        sAac = 0

    # Fig1A
    time = [20, 40, 60, 80]  # time in years
    sft = 0  # time of selection stopping; fig1A
    sff = 0.33  # final freq of allele; fig1A

    # Fig1B
    migp = np.arange(0.00, 0.01, 0.001)  # migration proportion
    selpB = np.arange(0.00, 0.06, 0.006)  # selection based on Fig1A
    sif_l = [0]
    sit = 60  # when selection was first started
    join_times = 80  # when the farms were established
    # time in coalescent units : (gens * time) / 4*Ne
    # time in years : (coal_U * 4*Ne) / gen

    if args.figA:
        fig1a(msms, Ne, pops, reps, s, rho, theta, sp, smu, sAAc, sAac, saa,
              sft, sff, gens, selpA, time, threads)
    if args.figB:
        fig1b(msms, Ne, pops, reps, s, rho, theta, sp, smu, sAAc, sAac, saa,
              sit, sif_l, selpB, migp, threads, join_times, gens)
