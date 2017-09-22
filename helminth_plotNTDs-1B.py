# -*- coding: utf-8 -*-
"""
@author: stsmall
"""
import sys
import pkg_resources
import numpy as np
import argparse
import subprocess
import pandas as pd
from libsequence.summstats import garudStats
if sys.version_info[0] > 2.9:
    raise Exception("Not compatible with python3")
pylib_v = pkg_resources.get_distribution("pylibseq").version
if pylib_v == '0.1.8':
    from libsequence.polytable import simData as simData
    from libsequence.summstats import polySIM as polySIM
else:
    from libsequence.polytable import SimData as simData
    from libsequence.summstats import PolySIM as polySIM

parser = argparse.ArgumentParser()
parser.add_argument('-N', '--effectivesize', type=int,
                    help='effective population size', required=True)
parser.add_argument('-p', '--populations', nargs='+', type=int, required=True,
                    help="population list")
parser.add_argument('-r', "--reps", type=int, required=True,
                    help="number of reps")
parser.add_argument('-rho', "--recombination", type=float, required=False,
                    default=0, help="population recombination rate, default 0")
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
        if line.startswith(b'//'):
            rep += 1
            trace = []
        elif line.startswith(b'Frequency'):
            line = next(msmsfile)
            while not line.strip():
                line = next(msmsfile)
            while not line.startswith(b'segsites'):
                trace.append(line.strip().split())
                line = next(msmsfile)
            trace = [map(float, x) for x in trace]
            freqtrace[str(rep)] = np.vstack(trace)
        elif line.startswith(b'positions'):
            pos = np.array(line.strip().split()[1:], dtype=np.float64)
            gt_array = np.zeros((nhap, pos.shape[0]), dtype=np.uint8)
            cix = 0
            while cix < nhap:
                line = next(msmsfile)
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
        elif line.startswith(b'OriginCount'):
            origcount[rep] = int(line.strip().split(b':')[1])
            # print("Orig:{}".format(line.strip().split(":")[1]))
        else:
            pass
    return(gtdict, posdict, origcount, freqtrace)


def hapbaxVmig_stats(gtdict, posdict, demesizelist, sp, origcount, sel, mig):
    """calculates the haplotype diversity of haplotypes carrying the resistant
       allele. Also: number of origins, total haplotype diversity, resistant
       haplotype congfig
    """
    pdist = []
    hapbax = []
    for rep in range(len(gtdict.keys())):
        rep = str(rep)
        smark = np.where(posdict[rep] == sp)[0]
        if len(smark) > 1:
            print("\nSkipping rep {}, smark gt 1\n".format(rep))
            continue
        piix = gtdict[rep][0:demesizelist[0]]  # index for the first pop
        riix = np.where(piix[:, smark] > 0)[0]  # location of the selected
        if riix.any():
            hapr = gtdict[rep][riix]
            uniqhaps_IBS = np.array([np.array(x) for x in set(tuple(x)
                                     for x in hapr)])
            hapfreq_IBS = np.array([len(hapr[np.all(hapr == x, axis=1)])
                                    for x in uniqhaps_IBS], dtype=int)
            uallel = hapr[:][:, smark]
            puallel = gtdict[rep][:, smark]
            hapr[:, smark] = 1  # sites that are IBS are read as uniqhaps
            # change all origins to 1 even if >1
            uniqhaps = np.array([np.array(x) for x in set(tuple(x)
                                 for x in hapr)])
            hapfreq = np.array([len(hapr[np.all(hapr == x, axis=1)])
                                for x in uniqhaps], dtype=int)
            print("\n#rep {}".format(rep))
            print("\nsel: {}\nmig: {}".format(sel, mig))
            print("#number of ORIGINS: {}".format(origcount[int(rep)]))
            print("#number of unique resistant ALLELES across ALL pops: {}".
                  format(len(np.unique(puallel[puallel > 0]))))
            print("#number of unique resistant ALLELES in SAMPLE pop: {}".
                  format(len(np.unique(uallel[uallel > 0]))))
            if uniqhaps_IBS.shape[0] > uniqhaps.shape[0]:
                print("#number of resistant HAPLOTYPES (IBS) in SAMPLE pop: "
                      "{}, frequencies {}".format(len(uniqhaps_IBS),
                                                  hapfreq_IBS))
                print("#number of hidden resistant HAPLOTYPES (IBS): {}".
                      format(uniqhaps_IBS.shape[0] - uniqhaps.shape[0]))
            print("#number of observable resistant HAPLOTYPES in SAMPLE pop: "
                  "{}, frequencies {}".format(len(uniqhaps), hapfreq))
            hapbax.append(len(uniqhaps))
            # full hap config
            n = sum(hapfreq)
            C_freq, C_count = np.unique(hapfreq, return_counts=True)
            C = np.zeros(piix.shape[0])
            C[C_freq - 1] = C_count
            # haplotype diversity
            Hd = 1 - sum([(((i+1)/float(n))**2) * c
                          for i, c in enumerate(C)])
            M = max(np.nonzero(C)[0]) + 1  # greatest non-zero position
            K = sum(C)  # number of haps
            # eveness from Chattopadhyay 2007
            lambda_e = sum([(float(hf)/sum(hapfreq))**2 for hf in hapfreq])
            Ds = 1.0/lambda_e
            Ev = Ds/uniqhaps.shape[0]
            rAfreq = riix.shape[0] / float(piix.shape[0])
            com = "#stats_r: Hd:{}\tKhaps:{}\tMaxAbsFreq:{}\tEv:{}\trFreq:{}\n"
            print(com.format(Hd, K, M, Ev, rAfreq))
            # fill pdist with below
            pdist.extend([np.count_nonzero(a != b) for i, a in enumerate(hapr)
                          for j, b in enumerate(hapr) if j > i])
            # popgen stats resistant
            gtpopr = [''.join(str(n) for n in y) for y in hapr]
            sdpopr = simData()
            sdpopr.assign_sep(posdict[rep], gtpopr)
            pspopr = polySIM(sdpopr)
            theta_r = pspopr.thetaw()
            tajd_r = pspopr.tajimasd()
            hprime_r = pspopr.hprime()
            pi_r = pspopr.thetapi()
            garudStats_r = garudStats(sdpopr)  # garud 2015
            # popgen stats all
            piix[:, smark] = 1
            gtpopa = [''.join(str(v) for v in y) for y in piix]
            sdpopa = simData()
            sdpopa.assign_sep(posdict[rep], gtpopa)
            pspopa = polySIM(sdpopa)
            theta_a = pspopa.thetaw()
            tajd_a = pspopa.tajimasd()
            hprime_a = pspopa.hprime()
            pi_a = pspopa.thetapi()
#            hapdiv_a = pspopa.hapdiv()
#            nhaps_a = pspopa.nhaps()
            garudStats_a = garudStats(sdpopa)  # garud 2015
            print("#popgen_r: theta_w:{}\ttheta_pi:{}\ttajD:{}\tfaywuH:{}".
                  format(theta_r, pi_r, tajd_r, hprime_r))
            print("#popgen_all: theta_w:{}\ttheta_pi:{}\ttajD:{}\tfaywuH:{}".
                  format(theta_a, pi_a, tajd_a, hprime_a))
            print("#garud2015_r: H12:{}\tH1:{}\tH2H1:{}".
                  format(garudStats_r['H12'], garudStats_r['H1'],
                         garudStats_r['H2H1']))
            print("#garud2015_all: H12:{}\tH1:{}\tH2H1:{}\n".
                  format(garudStats_a['H12'], garudStats_a['H1'],
                         garudStats_a['H2H1']))
        else:
            C = np.zeros(piix.shape[0])

    hapbaxm = np.mean(hapbax)
    hapbaxSE = np.std(hapbax) / len(hapbax)
    return(hapbaxm, hapbaxSE)


def hapbaxVmig_sims(msms, Ne, pops, reps, s, rho, theta, sp, smu, sAAc, sAac,
                    saa, sit, sif_l, selp, migp, threads, join_times, gens):
    """Simulation in msms, to explore the parameter space of migration and
    selection in determining the number of haplotype background carrying a
    resistance (derived) allele
    """
    nhap = sum(pops)
    demes = len(pops)
    hapbaxmean = []
    hapbaxSE = []

    for selco, saf in zip(selp, sif_l):
        sif = np.repeat(saf, demes)
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
            stats_out = hapbaxVmig_stats(gtdict, posdict, pops, sp, origcount,
                                         selco, m)
            hapbaxmean.append(stats_out[0])
            hapbaxSE.append(stats_out[1])
    # construct df for fig1b
    dfFig1b = pd.DataFrame({'mig': np.repeat(migp, len(selp)),
                            'sel': selp * len(migp),
                            'hapbaxmean': hapbaxmean,
                            'hapbaxSE': hapbaxSE,
                            })
    dfFig1b = dfFig1b.loc[:, ['sel', 'mig', 'hapbaxmean', 'hapbaxSE']]
    if sAAc == sAac:
        figname = "Fig1B_helminth-D-{}".format(rho)
    elif sAAc > sAac:
        figname = "Fig1B_helminth-A-{}".format(rho)
    elif sAac == 0:
        figname = "Fig1B_helminth-R-{}".format(rho)
    dfFig1b.to_csv(figname + ".csv")
    return(None)


if __name__ == '__main__':
    # universal options
    Ne = args.effectivesize
    pops = args.populations
    reps = args.reps
    msms = '/home/scott/programs_that_work/msms/bin/msms'
    threads = args.threads
    s = 34  # number of seg sites
    rho = args.recombination  # recombination rate, rho
    theta = 8.28  # UK 8.28, India 6.80, France 3.6, China 4.76, theta
    sp = 0.55555  # position of the selected locus
    smu = 0.01  # mutation rate from wildtype to derived, from theta/site
    gens = 12  # gens per year
    saa = 0  # selection coeff wildtype homo; fig1B

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

    # options for hapbaxVmig_sims
    migp = np.arange(0.00, 0.1, 0.01)  # migration proportion
    selpB = [0.0150, 0.005, 0.003]
    # selection based on estimates from freqVsel
    sif_l = [0, 0.01, 0.05]
    sit = 60  # when selection was first started
    join_times = 80  # when the farms were established
    # time in coalescent units : (gens * time) / 4*Ne
    # time in years : (coal_U * 4*Ne) / gen
    hapbaxVmig_sims(msms, Ne, pops, reps, s, rho, theta, sp, smu, sAAc, sAac,
                    saa, sit, sif_l, selpB, migp, threads, join_times, gens)
