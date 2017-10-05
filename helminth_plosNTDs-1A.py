# -*- coding: utf-8 -*-
"""
Created on Tue May  2 13:00:27 2017

@author: scott
"""
from __future__ import print_function
import pkg_resources
import numpy as np
import argparse
import subprocess
import pandas as pd
import random
from itertools import combinations
pylib_v = pkg_resources.get_distribution("pylibseq").version
if pylib_v == '0.1.8':
    from libsequence.polytable import simData as simData
    from libsequence.summstats import polySIM as polySIM
    from libsequence.fst import fst as fst
else:
    from libsequence.polytable import SimData as simData
    from libsequence.summstats import PolySIM as polySIM
    from libsequence.fst import Fst as fst
try:
    import configparser
except ImportError:
    import ConfigParser as configparser
parser = argparse.ArgumentParser()
parser.add_argument('config', metavar='config', type=str,
                    help='config file')
parser.add_argument('-o', "--outfile", type=str, required=True,
                    help='outfile name')
parser.add_argument('-m', "--msfile", default=None, type=str,
                    help='path to ms-formatted file')
parser.add_argument("--perm", action="store_true",
                    help="run permutation test")
args = parser.parse_args()


def nonblank_lines(f):
    """
    """
    for l in f:
        line = l.rstrip()
        if line:
            yield line


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
            if line.startswith('//'):
                line = next(msfile)
                if line.startswith('segsites'):
                    s = line.strip().split()[-1]
                    if int(s) > 0:
                        rep += 1
                    else:
                        continue
            elif line.startswith('positions'):
                pos = np.array(line.strip().split()[1:], dtype=np.float64)
                gt_array = np.zeros((nhap, pos.shape[0]), dtype=np.uint8)
                cix = 0
                while cix < nhap:
                    line = next(msfile)
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
        if line != b'':
            if line.startswith(b'//'):
                line = next(msfile)
                if line.startswith(b'segsites'):
                    s = line.strip().split()[-1]
                    if int(s) > 0:
                        rep += 1
                    else:
                        continue
            elif line.startswith(b'positions'):
                line = line.decode()
                pos = np.array(line.strip().split()[1:], dtype=np.float64)
                gt_array = np.zeros((nhap, pos.shape[0]), dtype=np.uint8)
                cix = 0
                while cix < nhap:
                    line = next(msfile)
                    line = line.decode()
                    line = list(line.strip())
                    gt_array[cix, :] = np.array(line, dtype=np.uint8)
                    cix += 1
                gtdict[str(rep)] = gt_array
                posdict[str(rep)] = pos
        else:
            break
    return(posdict, gtdict)


def calcfst(pops, posdict, gtdict, L, snp=True, hud=True):
    """
    """
    ix = 0
    pw = len(pops)
    fstarray = np.zeros([len(posdict.keys()), int((pw*(pw-1))/2)])
    # Observed FST
    pix = 0
    popdict = {}
    for p in pops:
        popdict[pix] = list(range(ix, ix + p))
        ix += p
        pix += 1
    for r in gtdict.keys():
        fst_obs = []
        if snp:
            pos_snp = [(np.random.choice(posdict[r]))]
            gt_snp = np.where(posdict[r] == pos_snp)[0]
        else:
            pos_snp = list(posdict[r])
        for x, y in combinations(popdict.keys(), 2):
            popX = gtdict[r][popdict[x]]
            popY = gtdict[r][popdict[y]]
            sdfst = simData()
            geno = np.vstack([popX, popY])
            geno_fst = geno[:, gt_snp]
            gtpop = [''.join(str(n) for n in y) for y in geno_fst]
            gtpop_fst = [i.encode() for i in gtpop]
            sdfst.assign_sep(pos_snp, gtpop_fst)
            size = [popX.shape[0], popY.shape[0]]
            f1 = fst(sdfst, size)
            if hud:
                fst_obs.append(f1.hsm())
            else:
                fst_obs.append(f1.slatkin())
                # fst_obs.append(f1.hbk())
#        # pi
        for z in popdict.keys():
            sdpop = simData()
            popZ = gtdict[r][popdict[z]]
            geno = popZ[:, gt_snp]
            gtpop = [''.join(str(n) for n in y) for y in geno]
            gtpop_pi = [i.encode() for i in gtpop]
            sdpop.assign_sep(pos_snp, gtpop_pi)
            pspopr = polySIM(sdpop)
            pi = pspopr.thetapi()
            print("pop {} pi:{}".format(z, pi/L))
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
#        fst_t.append(f1.slatkin())
        fst_t.append(f1.hsm())
#        fst_t.appen(f1.hbk())
    # mark significant FST
    fst_tnp = np.array(fst_t)
    Fstdist = [len(np.where(f > fst_tnp)[0]) for f in fstarray]
    return([1-(f/float(n_perm)) for f in Fstdist])


def runmssims(ms, Ne, migp, pops, reps, theta, rho, length, gens, time,
              outfile):
    """
    """
    nhap = sum(pops)
    demes = len(pops)
    fstpd = []
    migpd = []
    timepd = []
    for m in migp:
        for t in time:
            tgen = (gens * t) / (4.0 * Ne)
            ms_params = {
                        'ms': ms,
                        'nhaps': nhap,
                        'nreps': reps,
                        'theta': theta,
                        'rho': rho,
                        'L': length,
                        'demes': "{} {}".format(demes,
                                                " ".join(map(str, pops))),
                        'Nm': m * Ne * 4,
                        'join': "".join(["-ej {} {} {} ".format(tgen, i, i+1)
                                         for i in range(1, demes)])
                        }
            if rho > 0:
                msms_base = ("{ms} {nhaps} {nreps} -t {theta} "
                             "-r {rho} {L} -I {demes} {Nm} {join}")
            else:
                msms_base = ("{ms} {nhaps} {nreps} -t {theta} "
                             "-I {demes} {Nm} {join}")
            mscmd = msms_base.format(**ms_params)
            print(mscmd)
            msout = subprocess.Popen(mscmd, shell=True, stdout=subprocess.PIPE)
            # parse
            posdict, gtdict = readms2(msout, pops)
            # calc FST
            fstarray = calcfst(pops, posdict, gtdict, length)
            fstpd.append(np.nanmean(fstarray, axis=1))
            migpd.extend(np.repeat(m, fstarray.shape[0]))
            timepd.extend(np.repeat(t, fstarray.shape[0]))
    dfFig1a = pd.DataFrame({'mig': migpd,
                            'time': timepd,
                            'fst': np.concatenate(fstpd).ravel()
                            })
    dfFig1a = dfFig1a.loc[:, ['time', 'mig', 'fst']]
    dfFig1a.to_csv("{}.csv".format(outfile))
    return(None)


if __name__ == '__main__':
    config = configparser.ConfigParser()
    outfile = args.outfile
    config.read(args.config)
    sh = 'simulation'
    Ne = config.getint(sh, 'effectivesize')
    pops = list(map(int, config.get(sh, 'popsizes').split(",")))
    reps = config.getint(sh, "reps")
    ms = '/usr/bin/ms'  # default path
    gens = config.getint(sh, 'gens')
    mu = config.getfloat(sh, 'mu')
    L = config.getint(sh, 'length')
    theta = 4 * Ne * mu * L
    rhorat = config.getfloat(sh, 'rho')
    rho = rhorat * theta
    mig = list(map(float, config.get(sh, 'mig').split(",")))
    if len(mig) > 1:
        if mig[-1] < mig[1]:
            migp = np.arange(mig[0], mig[1], mig[2])
    else:
        migp = mig
    time = list(map(int, config.get(sh, 'join_times').split(",")))
    if args.msfile is not None:
        posdict, gtdict = readms(args.msfile, pops)
        fstarray = calcfst(pops, posdict, gtdict, L)
        a = np.round(np.nanmean(fstarray, axis=0), 2)
        print("Fst average value: {}".format(a[0]))
        if args.perm:
            fstpvalue = permtest(gtdict, posdict, pops, args.n_perm, a)
            print("[%s]" % ", ".join(map(str, fstpvalue)))
    else:
        runmssims(ms, Ne, migp, pops, reps, theta, rho, L, gens, time, outfile)
