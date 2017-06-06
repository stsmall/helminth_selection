#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon May 22 13:29:56 2017

@author: scott
"""
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('INfig', metavar="INfig", type=str,
                    help='path to vcf file')
args = parser.parse_args()


def renamefig(fig):
    """
    """
    f = open('{}.rename'.format(fig), 'w')
    with open(fig, 'r') as figb:
        for line in figb:
            if line.startswith(","):
                x = line.split(",")
                x[0] = "Index"
                f.write("{}\n".format(",".join(x)))
            else:
                if line.split(",")[3] == "R":
                    Rcount = 0
                    while line.strip().split(",")[3] == "R":
                        x = line.strip().split(",")
                        x[3] = "R{}".format(Rcount)
                        f.write("{}\n".format(",".join(x)))
                        line = figb.next()
                        Rcount += 1
                    if line.split(",")[3] == "S":
                        x = line.strip().split(",")
                        f.write("{}\n".format(",".join(x)))
    f.close()
    return(None)

if __name__ == '__main__':
    renamefig(args.INfig)
