#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon May 22 13:29:56 2017

@author: scott
"""

f = open('figbrename', 'w')

with open('Fig1B_helminth-R.csv', 'r') as figb:
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
