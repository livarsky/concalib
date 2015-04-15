#!/usr/bin/env python
import numpy as np 
import sys 
import math
from collections import defaultdict
import matplotlib.pyplot as plt

def readGramsWithValues(filename):
    d = dict()
    lines = 0
    with open(filename) as f:
        for line in f:
            row = line.strip().split(",")
            qgram = row[0]
            phred = row[1]
            occur = int(row[3])
            fm = int(row[4])
            fmm = int(row[5])
            fer = 1.0 * fmm / occur #row[6]
            d[(qgram, phred)] = (occur, fm, fmm, fer)
            lines += 1
#           print (qgram, values)
    print ("lines read from arg1:", lines)
    return d

if __name__ == '__main__':
    args = sys.argv
    d = readGramsWithValues(args[1])
    rows = []
    metric_log = 0.0
    metric = 0.0
    metric_real = 0.0
    metric_const = 0.0
    metric_random = 0.0
    all_occur = 0
    average_phred = 0.0
    average_fer = 0.0
    
    phred2fer = defaultdict(list)
    for qgram, phred in d.keys():
        occur, fm, fmm, fer = d[(qgram, phred)]
        if fmm < 5:
            continue
        all_occur += occur
        phred_real = math.pow(10.0, -float(phred) / 10.0)
        #average_phred += phred_real * occur
        phred2fer[phred].append(-10.0 * math.log(fer, 10) if fer !=0 else 0)
        average_fer += fmm
    
    #plt.plot(sorted(phred2fer.keys()), [sum(phred2fer[x]) / len(phred2fer[x]) for x in sorted(phred2fer.keys())], 'bo', label="avg. fer", markersize=2)
    #plt.plot(sorted(phred2fer.keys()), [float(x) for x in sorted(phred2fer.keys())], 'r', label="phred", markersize=1)
    
    #plt.grid(True)
    #plt.legend(loc=2)
    #plt.xlabel('phred')
    #plt.ylabel('-avg 10 log fer')
    #average_phred /= all_occur
    #print average_fer, all_occur
    average_fer /= all_occur
    const = float(args[2]) if len(args) > 2 else average_fer
    all_fmm = 0
    all_fm = 0
    metric1 = 0
    metric2 = 0
    def my_random():
        mu, sigma = 0.0, 0.02 # mean and standard deviation
        return np.random.normal(mu, sigma)
    
    for qgram, phred in d.keys():
        occur, fm, fmm, fer = d[(qgram, phred)]
        if occur < 100 or fmm < 2:
            continue
        phred_real = math.pow(10.0, -float(phred) / 10.0)
        metric += fm * math.pow(phred_real, 2) + fmm * math.pow(1.0 - phred_real, 2)
    #all_fmm += fmm
    #all_fm += fm        
        metric_real += fm * math.pow(fer - 0.0, 2) + fmm * math.pow(1.0 - fer, 2)
        metric_const += fm * math.pow(const - 0.0, 2) + fmm * math.pow(1.0 - const, 2)
        myrand = math.fabs(my_random())
        metric_random += fm * math.pow(myrand - 0.0, 2) + fmm * math.pow(1.0 - myrand, 2)
        
        all_occur += 1
        if fer == 0.0:
            fer = 1e-150
        #print fer, math.log(fer, 10)
        
        metric_log += fm * math.pow(math.log(fer, 10) + 150, 2) + fmm * math.pow(0 - math.log(fer, 10), 2)
    
    
    #m = math.sqrt(metric1 / all_fm + metric2 / all_fmm)
    metric_const = math.sqrt(metric_const)
    metric_real = math.sqrt(metric_real)
    metric = math.sqrt(metric)
    metric_random = math.sqrt(metric_random)
    metric_log = math.sqrt(metric_log)
    
    print "average fer =", average_fer
    print "if error rate =", const, ", metric(const) =", metric_const
    print "if error rate = fer, metric(real) =", metric_real
    print "if error rate ~ phred, metric(phred) =", metric
    print "if error rate ~ random, metric(random) =", metric_random
    
    print "metric_phred > metric_real", metric > metric_real
    print "metric_const > metric_real", metric_const > metric_real
    
    print metric_log 
    #plt.show()
    
        
