from __future__ import print_function
import numpy as np 
import sys 
import math
import rpy2.robjects as robjects

def readGramsWithValues(filename):
    d = dict()
    lines = 0
    with open(filename) as f:
        for line in f:
            row = line.strip().split(",")
            qgram = row[0]
            phred = int(row[1])
            occur = int(row[3])
            fm = int(row[4])
            fmm = int(row[5])
            fer = 1.0 * fmm / occur #row[6]
            newPhred = int(round(-10.0 * math.log(fer, 10), 0)) if fer != 0 else phred
            d[(qgram, phred)] = (occur, fm, fmm, fer, newPhred)
            lines += 1
#           print (qgram, values)
    print ("lines read from arg1:", lines)
    return d

def readGramsAndSwapPhred(filename, phredMap):
    d = dict()
    lines = 0
    with open(filename) as f:
        for line in f:
            row = line.strip().split(",")
            qgram = row[0]
            phred = int(row[1])
            occur = int(row[3])
            fm = int(row[4])
            fmm = int(row[5])
            fer = 1.0 * fmm / occur #row[6]
            newPhred = phredMap[(qgram, phred)][4] if (qgram, phred) in phredMap else phred            
            #print (fmm, occur, fer, newPhred)
            if (qgram, newPhred) in d:
                tmp = d[(qgram, newPhred)]
                d[(qgram, newPhred)] = map(lambda x,y:x+y, tmp, [occur, fm, fmm])
            else:
                d[(qgram, newPhred)] = [occur, fm, fmm]
            lines += 1
#           print (qgram, values)
    print ("lines read from arg2:", lines)
    return d 

def get_pbinom(t, n, p, logV, lowerV):
    f = robjects.r['pbinom']
    #print (t, lowerV)
    t = t - 1 if t > 0 and not lowerV else t
    #print (t, lowerV)
    return tuple(f(t, n, p, log=logV, lower=lowerV))[0]

if __name__ == '__main__':
    args = sys.argv
    filename_out = args[3]
    d = readGramsAndSwapPhred(args[2], readGramsWithValues(args[1]))
    rows = []
    for qgram, phred in d.keys():
        occur, fm, fmm = d[(qgram, phred)]
        phred_real = math.pow(10.0, -float(phred) / 10.0)
        pvalf_left = get_pbinom(fmm, occur, phred_real, True, True)
        pvalf_right = get_pbinom(fmm, occur, phred_real, True, False)
        pvalf_lr = min(pvalf_left, pvalf_right)
        rows.append((qgram, phred, phred_real, occur, fm, fmm, float(fmm)/occur if occur != 0 else "inf", pvalf_left, pvalf_right, pvalf_lr))

    with open(filename_out, "w") as f:
        for row in sorted(rows, key=lambda row: row[9]):
            print (reduce(lambda x,y: str(x) + "," + str(y), row), file=f)

