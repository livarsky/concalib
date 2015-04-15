#!/usr/bin/env python
import numpy as np 
import sys 
import math, HTSeq
from collections import defaultdict
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

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
            d[(qgram, phred)] = (occur, fm, fmm, fer)
            lines += 1
#           print (qgram, values)
    print ("lines read from arg1:", lines)
    return d

bases = {0:'A', 1:'C', 2:'G', 3:'T'}
def reverse_complement(s, rev=True):
    """Return the reverse complement of a DNA sequence s"""
    _complement = dict(A="T", T="A", C="G", G="C", N="N")
    t = reversed(s) if rev else s
    try:
        rc = (_complement[x] for x in t)  # lazy generator expression
    except KeyError:
        return s
        
    return "".join(rc)  # join the elements into a string

def bestLike(likes):
    sorted_likes = sorted(likes)
    index_max1 = likes.index(sorted_likes[-1])
    index_max2 = likes.index(sorted_likes[-2])
    #print sorted_likes
    ratio = math.log(sorted_likes[-1] / sorted_likes[-2])
    return [ratio, bases[index_max1], bases[index_max2]]

def calcLikeHood(likehoods, qgram, phred_real):
    #A = 0, C = 1, G = 2, T = 3           
    #print qgram, phred_real, likehoods
    likehoods[0] = likehoods[0] + math.log(1.0 - phred_real, 10) if qgram[-1] == 'A' else likehoods[0] + math.log(phred_real / 3.0, 10)
    likehoods[1] = likehoods[1] + math.log(1.0 - phred_real, 10) if qgram[-1] == 'C' else likehoods[1] + math.log(phred_real / 3.0, 10)
    likehoods[2] = likehoods[2] + math.log(1.0 - phred_real, 10) if qgram[-1] == 'G' else likehoods[2] + math.log(phred_real / 3.0, 10)
    likehoods[3] = likehoods[3] + math.log(1.0 - phred_real, 10) if qgram[-1] == 'T' else likehoods[3] + math.log(phred_real / 3.0, 10)
    return likehoods

def get_genome(ref_path):
    """Return genome for random access"""
    seq = None
    i = 0
    for s in HTSeq.FastaReader(ref_path):
        seq = str(s)
        break
    return seq

if __name__ == '__main__':
    args = sys.argv
    d = readGramsWithValues(args[2])
    likehoods = [0, 0, 0, 0]
    fer_likes = [0, 0, 0, 0]
    pos = -1
    genome = {}
    genome_str = get_genome(args[1])
    fail_x, fail_y, good_x, good_y = [], [], [], []
    with open(args[3]) as f:
        for line in f:
            if len(line.split(":")) == 2:
                #print pos
                if pos >= 0:
                    #print pos
                    genome[pos] = (bestLike(likehoods), bestLike(fer_likes))
                    #print likehoods, fer_likes
                    if genome_str[pos] != genome[pos][0][1] or genome_str[pos] != genome[pos][1][1]:
                        print pos, genome[pos], genome_str[pos]
                        fail_x += [genome[pos][0][0]]
                        fail_y += [genome[pos][1][0]]
                    else:
                        #print pos, genome[pos], genome_str[pos]
                        good_x += [genome[pos][0][0]]
                        good_y += [genome[pos][1][0]]
                        
                pos = int(line.split(":")[1])
                likehoods = [0, 0, 0, 0]
                fer_likes = [0, 0, 0, 0]
            if len(line.split(",")) == 4:
                qgram, cigar, phred, strand = tuple(line.split(","))
                phred = int(phred)
                phred_real = math.pow(10.0, -float(phred) / 10.0)
                fer = d[(qgram, phred)][3]
                if fer == 0.0:
                    fer = 1e-10
                    #continue
                if fer == 1.0:
                    fer = 1.0 - 1e-10
                    #continue 
                strand = strand.strip()
                if strand == "r":
                    qgram = reverse_complement(qgram, False)
                likehoods = calcLikeHood(likehoods, qgram, phred_real)
		fer_likes = calcLikeHood(fer_likes, qgram, fer)

    #plt.plot(fail_x, fail_y, 'bo', label="fail", markersize=2)
    plt.plot(good_x, good_y, 'ko', label="good", markersize=2)
    plt.grid(True)
    plt.legend(loc=2)    
    plt.xlabel('likehood arg1')
    plt.ylabel('likehood arg2')
    #average_phred /= all_occur
    plt.savefig("graph.png")
    #plt.show()
    
        
