#!/usr/bin/env python
# -*- coding: utf-8 -*-
import math, sys, pysam

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
            if occur > 10000:
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

def float2symbolic(phred, fer):
	newPhred = int(round(-10.0 * math.log(fer, 10), 0)) if fer != 0 else phred
	if newPhred > phred:
           return newPhred
       else:
           return phred

def replace(s, ch, index):
    tmp = list(s)
    tmp[index] = ch
    return ''.join(tmp)

def process(bampath, filename, phredMap):
	infile = pysam.Samfile(bampath, "rb")
	outfile = pysam.Samfile(filename, "wb", template=infile)
	for read in infile.fetch():
	    if read.is_unmapped:
                outfile.write(read)
		continue
	    if read.is_reverse:
		for i in xrange(len(read.seq)):
                    qgram = read.seq[i : i + NN]
                    if len(qgram) != NN or qgram.count('N') > 0:
		         continue
  	   	    qgram = reverse_complement(qgram[:])
		    first_p_phred = ord(read.qual[i]) - 33
		    new_first_p_phred = chr(float2symbolic(first_p_phred, phredMap[(qgram, first_p_phred)][3]) + 33) if (qgram, first_p_phred) in phredMap else read.qual[i]
	     	    read.qual = replace(read.qual, new_first_p_phred, i)
	    else:
	        for i in xrange(len(read.seq)):
                    qgram = read.seq[i - NN + 1: i + 1]
                    if len(qgram) != NN or qgram.count('N') > 0:
		        continue
	            last_p_phred = ord(read.qual[i]) - 33
		    new_last_p_phred = chr(float2symbolic(last_p_phred, phredMap[(qgram, last_p_phred)][3]) + 33) if (qgram, last_p_phred) in phredMap else read.qual[i]
		    read.qual = replace(read.qual, new_last_p_phred, i)
	    outfile.write(read)

NN = 4
		
if __name__ == '__main__':
    args = sys.argv
    d = readGramsWithValues(args[3])
    process(args[1], args[2], d)

		
		
		
		
		
		
		
		
		
		
