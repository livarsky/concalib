#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
%prog <REF> <BAM> <INT q> <INT n> 

Discover and print (stdout) Context-Specific Sequencing Errors.

<REF>        Reference genome 
<BAM>        BAM file
<INT q>      Length of <q>-gram to search for
<INT n>      Maximal allowed number of Ns that the <q>-gram contains

For more details, see:

Manuel Allhoff, Alexander Schoenhuth, Marcel Martin, Ivan G. Costa, 
Sven Rahmann and Tobias Marschall.
Discovering motifs that induce sequencing errors. 
BMC Bioinformatics, 2013. 14(Suppl 5):S1, 
DOI: 10.1186/1471-2105-14-S5-S1, 
http://www.biomedcentral.com/1471-2105/14/S5/S1

Author: Manuel Allhoff
"""

from __future__ import print_function
from optparse import OptionParser
import math, sys, HTSeq, pysam, re
import scipy.misc as sc
import rpy2.robjects as robjects
import array
from collections import defaultdict
import itertools

class HelpfulOptionParser(OptionParser):
    """An OptionParser that prints full help on errors."""
    def error(self, msg):
        self.print_help(sys.stderr)
        self.exit(2, "\n%s: error: %s\n" % (self.get_prog_name(), msg))


def get_annotate_qgram(genome, genome_annotate, q):
    """Compute for each q-gram in the genome its (composed) strand bias table. 
    Consider therefore the q-gram as well as its reverse complement."""
    #consider separately pileups of q-grams last and first positions
    qgram_last = {}
    qgram_first = {}
    j = 0 #counter for status info
    k = 0
    l = 0
    #pass through entire genome to analyse each q-grams' pileup
    for i in range(len(genome) - q):
        j += 1
        if j % 20000000 == 0: 
            print('%s / %s positions considered for q-gram annotation' %(j, len(genome)), file=sys.stderr)
        
        qgram = genome [i : i+q]
        qgram = qgram.upper()
        

        if len(qgram) != qgram.count('A') + qgram.count('C') + qgram.count('G') + qgram.count('T'):
            #print("Warning: q-gram contains other letters than A,C,G and T, ignore q-gram" ,file=sys.stderr)
            #print(qgram, file=sys.stderr)
            k += 1
            continue
        else:
            l += 1

        qgram_effect_last = [genome_annotate[0][i + q], genome_annotate[1][i + q], genome_annotate[2][i + q], genome_annotate[3][i + q]]
        
        #q-grams on forward direction, analyse therefore their last positions
        qgram_last[qgram] = _add_listelements(qgram_last[qgram], qgram_effect_last) if qgram_last.has_key(qgram) else qgram_effect_last
        
        #q-gram on reverse strand, analyse therefore their first positions. Furthermore, switch read direction
        qgram_effect_first = [genome_annotate[1][i], genome_annotate[0][i], genome_annotate[3][i], genome_annotate[2][i]]
        qgram_first[qgram] = _add_listelements(qgram_first[qgram], qgram_effect_first) if qgram_first.has_key(qgram) else qgram_effect_first
        
    #combine the q-grams on the forward and reverse strand
    qgram_annotate = {}
    for qgram in qgram_last:
        qgram_annotate[qgram] = qgram_last[qgram]
        qgram_rev = reverse_complement(qgram)
        if qgram_rev in qgram_first:
            result = qgram_annotate[qgram][:]
            result = _add_listelements(result, qgram_first[qgram_rev])
            qgram_annotate[qgram] = result
    print("Warning: %s q-grams of %s contain other letters than A,C,G and T, ignore these q-grams" %(k, l) ,file=sys.stderr)
    return qgram_annotate


def reverse_complement(s, rev=True):
    """Return the reverse complement of a DNA sequence s"""
    _complement = dict(A="T", T="A", C="G", G="C", N="N")
    t = reversed(s) if rev else s
    try:
        rc = (_complement[x] for x in t)  # lazy generator expression
    except KeyError:
        return s
        
    return "".join(rc)  # join the elements into a string


def get_pbinom(t, n, p, logV, lowerV):
    f = robjects.r['pbinom']
    #print (t, lowerV)
    t = t - 1 if t > 0 and not lowerV else t
    #print (t, lowerV)
    return tuple(f(t, n, p, log=logV, lower=lowerV))[0]

def get_pvalue(forward_match, reverse_match, forward_mismatch, reverse_mismatch):
    """Return p-value of given Strand Bias Table"""
        #decide whether Fisher's exact Test or ChiSq-Test should be used
    if forward_match > 5000 or reverse_match > 5000 or forward_mismatch > 5000 or reverse_mismatch > 5000:
        f = robjects.r['chisq.test']
        test = 'chisq'
    else:
        f = robjects.r['fisher.test']
        test = 'fisher'
    
    matrix = [forward_match, reverse_match, forward_mismatch, reverse_mismatch]
    table = robjects.r.matrix(robjects.IntVector(matrix), nrow=2)
    p_value_tmp = f(table)[0] if test == 'fisher' else f(table)[2]
    p_value = tuple(p_value_tmp)[0] #some necessary magic for r object
    return p_value


def get_genome(ref_path, learn_chrom):
    """Return genome for random access"""
    seq = None
    i = 0
    for s in HTSeq.FastaReader(ref_path):
        if s.name != learn_chrom:
            i += 1
            continue
        else:
            seq = str(s)
            break

    if i == 1:
        learn_chrom = s.name
        seq = str(s)
    elif seq is None:
        parser.error("Sorry, the Chromosome that is using for training (%s) is not contained \
        in the reference genome (%s)! Please use -c option!" %(learn_chrom, ref_path))
    
    return seq, learn_chrom


def get_annotate_genome(genome, bampath, learn_chrom):
    """Return strand bias table for each genome position on the base of the alignment. 
    The table is represented as a list: 
    [forward match, reverse match, forward mismatch, reverse mismatch]."""
    cigar_codes = {'M':0, 'I':1, 'D':2, 'N':3, 'S':4, 'H':5, 'P':6}
    samfile = pysam.Samfile(bampath, "rb")

    #initialize
    fm = array.array('H')
    fmm = array.array('H')
    rm = array.array('H')
    rmm = array.array('H')
    dict_rm = defaultdict(int)
    dict_rmm = defaultdict(int)
    dict_fm = defaultdict(int)
    dict_fmm = defaultdict(int)
	
    for i in range(len(genome)):
        fm.append(0)
        fmm.append(0)
        rm.append(0)
        rmm.append(0)
    
    j = 0 #counter for status info
    #consider each read
    for read in samfile.fetch():
        if read.is_unmapped or samfile.getrname(read.rname) != learn_chrom:
            continue
    #print(samfile.getrname(read.rname), file=sys.stderr)
        j += 1
        ref_pos = read.pos
        if j % 1000000 == 0: 
            print('%s reads considered for genome annotation ' %j, file=sys.stderr)
        
        #analyse CIGAR string to consecutively compute strand bias tables
        if read.cigar is not None:
            bias, current_pos_ref, current_pos_read = 0, 0, 0 
            for code, length in read.cigar:
                if code is cigar_codes['S']:
                    current_pos_read += length
                elif code is cigar_codes['M']:
                    if read.is_reverse:
                        for i in range(NN - 1, length):
                            pos = ref_pos + bias + current_pos_ref + i
                            qgram = read.seq[i + current_pos_read + bias - NN + 1 : i + current_pos_read + bias + 1]
			    flag_continue = False
                            for j in range(0, NN - 1):
                                if read.seq[i + current_pos_read + bias - j - 1] == 'N':#!= genome[pos - j - 1]:
                                    flag_continue = True
                                    break
                            if flag_continue is True:
                                continue
                            first_p_phred = ord(read.qual[i + current_pos_read + bias - NN + 1]) - 33
#			    print ("reverse", qgram, last_p_phred)
			    if read.seq[i + current_pos_read + bias] == genome[pos]:
                                fm[pos] += 1
				dict_fm[(reverse_complement(qgram), first_p_phred)] += 1
                            else:
                                fmm[pos] += 1
				dict_fmm[(reverse_complement(qgram), first_p_phred)] += 1
                    else:
                        for i in range(NN - 1, length):
                            pos = ref_pos + bias + current_pos_ref + i
                            qgram = read.seq[i + current_pos_read + bias - NN + 1 : i + current_pos_read + bias + 1]
			    flag_continue = False
                            for j in range(0, NN - 1):
                                if read.seq[i + current_pos_read + bias - j - 1] == 'N':#!= genome[pos - j - 1]:
                                    flag_continue = True
                                    break
                            if flag_continue is True:
                                continue
                            last_p_phred = ord(read.qual[i + current_pos_read + bias]) - 33
#			    print ("forward", qgram, last_p_phred)
			    if read.seq[i + current_pos_read + bias] == genome[pos]:
                                fm[pos] += 1
				dict_fm[(qgram, last_p_phred)] += 1
                            else:
                                fmm[pos] += 1
				dict_fmm[(qgram, last_p_phred)] += 1
                    bias += length
                elif code is cigar_codes['I']: 
                    current_pos_read += length
                elif code is cigar_codes['D']: 
                    current_pos_ref += length
                else:
                    print(code, length, file=sys.stderr)
                    
    #return (fm, rm, fmm, rmm)
    return (dict_fm, dict_rm, dict_fmm, dict_rmm)


def add_n(qgram_annotate, n, q):
    """Extend qgram_annotate by adding q-grams which contain Ns."""
    to_add = {}
    
    if n == 0:
        print('No q-grams with Ns to add' , file = sys.stderr)
        return qgram_annotate

    i = 0 #counter for status info    
    #consider each possible q-gram for the given length q and number n
    for qgram_with_n in _get_all_qgrams(['A','C','G','T','N'], ['A','C','G','T','N'], q, 1, n):
        qgram_with_n = qgram_with_n[0] #q-gram that contain at least one and at most n Ns
        i += 1
        if i % 20000000 == 0: 
            print(' %s q-grams with N considered' %(i), file = sys.stderr)

        #compute all concrete q-grams of n_qgram
        possible_qgrams = get_qgramlist(qgram_with_n)
        sb_table = [0,0,0,0] #initialize strand bias table
        for p_qgram in possible_qgrams:
            if qgram_annotate.has_key(p_qgram):
                sb_table = _add_listelements(sb_table, qgram_annotate[p_qgram][:])
        
        #does n containing q-gram corresponds to a combosed strand bias table?
        if sb_table != [0,0,0,0]: 
            to_add[qgram_with_n] = sb_table

    #extend qgram_annotate
    qgram_annotate.update(to_add)
    
    return qgram_annotate


def _add_listelements(a, b):
    """Add i-th element of list a and list b"""
    for i in range(len(a)):
        a[i] += b[i]
    
    return a


def get_qgramlist(qgram_with_n):
    """Return list of all possible q-grams for q-grams which contain Ns"""
    if qgram_with_n.count("N")==0:
        return [qgram_with_n]
    return _get_qgramlist_help(qgram_with_n, [])


def _get_qgramlist_help(qgram, result):
    """Substitute N with A, C, G or T and return resulting q-grams"""
    if qgram.count('N') == 0:
        result.append(qgram)
    else:
        i = qgram.index("N")
        for n in ['A','C','G','T']:
            new_qgram = qgram[:i] + n + qgram[i+1:]
            _get_qgramlist_help(new_qgram, result)
        return result


def _get_all_qgrams(alphabet, erg, length, level, n):
    """Return all possible q-grams of the given length, over the given alphabet
    and with at least one and at most n Ns"""
    if length == level:
        yield erg
    else:
        for letter in alphabet:
            for el in erg:
                #not too many Ns
                if letter == 'N' and el.count('N') >= n:
                    continue
                #not too less Ns
                if length - level <= 1 and el.count('N') == 0 and letter != 'N':
                    continue

                for r in _get_all_qgrams(alphabet, [el + letter], length, level+1, n):
                    yield r


def get_sb_score(qgram_annotate):
    """Calculate Strand Bias score (based on p-value) for each q-gram"""
    results = []
    i = 0 #counter for status info
    print('Start Strand Bias Score calculation', file=sys.stderr)
    for k in qgram_annotate.keys():
        i += 1
        if i % 1000000 == 0: 
            print(" %s / %s Strand Bias Scores calculated" %(i, len(qgram_annotate.keys())))
        
        #get p-value for the strand bias table of q-gram k
        p_value = get_pvalue(qgram_annotate[k][0], qgram_annotate[k][1], qgram_annotate[k][2], qgram_annotate[k][3])
        
        #compute negative logarithm (base 10) of p-value or, if necessary, set to maxint
        strand_bias_score = sys.maxint if p_value < 1/10.0**300 else -math.log(p_value, 10)
        
        results.append((k, qgram_annotate[k][0], qgram_annotate[k][1], qgram_annotate[k][2], qgram_annotate[k][3], strand_bias_score))
    
    return results


def output(results, genome):
    """Output the results"""
    print("#Sequence", "Occurrence", "Forward Match", "Backward Match", "Forward Mismatch", "Backward Mismatch", "Strand Bias Score", "FER (Forward Error Rate)",
          "RER (Reverse Error Rate), ERD (Error rate Difference)", sep = '\t')
    
    for seq, forward_match, reverse_match, forward_mismatch, reverse_mismatch, sb_score, fer, rer, erd in results:
        occ = count(seq, genome)
        print(seq, occ, forward_match, reverse_match, forward_mismatch, reverse_mismatch, sb_score, fer, rer, erd, sep = '\t')


def get_motifspace_size(q,n):
    """return length of search space according to equation which is mentioned in Section 3.1 of the paper"""
    return reduce(lambda x, y: x + (int(sc.comb(q, y, exact=True)) * 4**(q-y)), [i for i in range(1, n+1)], int(sc.comb(q, 0, exact=True)) * 4**(q-0))

def count(qgram, genome):
    """Count number of q-grams and its reverse complement in genome"""
    rev = reverse_complement(qgram)
    rev = rev.replace('N', '.')
    qgram = qgram.replace('N', '.')
    
    return  len([m.start() for m in re.finditer(r'(?=(%s))' %qgram, genome)] + [m.start() for m in re.finditer(r'(?=(%s))' %rev, genome)])

def ident(genome, genome_annotate, q, n, alpha=0.05, epsilon=0.03, delta=0.05):
    """Identify critical <q>-grams (with <n> Ns) with reference to significance and error rate"""
    results = []
    
    motifspacesize_log = math.log(get_motifspace_size(q, n), 10)
    alpha_log = math.log(float(alpha), 10)
    
    qgram_annotate = get_annotate_qgram(genome, genome_annotate, q) #annotate each q-gram with Strand Bias Table
    add_n(qgram_annotate, n, q) #extend set of q-grams with q-grams containing Ns
    
    all_results = get_sb_score(qgram_annotate) #annotate each q-gram with Strand Bias Score

    sig_results = filter(lambda x: x[5] > motifspacesize_log - alpha_log, all_results) #filter statistically significant motifs (Bonferroni Correction)

    #filter motifs with regards to background error rate (epsilon) error rate difference (delta)
    for seq, forward_match, reverse_match, forward_mismatch, reverse_mismatch, p_value_score in sig_results:
        fer = float(forward_mismatch) / (forward_match + forward_mismatch) #forward error rate
        rer = float(reverse_mismatch) / (reverse_match + reverse_mismatch) #reverse error rate
        if rer < epsilon: #filter sequences with too high epsilon (background error rate)
            erd = fer - rer #error rate difference
            if erd >= delta: #filter sequences with too low delta (error rate difference cutoff)
                results.append((seq, forward_match, reverse_match, forward_mismatch, reverse_mismatch, p_value_score, fer, rer, erd)) 
            
    results.sort(key=lambda x: x[8],reverse=True) #sort by erd (error rate difference)
    
    output(results, genome)


if __name__ == '__main__':
    global NN
    parser = HelpfulOptionParser(usage=__doc__)
    
    parser.add_option("-a", dest="alpha", default=0.05, type="float", help="FWER (family-wise error rate) alpha, default: 0.05")
    parser.add_option("-e", dest="epsilon", default=0.03, type="float", help="background error rate cutoff epsilon, default: 0.03")
    parser.add_option("-d", dest="delta", default=0.05, type="float", help="error rate difference cutoff delta, default: 0.05")
    parser.add_option("-c", dest="learn_chrom", default="chr1", help="chromosome that is used to derive Context Specific Errors, default: chr1")
    parser.add_option("-v", dest="version", default=False, action="store_true", help="show script's version")
    
    (options, args) = parser.parse_args()
    
    if options.version:
            version = "version \"0.32\""
            print("")
            print(version)
            sys.exit()
    
    if len(args) != 4:
        parser.error("Sorry, exactly four parameters are required.")  
    
    #map arguments
    refpath = args[0]
    bampath = args[1]
    NN = int(args[2])
    n = int(args[3])
    #print (NN)
    genome, options.learn_chrom = get_genome(refpath, options.learn_chrom)
    fm, rm, fmm, rmm  = get_annotate_genome(genome, bampath, options.learn_chrom)
    nucls = [('A', 'C', 'G', 'T')]
    sets = []
    #for x in data:
    #    print (x)
    rows = []
    for i in range(NN):
        sets += nucls
    #print ("qgram, phred, phred_real, f occur, fm, fmm, fer, left, right, min")
    for s in itertools.product(*sets):
        for phred in range(0, 93):
            x = reduce(lambda x,y: x+y, s)
            key = (x, phred)
            #print (key, file = sys.stderr) 
            if fm[key] + fmm[key] > 0: #rm[key] + fmm[key] + rmm[key] > 0:
                focc = fm[key] + fmm[key]
                #rocc = rm[key] + rmm[key]
                phred_real = math.pow(10.0, -float(phred)/10.0)
                pvalf_left  = get_pbinom(fmm[key], focc, phred_real, True, True)
                pvalf_right  = get_pbinom(fmm[key], focc, phred_real, True, False)
                pvalf_lr = min(pvalf_left, pvalf_right)
                #pvalr_left  = get_pbinom(rmm[key], rocc, phred_real)
                #pvalr_right  = 1.0 - get_pbinom(rmm[key] - 1, rocc, phred_real)
                #pvalr_lr = min(pvalr_left, pvalr_right)
                #print (x, ",", phred, ",", phred_real, ",", focc, ",", fm[key],
#                          ",", fmm[key], ",",
#                          float(fmm[key])/focc if focc != 0 else "inf", ",", pvalf_left, ",", pvalf_right, ",", pvalf_lr#,
                          #",", rocc, ",", rm[key], 
                          #",", rmm[key], ",",
                          #float(rmm[key])/rocc if rocc != 0 else "inf", ",", pvalr_left, ",", pvalr_right, ",", pvalr_lr
#                       )
                rows.append((x, phred, phred_real, focc, fm[key], fmm[key], float(fmm[key])/focc if focc != 0 else "inf", pvalf_left, pvalf_right, pvalf_lr))
    with open("report.csv", "w") as f:
        for row in sorted(rows, key=lambda row: row[9]):
           print(reduce(lambda x,y: str(x) + "," + str(y), row), file=f)
    #for x in genome_annotate:
    #    print (x)
    #ident(genome, genome_annotate, q, n, options.alpha, options.epsilon, options.delta)

