import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.axes as ax
import sys
import math

def common_entries(*dcts):
    for i in set(dcts[0]).intersection(*dcts[1:]):
        yield (i,) + tuple(d[i] for d in dcts)

def readGramsWithValues(filename):
    map = dict()
    lines = 0
    with open(filename) as f:
        for line in f:
            row = line.strip().split(",")
            qgram = row[0]
            phred = row[1]
            values = row[6:]
            map[(qgram, phred)] = values
            lines += 1
#           print (qgram, values)
    print ("lines read:", lines)
    return map

if __name__ == '__main__':

    args = sys.argv
    filename_out = args[3]
    useLog = int(args[4]) if len(args) > 4 else 0
    koeff = float(args[5]) if len(args) > 5 else 10000
    off = float(args[6]) if len(args) > 6 else 10000
    sign = int(args[7]) if len(args) > 7 else 1 
    coords = list(common_entries(readGramsWithValues(args[1]), readGramsWithValues(args[2])))
#    for x in coords:
#        print x
    #N = 50
    #fer = 0
    #p_left = 1
    #p_r = 2
    #p_min = 3
    fer1 = [float(row[1][0]) for row in coords]
    fer2 = [float(row[2][0]) for row in coords]
    log_fer1 = [0.01 if f < 1/10.0**100 else math.log(f, 10) for f in fer1]
    log_fer2 = [0.01 if f < 1/10.0**100 else math.log(f, 10) for f in fer2]
    x = [f for f in log_fer1] if useLog == 1 else [f for f in fer1]
    y = [f for f in log_fer2] if useLog == 1 else [f for f in fer2] #math.log(float(row[2][0]), 10) for row in coords]
    outliers_x = []
    outliers_y = []
    i = 0
    if koeff != 10000 and off != 10000:
        while i < len(x):
            #if (y[i] - off)/x[i] > koeff:
            #    print ("ololo", sign) 
            if (y[i]-off)/x[i] * sign > koeff * sign:
                outliers_x.append(x[i])
                outliers_y.append(y[i])
                print (coords[i][0], coords[i][1], coords[i][2], x[i], y[i])
                #x.pop(i)
                #y.pop(i)
                #i -= 1
            i += 1
    
    #print ("total points", len(x), len(y))
#    for q,w in zip(x,y):
#        print (q,w)    
#print (x,y)
    #area = np.pi * (15 * np.random.rand(N))**2 # 0 to 15 point radiuses

    #plt.autoscale(False)
    #plt.set_aspect('equal')
    if koeff == 10000 and off == 10000:
        plt.scatter(x, y, s= 2, alpha = 0.5)
    else:
        plt.plot(x, y, 'bo', outliers_x, outliers_y, "go", alpha=0.5)
#    plt.yscale('log')
#    plt.xscale('log')
    plt.axis('equal')
    #curr = plt.axis()
    #print (curr)
    #plt.axis([0, 100000, 0, curr[3]])
    #print (plt.axis())
    #plt.axis.set_xlim(left=0)
    #plt.axes().set_aspect('equal', 'box')
    #plt.autoscale(False)
    plt.savefig(filename_out)

