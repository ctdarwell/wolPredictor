#evaluate performance across ea iteration according to set thresholds
import pandas as pd
import numpy as np
from numpy import genfromtxt
import csv, re, sys

filename = sys.argv[1]
taxDeg = pd.read_csv(sys.argv[2], header=0)
with open(filename) as f:
    reader = csv.reader(f)
    columns = next(reader)
    colmap = dict(zip(range(len(columns)), columns))
f = None

dat2 = genfromtxt(filename, delimiter=',', dtype = 'U20', skip_header = 1) 
wsp = np.unique(dat2[:, 1])
#g=open('outFile.csv','w') #uncomment out (and lines 33 & 35) to write file

noWols, posStr, tot = sys.argv[3], sys.argv[4], sys.argv[5]  # RECOMMENDED set noWols = -1, pos strains = <positive value>, total = <positive value>
print('[ iteration //noWol//match//tot//%]')
res = [] #TABLE: iteration, 'noWol' counts, matched strain counts
for i in range(2, len(colmap)):
    tmp = []
    for w in wsp:
        qw = np.where(dat2[:, 1] == w)
        qw2 = np.where(dat2[:, i] == w)
        tmp.append(len(np.intersect1d(qw, qw2)))
    noWol = np.where(wsp == 'noWol')[0][0]
    if tmp[noWol] > noWols and sum(tmp[:noWol]) + sum(tmp[noWol + 1:]) > posStr and tmp[noWol] + sum(tmp[:noWol]) + sum(tmp[noWol + 1:]) > tot:
        tax_deg = [string for string in list(taxDeg) if re.match(re.compile(colmap.get(i).split('_')[0] + '.'), string)]
        nSpp = tax_deg[0].split('_nSpp')[1]
        print([colmap.get(i), tmp[noWol], sum(tmp[:noWol]) + sum(tmp[noWol + 1:]), sum(tmp), float("{0:.2f}".format(100*sum(tmp)/np.shape(taxDeg)[0])), 'nSpecies = ', nSpp])
        #g.write(str(colmap.get(i)) + ',' + str(tmp[noWol]) + ',' + str(sum(tmp[:noWol]) + sum(tmp[noWol + 1:])) + ',' + str(sum(tmp)) + '\n')

#g.close()
