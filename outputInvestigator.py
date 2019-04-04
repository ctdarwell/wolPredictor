#sys: 1-filename - root = correctedWolPreds 2-min. no. of negative strains 3-min. no. of positive strains

import numpy as np
from numpy import genfromtxt
import csv, sys

filename = sys.argv[1]
with open(filename) as f:
    reader = csv.reader(f)
    columns = next(reader)
    colmap = dict(zip(range(len(columns)), columns))
f = None

dat2 = genfromtxt(filename, delimiter=',', dtype = 'U20', skip_header = 1) 

wsp = np.unique(dat2[:, 1])

res = [] #TABLE: iteration, 'noWol' counts, matched strain counts
for i in range(2, len(colmap)):
    tmp = []
    for w in wsp:
        qw = np.where(dat2[:, 1] == w)
        qw2 = np.where(dat2[:, i] == w)
        tmp.append(len(np.intersect1d(qw, qw2)))
    if tmp[0] > int(sys.argv[2]) and sum(tmp[1:]) > int(sys.argv[3]): res.append([colmap.get(i), tmp[0], sum(tmp[1:])])

print('[--- iteration --- //noWol//matched]')
for r in res: print(r)
