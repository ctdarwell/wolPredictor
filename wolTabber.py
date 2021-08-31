#evaluate performance across ea iteration according to set thresholds
import pandas as pd
import numpy as np
import csv, sys

try:
    filename = sys.argv[1]
    td = sys.argv[2]
    mtchs = sys.argv[3]
    prfx = sys.argv[4]
except:
    filename = 'manual_correctedWolPreds_n200_sppDelims.csv' #'sppComm_correctedWolPreds_incr2x50_prg10.csv'
    td = 'manual_taxonDesignations_n200_sppDelims.csv' #'sppComm_taxonDesignations_2x50_prg10.csv'
    mtchs = 'manual_sppDelimMatches.csv' #'sppCommXelevationMatches.csv'
    prfx = 'manual'

taxDeg = pd.read_csv(td, header=0)
matches = pd.read_csv(mtchs, header=0)

with open(filename) as f:
    reader = csv.reader(f)
    columns = next(reader)
    colmap = dict(zip(range(len(columns)), columns))
f = None

dat2 = np.genfromtxt(filename, delimiter=',', dtype = 'U20', skip_header = 1)
wsp = np.unique(dat2[:, 1])


res = [] #TABLE: iteration, 'noWol' counts, matched strain counts
for i in range(2, dat2.shape[1]):
    tmp = []
    for w in wsp:
        qw = np.where(dat2[:, 1] == w)
        qw2 = np.where(dat2[:, i] == w)
        tmp.append(len(np.intersect1d(qw, qw2)))
    res.append([columns[i], tmp[0], sum(tmp[1:])])

new_res = []
ndiff, pdiff, tot, s = 0, 0, 0, 0
for line in res:
    if 'noPurge' in line[0]:
        if s == 1:
            div = res[res.index(line) - 1][0].split('_')[0]
            i1 = matches.ecoDegsCol[matches.ecoDegsCol.str.startswith(div)].iloc[0]
            i2 = int(matches.match[matches.ecoDegsCol.str.startswith(div)])
            new_res.append([i1, i2, ndiff, pdiff, tot, n, p])
        
        n, p = line[1], line[2]
        ndiff, pdiff, tot, s = 0, 0, 0, 0
    else:
        s = 1
        ndiff = max(ndiff, line[1] - n) #max noWol improv
        tot = max(tot, line[2] + line[1]) #max tot
        if line[2] == p: pdiff = max(pdiff, line[1] - n) #max noWol wo pos decr

div = res[res.index(line) - 1][0].split('_')[0]
i1 = matches.ecoDegsCol[matches.ecoDegsCol.str.startswith(div)].iloc[0]
i2 = int(matches.match[matches.ecoDegsCol.str.startswith(div)])
new_res.append([i1, i2, ndiff, pdiff, tot, n, p])

df = pd.DataFrame(new_res)
df.columns = ['taxDeg','match','pgeIncr','pgeIncr_wo_loss','maxScore','noPge_noWol','noPge_posWol']

def myfunc(a): return int(a.split('Spp')[-1])
vfunc = np.vectorize(myfunc)
try: df['match_nSpp'] = vfunc(df.taxDeg) #depending on version it adds an extra column indicating n species at each division
except: pass

df.to_csv(f'{prfx}_assocVars.csv', index=False)

print(f'Outputted: "{prfx}_assocVars.csv"')
