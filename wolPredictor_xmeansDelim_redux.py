#this version allows members of the same designated clade occupying different communities to have different matched strains
#spp delim via incremental cycling thru Xmeans

import numpy as np
from numpy import genfromtxt
import csv
import subprocess
import pandas as pd
from timeit import default_timer as timer
import itertools
from pyclustering.cluster.xmeans import xmeans
from pyclustering.cluster.center_initializer import kmeans_plusplus_initializer
import re

###### SET SOME PARAMETERS!!!!! ######
increment = 100 #how many divisions btwn upper & lower to split "species"
purge = 80 #upper % for examining purging
pge_incr = 10 #purge increment - default = 1
prefix = 'TEST' #add a prefix to file names - do NOT use 'x'
min_nSpp = 6 #min no. species from Xmeans - must be >=2
max_nSpp = 50 #max no. species from Xmeans
nSpp_incr = 1 #nSpp increment

wolbachia = 'wspClade' #hypo_wspClade
comm_column = 'community'
filename = 'testData.csv' #waspsExaDat_20190624 waspsRaxDat_20190624
tree = 'testPhylo.tre' # wasps20190624.exa.edit wasps20190624.rax.edit
NameOnPhylo = 'NameOnPhylo' #exaNames editTreeNames
out_dir = '.'

with open(filename) as f:
    reader = csv.reader(f)
    columns = next(reader)
    colmap = dict(zip(columns, range(len(columns))))
f = None

dat = genfromtxt(filename, delimiter=',', dtype = 'U20', skip_header = 1) # working CSV = datMLformatted_geiger2.csv
dat = dat[dat[:, colmap.get(NameOnPhylo)].argsort()]
taxa = np.squeeze(dat[:, [colmap.get(NameOnPhylo)]])
wols = np.unique(dat[:, [colmap.get(wolbachia)]])
wols = np.delete(wols, np.where(wols == 'noWol'), axis=0)
comms = np.unique(dat[:, [colmap.get(comm_column)]])
cols, taxaCols = ['taxa', 'wspClade'], ['taxa']
z = len(str(purge)) #fao CSV column name house-keeping (i.e. zfill length)
path2script = 'cophen4py.R'

def main():
    '''
    Main pgm control function
    '''
    if 'x' in prefix:
        print('RUN STOPPED: Do NOT use a small "x" in "prefix" variable!')
        return
    f = '{}-{}x{}_prg{}'.format(min_nSpp, max_nSpp, increment, purge) #STRING CONSTRUCTION!!!!!!!!!!!!!!!!
    print('Running "wolPredictor" - params:', prefix, f)
    R_cophen(tree, path2script)
    cophen = genfromtxt('{}_cophen.csv'.format(tree), dtype = 'float32', delimiter=',', skip_header = 1)
    x1, x2 = 0, 0 #calc array dims
    for _ in range(min_nSpp, max_nSpp, nSpp_incr): x1 += 1
    for _ in range(0, purge + 1, pge_incr): x2 += 1
    assigned, taxonDesignations, tmpDF = np.empty((len(taxa), (x1 * (x2 + 1)) + 2), dtype='U20' ), np.empty((len(taxa), (x1) + 1), dtype='U20' ), np.empty((len(taxa), 2), dtype='U20' )
    print('Matx shape', assigned.shape, 'TaxaShape', taxonDesignations.shape)
    assigned[:, 0] = np.squeeze(dat[:, [colmap.get(NameOnPhylo)]])
    assigned[:, 1] = np.squeeze(dat[:, [colmap.get(wolbachia)]])
    taxonDesignations[:, 0] = np.squeeze(dat[:, [colmap.get(NameOnPhylo)]])
    start = timer()

    for nSpp in range(min_nSpp, max_nSpp, nSpp_incr):
        if nSpp % (increment/10) == 0: print('Matrix iteration: ', str(nSpp), ' - time: ', str(timer() - start)[:-8])
        res = calcX(cophen, nSpp) #get optimal Xmeans cluster
        df = addPredict(res, nSpp, tmpDF) #spDelim is delimited clusters in actual host comms
        taxonDesignations[:, len(taxaCols) - 1] = res
        assigned[:, len(cols) - 1] = df[:, 1]
        tupl_df = tuple(df)
        for thresh2 in range(0, purge + 1, pge_incr): assigned = wolPurger(assigned, np.array(tupl_df), thresh2, nSpp, cophen)

    output1 = pd.DataFrame(assigned)
    output1.columns = [cols]
    output1.to_csv('{}/{}_wolPreds_incr{}.csv'.format(out_dir, prefix, f), index = False)
    output2 = pd.DataFrame(taxonDesignations)
    output2.columns = [taxaCols]
    output2.to_csv('{}/{}_taxonDesignations_{}.csv'.format(out_dir, prefix, f), index = False)
    
    assigned = matchStrains(assigned, taxonDesignations, start)
    
    output3 = pd.DataFrame(assigned)
    output3.columns = [cols]
    output3.to_csv('{}/{}_correctedWolPreds_incr{}.csv'.format(out_dir, prefix, f), index = False)
    print("Time:", str(timer() - start)[:-8])

def addPredict(res, nSpp, tmpDF):
    '''
    predict Wolbachia strain accoring to rules based on clusters on indvs within a community
    indvs will be given different strains if there are >=2 different taxon designations in a community 
    '''
    wolClades, tmpTaxa = [], []
    labsDik = dict(zip(np.squeeze(taxa), res))
    for comm in comms: #go thru fig host comms
        indvs = dat[np.where(dat[:, colmap.get(comm_column)] == comm)][:, colmap[NameOnPhylo]]
        [tmpTaxa.append(indv) for indv in indvs]
        clusters, grps = [], []
        [clusters.append(labsDik.get(indv)) for indv in indvs] #what barcode grps in that community?
        setClusts = list(set(clusters))
        if len(setClusts) > 1: [grps.append('{}_w{}'.format(comm[:3], str(setClusts.index(clusters[np.where(indvs == indv)[0][0]]) + 1))) for indv in indvs]
        else: grps = ['noWol'] * len(indvs)
        [wolClades.append(grp) for grp in grps]
        
    cols.append('nSpp{}_noPurge'.format(str(nSpp).zfill(z)))
    taxaCols.append('nSpp{}_nSpp{}'.format(str(nSpp).zfill(z), str(len(set(res)))))
    tmpDF[:, 0], tmpDF[:, 1] = np.array(tmpTaxa), np.array(wolClades)
    tmpDF = tmpDF[tmpDF[:, 0].argsort()]

    return tmpDF #prev returned assigned

def wolPurger(assigned, dfm, thresh2, nSpp, cophen):
    '''
    Decide whether divergence between species clusters warrants Wolbachia purging
        -at low thresh2 vals - most strains are purged (as most distances between clades are greater than thresh2)
        -at hi thresh2 vals - most are retained
    '''
    strainComms = list(set([x.split('_')[0] for x in dfm[:, 1]])) #all communities with strains
    strainComms = np.delete(strainComms, np.where(strainComms == 'noWol'), axis=0)
    to_purge, indv_strains, pairs = [], [], []
    for sc in strainComms: #eg ['arf', 'mic', 'tri']
        tmpStrains = []
        [tmpStrains.append(ss) for ss in dfm[:, 1] if sc in ss ] #matching assigned strains
        lst = list(itertools.combinations(set(tmpStrains), 2)) #all unique strains in this community
        for l in lst: #get indv clade pairs where WOL is NOT required
            grp1, grp2 = dfm[:, 0][(np.where((dfm[:, 1] == l[0])))], dfm[:, 0][(np.where((dfm[:, 1] == l[1])))]
            interDists = [cophen[np.where(dfm==g1)[0][0], np.where(dfm==g2)[0][0]] for g1 in grp1 for g2 in grp2]
            if min(interDists) >= thresh2/increment: indv_strains.append(l[0]), indv_strains.append(l[1]), pairs.append(sorted([l[0], l[1]])) #pairs indicates where wol is not required
 
    setStrains, sortPairs, strainComms = sorted(list(set(indv_strains))), [sorted(i) for i in pairs], list(set([x.split('_')[0] for x in indv_strains])) #strainComms = all communities with strains

    for sc in strainComms: #'tri', 'mic'
        matches = [string for string in setStrains if re.match(re.compile(sc + '.'), string)] #count no. clades/strains in this community
        for strain in setStrains: #['mic_w1', 'mic_w2', 'tri_w1', 'tri_w2', 'tri_w3']
            if sc in strain:
                cnt1 = sum([1 for pair in sortPairs if strain in pair]) #How many times does the clade/strain feature?
                if len(matches) - cnt1 < 2: to_purge.append(strain) #STILL NOT SURE ABOUT THIS!!!!

    for p in list(set(to_purge)): dfm[:, 1][np.where(dfm[:, 1] == p)] = 'noWol'

    cols.append('nSpp{}_Purge{}'.format(str(nSpp).zfill(z), str(thresh2).zfill(z)))
    assigned[:, len(cols) - 1] = dfm[:, 1]

    return assigned

def matchStrains(assigned, taxonDesignations, start):
    '''
    match the predictions with the empirically (yet arbritarily) named strains as much as possible
    ISSUE: the 2nd block may not find a solution (e.g. thresh16_noPurge):
        -because it does not remove previously attempted suboptimal solutions (attempted if optimal soln doesn't solve)
        -i.e. 'replace_with' is retained and attempted again at the next iteration if the solution doesn't work
        -currently fudged to give up at 20 attempts (basically a heuristic search)
    '''
    #identify best matches
    for column in range(2, assigned.shape[1]):
        tax_deg = [string for string in taxaCols if re.match(re.compile(cols[column].split('_')[0] + '.'), string)]
        df = np.concatenate([assigned[:, 1], taxonDesignations[:, taxaCols.index(tax_deg[0])], assigned[:, column]], axis = 0).reshape(3, len(taxa)).T
        unq, counts = np.unique(df, axis=0, return_counts=True)
        tab = np.column_stack((unq, counts))

        #select best wsp strain vs. predicted strain combinations
        sel_rows = []
        for wol in wols:
            indxs = np.where((tab[:, 0] == wol) & (tab[:, 2] != 'noWol')) #table indexes of wol not featuring 'noWol' #ex: np.where(tab == wol)
            if len(indxs[0]) == 0: continue
            tmpTab = np.squeeze(tab[indxs, :])
            if len(indxs[0]) == 1:
                sel_rows.append(tmpTab)
                continue
            for comm in comms:
                tmpIndxs = np.array([list(tmpTab[:, 2]).index(string) for string in list(tmpTab[:, 2]) if re.match(re.compile(comm[:3] + '.'), string)])
                if len(tmpIndxs) == 0: continue
                tmpTab2 = np.squeeze(tmpTab[tmpIndxs, :])
                if len(tmpIndxs) == 1:
                    sel_rows.append(tmpTab2)
                    continue
                b = tmpTab2[(np.where((tmpTab2[:, 0] == wol) & (tmpTab2[:, 2] != 'noWol')))[0], 3].astype(int) #counts of table of wol not featuring 'noWol' #ex: tab[np.where(tab == wol)[0], 3].astype(int)
                best = np.random.choice(np.flatnonzero(b == b.max()))
                sel_rows.append(tmpTab2[best, :])

        workTab = np.array(sel_rows)
        if len(workTab) < 2: continue

        #remove predicted strain clashes
        final_rows = []
        for strain in np.unique(workTab[:, 2]):
            tmpIndxs2 = np.where(workTab[:, 2] == strain)
            tmpTab3 = np.squeeze(workTab[tmpIndxs2, :])
            if len(tmpIndxs2[0]) == 1:
                final_rows.append(tmpTab3)
                continue
            b = tmpTab3[:, 3].astype(int)
            best = np.random.choice(np.flatnonzero(b == b.max()))
            final_rows.append(tmpTab3[best, :])

        finalTab = np.array(final_rows)

        for selected in finalTab:
            tba = assigned[:, 0][(assigned[:, 1] == selected[0]) & (assigned[:, column] == selected[2])]
            for taxon in tba: assigned[:, column][(assigned[:, 0] == taxon)] = selected[0] #all taxa in sub.ass with strain eg mic_w1

    return assigned

def R_cophen(tree, path2script):
    cmd = 'Rscript ' + path2script + ' ' + tree #Build subprocess command
    with open("PCPS_4_py.err", "wb") as err:
        subprocess.run(cmd, stderr=err)

def calcX(mat, nSpp):
    '''
    calculate clusters of individuals according to Xmeans
    '''
    initial_centers = kmeans_plusplus_initializer(mat, nSpp).initialize()
    xmeans_instance = xmeans(mat, initial_centers, nSpp)
    xmeans_instance.process()
    clusters = xmeans_instance.get_clusters()
    h, start, l = np.hstack(clusters), 0, [[] for _ in range(np.shape(mat)[0])]
    for cluster in clusters:
        for indv in h[start:start+len(cluster)]: l[indv] = clusters.index(cluster) #give groupings distinct cluster id's 
        start = start + len(cluster)

    return np.array(l)

if __name__ == '__main__':
    main()
