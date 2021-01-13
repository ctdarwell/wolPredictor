#This version allows the option to control whether 
#members of the same designated clade occupying different
#communities have different matched strains
#NB spp delim performed via incremental cycling thru Xmeans

import numpy as np
from numpy import genfromtxt
import csv, subprocess, sys, getopt, re, itertools
import pandas as pd, collections as cl
from timeit import default_timer as timer 
from pyclustering.cluster.xmeans import xmeans
from pyclustering.cluster.center_initializer import kmeans_plusplus_initializer
import makePDF as mp

###### SET SOME PARAMETERS!!!!! ######
shuffle = 0 #####
cntrl = 0 #####
prefix = 'TEST' #add a prefix to file names - do NOT use a small 'x'
min_nSpp = 2 #min no. species from Xmeans - must be >=2
max_nSpp = 50 #max no. species from Xmeans
nSpp_incr = 1 #nSpp increment
increment = 100 #how many divisions btwn upper & lower to split "species". LEAVE AS 100
pdf = 1 #switch on pdf maker
gap = 10 #gap btwn ticks on figure x-axis

filename = 'testData.csv' #data file
wolbachia = 'wspClade' #hypo_wspClade
comm_column = 'community' #community column in data file
NameOnPhylo = 'NameOnPhylo' #phylogeny taxa names column
tree = 'testPhylo.tre' # tree file
dat_dir, out_dir = '.', '.' #output directory
path2script = 'cophen4pyOut.R' #R script

#read any new params
if '-h' in sys.argv:
    print("\n --Flag: Variable Name - Explanation - Default Value\n\n '-p': prefix - add a prefix to file names - default = '{}' (NB do NOT use a small 'x')\n '-m': min_nSpp - min no. species to iterate thru Xmeans; NB -m >= 2 - default = '{}'\n '-M': max_nSpp - max no. species to iterate thru Xmeans - default = '{}'\n '-i': nSpp_incr - incremental increase btwn min. & max. Xmeans 'number of species' clustering - default = '{}'\n '-f': filename - main input data file - default = '{}'\n '-w': wolbachia - column name for empirically derived Wolbachia strains - default = '{}'\n '-c': comm_column - column name for host communities - default = '{}'\n '-n': NameOnPhylo - column name for sample individual names - default = '{}'\n '-t': tree - input phylogenetic tree - default = '{}'\n '-o': out_dir - output directory - default = '{}'\n '-d': dat_dir - data files directory - default = '{}'\n '-q': pdf - make figure (Off/On: 0/1) - default = '{}'\n '-g': gap between ticks on figure x-axis - default = '{}'\n '-s': shuffle wsp clades (Off/On: 0/1) - default = '{}'\n '-C': control strains in multiple communities (Off/On: 0/1) - default = '{}'\n"
          .format(prefix, min_nSpp, max_nSpp, nSpp_incr, filename, wolbachia, comm_column, NameOnPhylo, tree, out_dir, dat_dir, pdf, gap, shuffle, cntrl))
    sys.exit(2)
try: #input sys args
    myopts, args = getopt.getopt(sys.argv[1:],"p:m:M:i:w:c:f:t:n:o:d:q:g:s:C:")
except:
    print('I/O ERROR!')
    sys.exit(2)

#Sort non-default params and check for input errors
print('\n### Check inputted params ###')
for a, b in myopts:
    print(a, b)
    if a == '-p': prefix = b
    if a == '-m': min_nSpp = int(b)
    if a == '-M': max_nSpp = int(b)
    if a == '-i': nSpp_incr = int(b)
    if a == '-w': wolbachia = b
    if a == '-c': comm_column = b
    if a == '-f': filename = b
    if a == '-t': tree = b
    if a == '-n': NameOnPhylo = b
    if a == '-o': out_dir = b
    if a == '-d': dat_dir = b
    if a == '-q': pdf = int(b)
    if a == '-g': gap = int(b)
    if a == '-s': shuffle = int(b)
    if a == '-C': cntrl = int(b)

swch = 0
for sysarg in sys.argv[1:]:
    if sysarg[0] == '-': swch = 1
    if sysarg[0] != '-' and swch == 1:
        swch = 0
        continue
    if sysarg[0] != '-' and swch == 0:
        print('No flag added for: ', sysarg, sep = '')
        swch = 0
        sys.exit(2)

#read main data file and organise data
with open('{}/{}'.format(dat_dir, filename)) as f:
    reader = csv.reader(f)
    columns = next(reader)
    colmap = dict(zip(columns, range(len(columns))))
f, myopts = None, None

dat = genfromtxt('{}/{}'.format(dat_dir, filename), delimiter=',', dtype = 'U20', skip_header = 1) # working CSV = datMLformatted_geiger2.csv
dat = dat[dat[:, colmap.get(NameOnPhylo)].argsort()]
taxa = np.squeeze(dat[:, [colmap.get(NameOnPhylo)]])
wols = np.unique(dat[:, [colmap.get(wolbachia)]])
wols = np.delete(wols, np.where(wols == 'noWol'), axis=0)
comms = np.unique(dat[:, [colmap.get(comm_column)]])
cols, taxaCols = cl.deque(['taxa', 'wspClade']), cl.deque(['taxa'])
max_nSpp += 1

if shuffle == 1: dat[:, colmap.get(wolbachia)] = np.random.permutation(dat[:, colmap.get(wolbachia)]) #randomly shuffle wsp strains

#Main control loop
def main():
    '''
    Main pgm control function
    '''
    if 'x' in prefix:
        print('RUN STOPPED: Do NOT use a small "x" in "prefix" variable!')
        return
    cophen, purge, pge_incr, z = R_cophen('{}/{}'.format(dat_dir, tree), path2script) #get Dist Mat from phylogeny in R
    f = '{}x{}_prg{}'.format(min_nSpp, max_nSpp - 1, purge) #original max_nSpp val
    print(f'\nRunning "wolPredictor_xmeansDelim" - params: {prefix} {f}')
    x1, x2 = 0, 0 #calc array dims then initialise arrays
    for _ in range(min_nSpp, max_nSpp, nSpp_incr): x1 += 1
    for _ in range(pge_incr, purge + 1, pge_incr): x2 += 1
    assigned, taxonDesignations, tmpDF = np.empty((len(taxa), (x1 * (x2 + 1)) + 2), dtype='U20' ), np.empty((len(taxa), (x1) + 1), dtype='U20' ), np.empty((len(taxa), 2), dtype='U20' )
    print('Matx shape', assigned.shape, 'TaxaShape', taxonDesignations.shape)
    assigned[:, 0] = np.squeeze(dat[:, [colmap.get(NameOnPhylo)]])
    assigned[:, 1] = np.squeeze(dat[:, [colmap.get(wolbachia)]])
    taxonDesignations[:, 0] = np.squeeze(dat[:, [colmap.get(NameOnPhylo)]])
    start = timer()

    for nSpp in range(min_nSpp, max_nSpp, nSpp_incr): #Main pgm loop
        if nSpp % (increment/10) == 0: print('Matrix iteration: ', str(nSpp), ' - time: ', str(timer() - start)[:-8])
        res = calcX(cophen, nSpp) #get optimal Xmeans cluster (i.e. spp delim)
        df = addPredict(res, nSpp, tmpDF, z) #add wol strain predictions at this spp delim level
        taxonDesignations[:, len(taxaCols) - 1] = res #add spp delim to df
        assigned[:, len(cols) - 1] = df[:, 1] #add predictions to df
        tupl_df = tuple(df)
        for thresh2 in range(pge_incr, purge + 1, pge_incr): assigned = wolPurger(assigned, np.array(tupl_df), thresh2, nSpp, cophen, z) #purge wol at incremental thresholds at this spp delim

    output1 = pd.DataFrame(assigned)
    output1.columns = cols
    output1.to_csv('{}/{}_wolPreds_incr{}.csv'.format(out_dir, prefix, f), index = False)
    output2 = pd.DataFrame(taxonDesignations)
    output2.columns = taxaCols
    output2.to_csv('{}/{}_taxonDesignations_{}.csv'.format(out_dir, prefix, f), index = False)
        
    assigned = matchStrains(assigned, taxonDesignations, start) #match predictions with empirical data
    
    output3 = pd.DataFrame(assigned)
    output3.columns = cols
    output3.to_csv('{}/{}_correctedWolPreds_incr{}.csv'.format(out_dir, prefix, f), index = False)
    if pdf == 1: #y/n write results figures
        mp.makePDF('{}/{}_correctedWolPreds_incr{}.csv'.format(out_dir, prefix, f), gap)
        print('\nFigures outputted to "', out_dir, '"\n', sep = '')
    print("Time:", str(timer() - start)[:-8])

def R_cophen(tree, path2script):
    '''
    Build cophenetic distance matrix from inputted tree and decide purging parameters according to max(cophen)
    '''
    cmd = 'Rscript ' + path2script + ' ' + tree #Build subprocess command
    with open("PCPS_4_py.err", "wb") as err:
        out = subprocess.check_output(cmd, shell=False, stderr=err).decode("utf-8")
 
    cophen = genfromtxt(out.split('\r\n')[1:-1], dtype = 'float32', delimiter=',')
    purge = int(np.max(cophen) * 100) + 1
    pge_incr = int(purge/6)
    z = len(str(purge)) #fao CSV column name house-keeping (i.e. zfill length)

    return cophen, purge, pge_incr, z

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

def addPredict(res, nSpp, tmpDF, z):
    '''
    predict Wolbachia strain accoring to rules based on clusters on indvs within a community
    indvs will be given different strains if there are >=2 different taxon designations in a community 
    '''
    wolClades, tmpTaxa = cl.deque([]), cl.deque([])
    labsDik = dict(zip(np.squeeze(taxa), res))
    for comm in comms: #go thru fig host comms
        indvs = dat[np.where(dat[:, colmap.get(comm_column)] == comm)][:, colmap[NameOnPhylo]]
        [tmpTaxa.append(indv) for indv in indvs]
        clusters, grps = cl.deque([]), cl.deque([])
        [clusters.append(labsDik.get(indv)) for indv in indvs] #what barcode grps in that community?
        setClusts = list(set(clusters))
        if len(setClusts) > 1: [grps.append('{}_w{}'.format(comm[:3], str(setClusts.index(clusters[np.where(indvs == indv)[0][0]]) + 1))) for indv in indvs]
        else: grps = ['noWol'] * len(indvs)
        [wolClades.append(grp) for grp in grps]
        
    cols.append('nSpp{}_noPurge'.format(str(nSpp).zfill(z))) #!!!
    taxaCols.append('nSpp{}_nSpp{}'.format(str(nSpp).zfill(z), str(len(set(res))))) #!!!
    tmpDF[:, 0], tmpDF[:, 1] = np.array(tmpTaxa), np.array(wolClades)
    tmpDF = tmpDF[tmpDF[:, 0].argsort()]

    return tmpDF #prev returned assigned

def wolPurger(assigned, dfm, thresh2, nSpp, cophen, z):
    '''
    Decide whether divergence between species clusters warrants Wolbachia purging
        -at low thresh2 vals - most strains are purged (as most distances between clades are greater than thresh2)
        -at hi thresh2 vals - most are retained
    '''
    strainComms = list(set([x.split('_')[0] for x in dfm[:, 1]])) #all communities with strains
    strainComms = np.delete(strainComms, np.where(strainComms == 'noWol'), axis=0)
    to_purge, indv_strains, pairs = cl.deque([]), cl.deque([]), cl.deque([])
    for sc in strainComms: #eg ['arf', 'mic', 'tri']
        tmpStrains = cl.deque([])
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
                if len(matches) - cnt1 < 2: to_purge.append(strain)

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
        sel_rows = cl.deque([])
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
        final_rows = cl.deque([])
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
        
        if cntrl == 1: #remove taxon clashes
            ult_rows = cl.deque([])
            for taxon in np.unique(finalTab[:, 1]):
                tmpIndxs3 = np.where(finalTab[:, 1] == taxon)
                tmpTab4 = np.squeeze(finalTab[tmpIndxs3, :])
                if len(tmpIndxs3[0]) == 1:
                    ult_rows.append(tmpTab4)
                    continue
                b = tmpTab4[:, 3].astype(int)
                best = np.random.choice(np.flatnonzero(b == b.max()))
                ult_rows.append(tmpTab4[best, :])

            finalTab = np.array(ult_rows)

        for selected in finalTab:
            tba = assigned[:, 0][(assigned[:, 1] == selected[0]) & (assigned[:, column] == selected[2])]
            for taxon in tba: assigned[:, column][(assigned[:, 0] == taxon)] = selected[0] #all taxa in sub.ass with strain eg mic_w1

    return assigned

if __name__ == '__main__':
    main()
