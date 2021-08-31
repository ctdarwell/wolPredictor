#This version allows the option to control whether...
#members of the same designated clade occupying different...
#communities have different matched strains.
#NB spp delim performed via file of pre-defined designations 

#THIS IS PARALLELISED 
#AROUND 3X QUICKER ON MY 4 CORE 16GB RAM LAPTOP

import numpy as np
from numpy import genfromtxt
import csv, subprocess, sys, getopt, re, itertools
import pandas as pd, collections as cl
from timeit import default_timer as timer 
import makePDF3 as mp
import concurrent.futures

###### SET SOME PARAMETERS!!!!! ######
cntrl = 0 #####
prefix = 'multi' #add a prefix to file names - do NOT use a small 'x'
pdf = 0 #switch on pdf maker
gap = 10 #gap btwn ticks on figure x-axis
nPges = 4 #no. of purges at each sp delim iteration
increment = 100 #how many divisions btwn upper & lower to split "species". LEAVE AS 100

spp_delim = 'randomDegs_n200.csv' #'tmpDegs.csv'
filename = 'exaData_faoReview.csv' #data file
wolbachia = 'wspClade' #hypo_wspClade
comm_column = 'sp.complex' #community column in data file
NameOnPhylo = 'taxa' #phylogeny taxa names column
tree = 'exaData_faoReview.tre' # tree file
dat_dir, out_dir = '.', '.' #output directory
path2script = 'cophen4pyOut.R' #R script

#read any new params
if '-h' in sys.argv:
    print("\n --Flag: Variable Name - Explanation - Default Value\n\n '-p': prefix - add a prefix to file names - default = '{}' (NB do NOT use a small 'x')\n '-r': spp_delim - filename of pre-defined taxon designations - default = '{}'\n '-f': filename - main input data file - default = '{}'\n '-w': wolbachia - column name for empirically derived Wolbachia strains - default = '{}'\n '-c': comm_column - column name for host communities - default = '{}'\n '-n': NameOnPhylo - column name for sample individual names - default = '{}'\n '-t': tree - input phylogenetic tree - default = '{}'\n '-o': out_dir - output directory - default = '{}'\n '-d': dat_dir - data files directory - default = '{}'\n '-q': pdf - make figure (Off/On: 0/1) - default = '{}'\n '-g': gap between ticks on figure x-axis - default = '{}'\n '-C': control strains in multiple communities (Off/On: 0/1) - default = '{}'\n '-N': No. of purges at each species delim iteration (max=10) - default = '{}'\n"
          .format(prefix, spp_delim, filename, wolbachia, comm_column, NameOnPhylo, tree, out_dir, dat_dir, pdf, gap, cntrl, nPges))
    sys.exit(2)
try: #input sys args
    myopts, args = getopt.getopt(sys.argv[1:],"p:r:w:c:f:t:n:o:d:q:g:C:N:")
except:
    print('I/O ERROR!')
    sys.exit(2)

#Sort non-default params and check for input errors
print('\n### Check inputted params ###')
for a, b in myopts:
    print(a, b)
    if a == '-p': prefix = b
    if a == '-r': spp_delim = b
    if a == '-w': wolbachia = b
    if a == '-c': comm_column = b
    if a == '-f': filename = b
    if a == '-t': tree = b
    if a == '-n': NameOnPhylo = b
    if a == '-o': out_dir = b
    if a == '-d': dat_dir = b
    if a == '-q': pdf = int(b)
    if a == '-g': gap = int(b)
    if a == '-C': cntrl = int(b)
    if a == '-N': nPges = int(b)

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

if nPges > 10:
    print("Variable 'nPges' is greater than 10.......terminating!")
    sys.exit(2)

#read main data file and organise data
with open('{}/{}'.format(dat_dir, filename)) as f:
    reader = csv.reader(f)
    columns = next(reader)
    colmap = dict(zip(columns, range(len(columns))))
f, myopts = None, None

dat = genfromtxt('{}/{}'.format(dat_dir, filename), delimiter=',', dtype = 'U20', skip_header = 1) # working CSV = datMLformatted_geiger2.csv
dat = dat[dat[:, colmap.get(NameOnPhylo)].argsort()]

#get taxDeg colNames
with open('{}/{}'.format(dat_dir, spp_delim)) as f2:
    reader = csv.reader(f2)
    degCols = next(reader)

taxonDesignations = genfromtxt('{}/{}'.format(dat_dir, spp_delim), delimiter=',', dtype = 'U20', skip_header = 1) # working CSV = datMLformatted_geiger2.csv
taxonDesignations = taxonDesignations[taxonDesignations[:, 0].argsort()]

taxa = np.squeeze(dat[:, [colmap.get(NameOnPhylo)]])
wols = np.unique(dat[:, [colmap.get(wolbachia)]])
wols = np.delete(wols, np.where(wols == 'noWol'), axis=0)
comms = np.unique(dat[:, [colmap.get(comm_column)]])
cols, taxaCols = cl.deque(['taxa', 'wspClade']), cl.deque(['taxa'])

newDegCols = ['taxa'] #make an additional list for taxon degs column names
for deg in degCols[1:]: newDegCols.append(f"{deg}_nSpp{np.unique(taxonDesignations[:, degCols.index(deg)]).shape[0]}")

#vfuncs fao addPredict/wolPurger
def f1(comm, x): return f"{comm[:3]}_{x}"
vf1 = np.vectorize(f1, excluded=['comm']) #, otypes=['O']

def f2(x): return x.split('_')[0]
vf2 = np.vectorize(f2)

#Main control loop
def main():
    '''
    Main pgm control function
    '''
    if 'x' in prefix:
        print('RUN STOPPED: Do NOT use a lower case "x" in "prefix" variable!')
        return

    #get phylo, make df
    cophen, purge, pge_incr, z, _ = R_cophen('{}/{}'.format(dat_dir, tree), path2script) #get Dist Mat from phylogeny in R
    f = f"n{taxonDesignations.shape[1] - 1}_sppDelims" #'{}x{}_prg{}'.format(min_nSpp, max_nSpp - 1, purge) #original max_nSpp val
    print(f'\nRunning "wolPredictor_MANUAL" - params: {prefix}_{f}')

    assigned, tmpDF = np.empty((len(taxa), ((taxonDesignations.shape[1] - 1) * (1 + _)) + 2), dtype='U20' ), np.empty((len(taxa), 2), dtype='U20' )
    print('Matx shape', assigned.shape, 'TaxaShape', taxonDesignations.shape)

    assigned[:, 0] = np.squeeze(dat[:, [colmap.get(NameOnPhylo)]])
    assigned[:, 1] = np.squeeze(dat[:, [colmap.get(wolbachia)]])

    start = timer()

    #run addPredict
    with concurrent.futures.ProcessPoolExecutor() as executor:
        clades = [taxonDesignations[:, delim] for delim in range(1, taxonDesignations.shape[1])]
        results = executor.map(addPredict, clades, list(range(1, taxonDesignations.shape[1])), [tmpDF]*len(clades))

    #run wolPurger
    with concurrent.futures.ProcessPoolExecutor() as executor:
        res = [result for result in results] #results from addPredict
        thresh2 = [[thresh2 for thresh2 in range(pge_incr, purge + 1, pge_incr)]] * len(res)
        pgeRes = executor.map(wolPurger, res, thresh2, [cophen] * len(res)) #purge wol at incremental thresholds at this spp delim

    purges = [p for p in pgeRes] #results from wolPurger
    
    #build up assigned dataframe
    cntr, cntr2, cntr3 = 2, 3, 0
    for r in res:
        assigned[:, cntr] = r
        assigned[:, cntr2:cntr2 + len(thresh2[0])] = purges[cntr3]
        cntr2 += len(thresh2[0]) + 1
        cntr += len(thresh2[0]) + 1
        cntr3 += 1

    for qw in degCols[1:]:
        cols.append(f"{qw}_noPurge")
        for thresh in thresh2[0]:
            cols.append(f"{qw}_Purge{str(thresh).zfill(z)}")

    #write predictions & taxon delim file
    output1 = pd.DataFrame(assigned)
    output1.columns = cols
    output1.to_csv('{}/{}_wolPreds_{}.csv'.format(out_dir, prefix, f), index = False)

    output2 = pd.DataFrame(taxonDesignations)
    output2.columns = newDegCols
    output2.to_csv('{}/{}_taxonDesignations_{}.csv'.format(out_dir, prefix, f), index = False)

    print("Finished addPred/wolPge:", str(timer() - start))
    
    #run matchStrains
    with concurrent.futures.ProcessPoolExecutor() as executor:
        asses = [assigned[:, x] for x in range(2, assigned.shape[1])]
        assCols = [z.split('_')[0] for z in cols][2:]
        matches = executor.map(matchStrains, asses, assCols) #purge wol at incremental thresholds at this spp delim

    #update assigned dataframe with matches 
    mtchRes = [m for m in matches]
    
    cntr = 2
    for mtch in mtchRes:
        assigned[:, cntr] = mtch
        cntr += 1
    
    #write final matched df
    output3 = pd.DataFrame(assigned)
    output3.columns = cols
    output3.to_csv('{}/{}_correctedWolPreds_{}.csv'.format(out_dir, prefix, f), index = False)

    if pdf == 1: #y/n write results figures
        tix = 0
        print("Analyses complete. Making figures....")
        mp.makePDF('{}/{}_correctedWolPreds_{}.csv'.format(out_dir, prefix, f), gap, tix, 'div')
        print('\n\t....Figures outputted to "', out_dir, '"\n', sep = '')

    print("Time:", str(timer() - start))

def R_cophen(tree, path2script):
    '''
    Build cophenetic distance matrix from inputted tree and decide purging parameters according to max(cophen)
    '''
    cmd = 'Rscript ' + path2script + ' ' + tree #Build subprocess command
    with open("PCPS_4_py.err", "wb") as err: #call/run R script
        out = subprocess.check_output(cmd, shell=False, stderr=err).decode("utf-8")
 
    cophen = genfromtxt(out.split('\r\n')[1:-1], dtype = 'float32', delimiter=',') #get the distance matrix from R call
    purge = int(np.max(cophen) * 100) + 1 #calculate some house-keeping variables regarding purging
    if purge >= 10: pge_incr = int(purge / nPges)
    else: pge_incr = 2
    cnt = 0
    for i in range(pge_incr, purge + 1, pge_incr): cnt += 1
    z = len(str(purge)) #fao CSV column name house-keeping (i.e. zfill length)

    return cophen, purge, pge_incr, z, cnt

def addPredict(res, nSpp, tmpDF):
    '''
    predict Wolbachia strain accoring to rules based on clusters on indvs within a community
    indvs will be given different strains if there are >=2 different taxon designations in a community 
    '''
    wolClades, tmpTaxa = cl.deque([]), cl.deque([])
    for comm in comms: #go thru fig host comms
        indvs = dat[np.where(dat[:, colmap.get(comm_column)] == comm)][:, colmap[NameOnPhylo]] #taxon names in the comm
        tmpTaxa.extend(indvs)
        clusters = res[np.isin(taxa, indvs)] #get taxon delims for those indvs
        if np.unique(clusters).shape[0] > 1: grps = vf1(comm=comm, x=clusters) #if >1 species assign putative Wolbachia associations
        else: grps = np.full((indvs.shape[0]), 'noWol') #else assign 'noWol'
        wolClades.extend(grps) 
        
    tmpDF[:, 0], tmpDF[:, 1] = np.array(tmpTaxa), np.array(wolClades)
    tmpDF = tmpDF[tmpDF[:, 0].argsort()] #resort the df by taxon name

    return tmpDF[:, 1] #prev returned assigned

def wolPurger(df, thresh2, cophen):
    '''
        Decide whether divergence between species clusters warrants Wolbachia purging
            -at low thresh2 vals - most strains are purged (as most distances between clades are greater than thresh2)
            -at hi thresh2 vals - most are retained
    '''
    pgeDF = np.vstack([df] * len(thresh2)).T #create array of predictions (from addPredict)large enough for all purge thresholds
    
    for thresh in thresh2: #at each thresh
        dfm = np.column_stack((taxa, pgeDF[:, thresh2.index(thresh)])) #tmp array of taxa & predictions
        strainComms = np.unique(np.array(list(map(vf2, dfm[:,1])))) #id all comms with assignments
        strainComms = np.delete(strainComms, np.where(strainComms == 'noWol'), axis=0) #rm noWol (ie pos assigns only) 
        to_purge, indv_strains, pairs = cl.deque([]), cl.deque([]), cl.deque([])
        for sc in strainComms: #eg ['arf', 'mic', 'tri']
            tmpStrains = cl.deque([])
            [tmpStrains.append(ss) for ss in dfm[:, 1] if sc in ss ] #matching assigned strains
            lst = list(itertools.combinations(set(tmpStrains), 2)) #all unique strains in this community
            for l in lst: #get indv clade pairs where WOL is NOT required
                grp1, grp2 = dfm[:, 0][(np.where((dfm[:, 1] == l[0])))], dfm[:, 0][(np.where((dfm[:, 1] == l[1])))]
                interDists = [cophen[np.where(dfm==g1)[0][0], np.where(dfm==g2)[0][0]] for g1 in grp1 for g2 in grp2]
                if min(interDists) >= thresh/increment: indv_strains.append(l[0]), indv_strains.append(l[1]), pairs.append(sorted([l[0], l[1]])) #pairs indicates where wol is not required
     
        setStrains, sortPairs, strainComms = sorted(list(set(indv_strains))), [sorted(i) for i in pairs], list(set([x.split('_')[0] for x in indv_strains])) #strainComms = all communities with strains
    
        for sc in strainComms: #'tri', 'mic'
            matches = [string for string in setStrains if re.match(re.compile(sc + '.'), string)] #count no. clades/strains in this community
            for strain in setStrains: #['mic_w1', 'mic_w2', 'tri_w1', 'tri_w2', 'tri_w3']
                if sc in strain:
                    cnt1 = sum([1 for pair in sortPairs if strain in pair]) #How many times does the clade/strain feature?
                    if len(matches) - cnt1 < 2: to_purge.append(strain)

        for p in list(set(to_purge)): pgeDF[:, thresh2.index(thresh)][np.where(dfm[:, 1] == p)] = 'noWol'

    return pgeDF

def matchStrains(pred, col):
    '''
    match the predictions with the empirically (yet arbritarily) named strains as much as possible
    ISSUE: the 2nd block may not find a solution (e.g. thresh16_noPurge):
        -because it does not remove previously attempted suboptimal solutions (attempted if optimal soln doesn't solve)
        -i.e. 'replace_with' is retained and attempted again at the next iteration if the solution doesn't work
        -currently fudged to give up at 20 attempts (basically a heuristic search)
    '''
    #identify best matches - concatenate wspClade/species delim/predictions in one d.f.
    df = np.concatenate([dat[:, colmap.get(wolbachia)], taxonDesignations[:, degCols.index(col)], pred], axis = 0).reshape(3, len(taxa)).T
    unq, counts = np.unique(df, axis=0, return_counts=True)
    tab = np.column_stack((unq, counts)) #get counts of all data combinations

    #select best wsp strain vs. predicted strain combinations
    sel_rows = cl.deque([])
    for wol in wols: #thru all wspClades
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
    if len(workTab) < 2: return pred

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
        tba = taxa[(df[:, 0] == selected[0]) & (df[:, 2] == selected[2])]
        for taxon in tba: df[:, 2][(taxa == taxon)] = selected[0] #all taxa in sub.ass with strain eg mic_w1

    return df[:, 2]

if __name__ == '__main__': main()
