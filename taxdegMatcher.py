import pandas as pd, sys
import warnings
warnings.simplefilter(action='ignore', category=Warning)

try:
    filename, ecofile = sys.argv[1], sys.argv[2]
    taxa, ecotaxa, arbDegs = sys.argv[3], sys.argv[4], sys.argv[5]
except:
    filename, ecofile = 'exaData_faoReview.csv', 'manual_taxonDesignations_n200_sppDelims.csv'
    taxa, ecotaxa, arbDegs = 'taxa', 'taxa', 'taxonDegs'

prfx = ecofile.split('_')[0]
file = pd.read_csv(filename, header=0).sort_values(taxa).reset_index(drop=True)
ecoDegs = pd.read_csv(ecofile, header=0).sort_values(ecotaxa).reset_index(drop=True)

mtchs = []
for col in ecoDegs.columns[1:]:
    if len(mtchs) % 100 == 0: print(f"iter {len(mtchs)}")
    tmp = pd.concat([file[arbDegs], ecoDegs[col]], axis=1)
    tmp2 = tmp.groupby([arbDegs, col]).size().reset_index().rename(columns={0:'cnt'}).sort_values('cnt', ascending=False)
    done, dfr = [], pd.DataFrame()
    for sp in tmp2[arbDegs].unique(): dfr = pd.concat([dfr, tmp2[tmp2[arbDegs] == sp].iloc[0,:]], axis=1)
    dfr = dfr.T
    unqTaxa = dfr[col].unique()
    dupes = dfr[dfr[col].duplicated()]
    if not dupes.empty:
        for i in dupes.index:
            try:
                qw = tmp2[tmp2[arbDegs] == dfr[arbDegs][i]][~tmp2[tmp2[arbDegs] == dfr[arbDegs][i]][col].isin(unqTaxa)].iloc[0,:]
                dfr.loc[i] = qw
                if tmp2[tmp2[arbDegs] == dfr[arbDegs][i]][~tmp2[tmp2[arbDegs] == dfr[arbDegs][i]][col].isin(unqTaxa)].shape[0] > 1: print(tmp2[tmp2[arbDegs] == dfr[arbDegs][i]][~tmp2[tmp2[arbDegs] == dfr[arbDegs][i]][col].isin(unqTaxa)])
            except:
                dfr = dfr.drop([i])
            unqTaxa = dfr[col].unique() #DO I NEED???

    mtchs.append(dfr.cnt.sum())

g = open(f'{prfx}_sppDelimMatches.csv','w')
cnt = 0
g.write("ecoDegsCol,match\n")
for col in ecoDegs.columns[1:]:
    g.write(f"{col},{mtchs[cnt]}\n")
    cnt += 1
g.close()

print(f'Outputted: "{prfx}_sppDelimMatches.csv"')