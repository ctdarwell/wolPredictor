#build CSV of eco-based taxon divisions
#for numpy implementation maybe able to use np.meshgrid for combinations: 
#https://stackoverflow.com/questions/1208118/using-numpy-to-build-an-array-of-all-combinations-of-two-arrays

#THIS NEEDS TO OUTPUT AS INT AND NOT FLOATS!
import numpy as np, sys, getopt
import pandas as pd

#Load data, define parameters
swch = 1 #[1] constrains pops to numerically adjacent values; [0] allows disjunct values as a pop: e.g. [100,200,300], option[1] does not allow [100,300] as a pop 
ecoVar = 'elevation'
community = 'community'
filename = 'exaData_faoReview.csv'
prfx = 'eco'

#check altered input params
try: #input sys args
    myopts, args = getopt.getopt(sys.argv[1:],"s:e:c:d:p:")
except:
    print('I/O ERROR!')
    sys.exit(2)

#Sort non-default params and check for input errors
print('\n### Check inputted params ###')
for a, b in myopts:
    print(a, b)
    if a == '-s': swch = int(b)
    if a == '-e': ecoVar = b
    if a == '-c': community = b
    if a == '-d': filename = b
    if a == '-p': prfx = b

#housekeeping
dat = pd.read_csv(filename, header=0) # default = datMLformatted_geiger2.csv
taxa = dat.taxa.sort_values() #sorted(list(set(dat['taxa'])))
comms = np.sort(dat[community].unique()).tolist() #sorted(list(set(dat[community])))
sfx = ""
if swch == 1: sfx = '_constr'
cladeDivNames = f"{prfx}Delim_{community}X{ecoVar}{sfx}" #add a prefix to file names - do NOT use 'x'

def main():
    taxonDesignations = pd.DataFrame(taxa, columns = ['taxa'])
    
    elev_dik = {}
    for comm in comms: elev_dik.update({comm :getCombs(dat[ecoVar][dat.community == comm].unique().tolist(), swch)})
    clade_divs = cladeDivyer(elev_dik)
    taxonDesignations = taxonDesignator(taxonDesignations, clade_divs)
    taxonDesignations.iloc[:, 1:] = taxonDesignations.iloc[:, 1:].astype(int)
    taxonDesignations.to_csv(f"{cladeDivNames}.csv", index = False)

def taxonDesignator(taxonDesignations, clade_divs):
    print(f"nPermuations: {len(clade_divs)}")
    cntr = 0
    for div in clade_divs: #eg [[[2700, 700]], [[1200, 100, 2700]]]
        if cntr % 100 == 0: print(f"iter {cntr}")
        cntr2, cntr3 = 0, 0
        taxonDesignations=taxonDesignations.assign(A = np.nan)
        taxonDesignations.rename(columns={list(taxonDesignations)[-1]: 'div' + str(cntr).zfill(6)}, inplace=True)
        for part in div: #eg [[2700, 700]...or...[[1200, 100, 2700]]
            qw=str(part).split('], [')# eg ['[[2700, 700]]']...or....['[[100, 2700', '1200]]']
            for q in qw:
                clade = q.replace('[','').replace(']','').split(',')
                for elevation in clade:
                    indvs = list(dat['taxa'][dat[ecoVar] == int(elevation)][dat[community] == comms[cntr3]])
                    for indv in indvs: taxonDesignations.loc[taxonDesignations['taxa'] == indv, 'div' + str(cntr).zfill(6)] = cntr2
                cntr2 += 1
            cntr3 += 1
        cntr += 1
    
    return taxonDesignations

def cladeDivyer(elev_dik):
    tots, tots2, posns, selns, clade_divs = [], [], [0] * len(comms), [], []
    for comm in comms: tots.append(len(elev_dik.get(comm))), tots2.append(len(elev_dik.get(comm)) - 1) #build counting system lists
    while tots != posns:
        for comm in comms: selns.append(elev_dik.get(comm)[posns[comms.index(comm)]]) #COMMENT??
        clade_divs.append(selns)
        selns = []
        if posns == tots2: break #can maybe rm & replace while loop w 'tots2 != posns'
        posns[-1] += 1
        for indx in range(-1, (len(posns) * -1) - 1, -1):
            if posns[indx] == tots[indx]:
                posns[indx - 1] =  posns[indx - 1] + 1
                for i in range(indx, 0): posns[i] = 0
            else:
                break

    return clade_divs

def getCombs(elevs, swch):
    #out all combinations of this set
    elevs = sorted(elevs)
    d = {item: idx for idx, item in enumerate(elevs)}
    
    combs = []
    for n, parts in enumerate(partition(elevs), 1): #https://stackoverflow.com/questions/19368375/set-partitions-in-python
        s = 0
        for part in sorted(parts):
            tmp = []
            [tmp.append(d.get(item)) for item in part]
            tmp = sorted(set(tmp))
            for t in tmp[:-1]:
                if tmp[tmp.index(t)+1] - t > 1: s = swch
                
        if s == 0: combs.append(parts)
    
    return combs

def partition(collection):
    if len(collection) == 1:
        yield [ collection ]
        return

    first = collection[0]
    for smaller in partition(collection[1:]):
        # insert `first` in each of the subpartition's subsets
        for n, subset in enumerate(smaller):
            yield smaller[:n] + [[ first ] + subset]  + smaller[n+1:]
        # put `first` in its own subset 
        yield [ [ first ] ] + smaller

if __name__ == '__main__':
    main()
