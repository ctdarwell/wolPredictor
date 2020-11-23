import numpy as np, math
import matplotlib.pyplot as plt
import pandas as pd

fldr = '.' #output folder
oviLayers = [1, 0.75, .5, 0.25, 0] #survival %s according to oviposition layers in fig - ie all eggs survive in 1st layer nearesr centre as away from parasitoids
fec = 1000 #potential foundress egg load
incrs = 100 #increments in fitness to assess (btwn 0-1)
nreps = 10 #nReplicates to take mean result from
start, jumps = 5, 5 #start value of consp proportions and increments to cycle thru

def fitnesses(x, y, z):
    propCons = x #prop of conspecifics in pop i.e. P(consp mating)
    fitness = [y, z] #fitness due to phenotype - hetero vs consp mating
    
    remain = int(fec * propCons) #viable eggs remaining under CI acounting for heterospp matings

    w_ci = np.full((remain), fitness[1]) #under CI all remaining eggs get consp mating fitness
    wo_ci = np.random.choice(fitness, fec, p=[1 - propCons, propCons]) #wo CI randomly get consp or heterosp fitnesses according to proportions of consp v het matings

    #work out remaining oviposition layer sites if egg load reduced by CI
    divs = fec / len(oviLayers)
    newOviLayers = oviLayers[:int(len(oviLayers) * propCons) +1]
    props, tot = [], 0
    for qw in newOviLayers:
        tot += divs
        if tot < fec * propCons: props.append(divs/(fec * propCons))
        else: props.append(1 - sum(props))

    #calc ovipos site fitnesses - under CI losses of egg load means incr % eggs go to optimal sites
    ovi_w_ci = np.random.choice(newOviLayers, remain, p=props)
    ovi_wo_ci = np.random.choice(oviLayers, fec, p=[1/len(oviLayers)]*len(oviLayers))

    return np.sum(w_ci * ovi_w_ci), np.sum(wo_ci * ovi_wo_ci) #mutiply egg fitnesses from con v het matings X ovi site survuvorship fitnesses and take sum (inclusive/cumulative fitness for foundress)

#Main
vals = []
for prop in range(start, incrs, jumps): #cycle thru at diff %s of conspecifics in pop
    arr = np.zeros((incrs, incrs))
    for hetfit in range(1, incrs + 1): #cycle thru 0-1 of heterosp offspring fitnesses
        for confit in range(1, incrs + 1): #cycle thru 0-1 of consp offspring fitnesses
            tot1, tot2 = 0, 0
            for i in range(nreps): #replication
                wCI, woCI = fitnesses(prop/incrs, hetfit/incrs, confit/incrs)
                tot1 += wCI
                tot2 += woCI
            arr[hetfit - 1, confit - 1] = (tot1/nreps) - (tot2/nreps) #assign mean val for consp - hetsp inclusive fitness at incremental values

    if np.where(arr > 0)[0].size > 0: #if >1 cell has better fitness for CI make figure
        print(f"\nCI wins in {np.where(arr > 0)[0].size/100}% of pixels at {prop}% conspecifics")
        vals.append(np.where(arr > 0)[0].size)
        
        #est colour bar range
        m = arr.min()
        M = arr.max()
        M = max(M, abs(m))
        M = math.ceil(M/50) * 50
        m = M * -1

        fig = plt.figure()
        x = np.linspace(incrs, 0, 6)
        y = np.linspace(0, incrs, 6)
        plt.imshow(arr, cmap='jet',
               vmin=m, vmax=M)
        plt.xticks(x, x/incrs)
        plt.yticks(y, y/incrs)
        plt.title(f"Conspecific pop freq: {prop}%")
        plt.xlabel('Conspecific mating (ω)')
        plt.ylabel('Heterospecific mating (ω)')
        plt.colorbar()
        fig.savefig(f"{fldr}/prop{str(prop).zfill(3)}.png", dpi = 600)

        plt.show()

#summary table
qw = np.array(vals)/100
tmp, x = [], start #see: for prop in range(5, incrs, 5)
for q in qw:
    tmp.append([x,q])
    x += jumps

df = pd.DataFrame(tmp)
df.columns=['% consp','CI favoured (%)']
df.to_csv(f'{fldr}/summary.csv', index=False)

#example at prop conspecifics = 90% with only slight diff in consp v het fitnesses
cnt = 0
for i in range(100):
    x, y = fitnesses(.9, .45, .55)
    if x > y: cnt += 1

print(f"\nCI wins {cnt} out of 100 with hetfit 0.45 and confit 0.55") # around >65% of the time CI gives better fitness
