import pandas as pd
import sys
import matplotlib.pyplot as plt

def makePDFv2(f):
    '''
    calculate predictability without considering 'noWol'
    '''
    print("Figure title:", f)
    props=[]
    rng = f.split('incr')[1].split('x')[0]
    low = int(rng.split('-')[0])
    up = int(rng.split('-')[1])
    print(low, up)
    incr = 100/float(f.split('x')[1].split('_')[0])
    prg = f.split('prg')[1]
    dat = pd.read_csv(f + ".csv", header=0)
    runs = list(dat)[2:]
    anals, tots, cnts, maximum = [], [], [], 0
    for run in runs:
        if list(dat[run]).count('noWol') + list(dat[run]).count('Vetoed') == len(list(dat[run])): continue
        cnt, cnt2, tot = 0, 0, 0
        for cell in dat[run]:
            if cell == dat.loc[cnt][1] and cell != 'noWol': tot += 1
            if dat.loc[cnt][1] != 'noWol': cnt2 += 1
            cnt += 1
        if cnt2 == 0: cnt2 = 1000000
        props.append((tot/cnt2)*100)
        anals.append(run), tots.append(tot), cnts.append(cnt2)
        if tot > maximum:
            print('Best run: ', run, tot, cnt2)
            maximum = tot
            bestThr = run[6:]

    print("maximum", maximum/cnt2*100,"%")
    maximum = str(maximum/cnt2*100)[:5]
    fig1 = plt.figure()
    plt.title('incr@'+str(incr)+'%; purge '+prg+'%; max '+maximum+'% @0.'+bestThr)
    plt.xticks([], [])
    plt.xlabel('Species delim. range (%): ' + str(low) +  ' to ' + str(up) + ' - NB non-linear x-axis')
    plt.plot(props)
    fig1.savefig(f + '_posStrains.pdf', dpi = 900)

    A = pd.DataFrame(cnts, columns = list("A"))
    B = pd.DataFrame(tots, columns = list("B"))
    C = pd.DataFrame(anals, columns = list("C"))

    tab = pd.concat([C, A, B], axis = 1)
    tab.columns = ['anals', 'totalPreds', 'correctCnts']
    tab.to_csv(f + '_posStrains_Summary.csv', index = False)

f = sys.argv[1]
makePDFv2(f)