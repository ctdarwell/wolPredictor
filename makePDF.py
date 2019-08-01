import numpy as np
import matplotlib.pyplot as plt
import csv
from numpy import genfromtxt

def makePDF(f, gap):
    with open(f) as file:
        reader = csv.reader(file)
        columns = next(reader)
        colmap = dict(zip(range(len(columns)), columns))
    file = None

    dat = genfromtxt(f, delimiter=',', dtype = 'U20', skip_header = 1) 
    wsp = np.unique(dat[:, 1])
    noWol = np.where(wsp == 'noWol')[0][0] #id posn of noWol

    neg, pos, posRes, negRes = [], [], [], []
    for j in range(2, len(colmap)):
        if np.shape(np.where(dat[:, j] == 'noWol'))[1] == np.shape(dat)[0]:
            neg.append(0), pos.append(0)
            continue
        tmp = []
        for w in wsp:
            qw = np.where(dat[:, 1] == w)
            qw2 = np.where(dat[:, j] == w)
            tmp.append(len(np.intersect1d(qw, qw2)))
        neg.append(tmp[noWol] * -1), pos.append(sum(tmp[:noWol]) + sum(tmp[noWol + 1:]))
        negRes.append(tmp[noWol] * -1), posRes.append(sum(tmp[:noWol]) + sum(tmp[noWol + 1:]))

    data = np.array([pos, neg])
    
    data_shape = np.shape(data)

    # Take negative and positive data apart and cumulate
    def get_cumulated_array(data, **kwargs):
        cum = data.clip(**kwargs)
        cum = np.cumsum(cum, axis=0)
        d = np.zeros(np.shape(data))
        d[1:] = cum[:-1]
        return d  

    cumulated_data = get_cumulated_array(data, min=0)
    cumulated_data_neg = get_cumulated_array(data, max=0)

    # Re-merge negative and positive data.
    row_mask = (data<0)
    cumulated_data[row_mask] = cumulated_data_neg[row_mask]
    data_stack = cumulated_data

    #build tick mark lists
    nSpp, sppPosns, nPges = [], [], [] #sppPosns comes up short sometimes
    M = columns[-1].split('_')[0].split('pp')[1]
    [nSpp.append(i) for i in range(gap, int(M) + 1, gap)]
    [nPges.append(1) for c in columns if M in c]
    x = len(pos) - (len(nPges) * (int(M) % gap))
    while x > 0:
        sppPosns.append(x)
        x -= (len(nPges) * gap)

    cols = ['C1', 'b']
    fig = plt.figure()
    ax = plt.subplot(111)
    plt.ylim([-150, 100])
    plt.xlabel('Number of species delimitation clusters')
    plt.ylabel('Number of correct predictions\n(negative/positive infections)')
    plt.xticks(sppPosns[::-1], nSpp) # posns: len(pos) - nCols per sp x 10
    plt.yticks([-150, -100, -50, 0, 50, 100], [150, 100, 50, 0, 50, 100])

    for i in np.arange(0, data_shape[0]):
        ax.bar(np.arange(data_shape[1]), data[i], bottom=data_stack[i], color=cols[i],)

    #plt.show()
    fig.savefig(f.replace('_correctedWolPreds', '')[:-4] + '.pdf', dpi = 900)
    fig.savefig(f.replace('_correctedWolPreds', '')[:-4] + '.png', dpi = 900)
