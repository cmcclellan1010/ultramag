import pickle
from copy import deepcopy
import numpy as np
from astropy.utils.console import ProgressBar
import matplotlib.pyplot as plt

starbins = pickle.load(open('starbins.p', 'rb'))
scalefactors = pickle.load(open('scalefactors.p', 'rb'))
scalebins = deepcopy(starbins)

for i in range(len(scalebins)):
    for star in scalebins[i]:
        star.data['evry_err'] = star.data['evry_err']*np.sqrt(scalefactors[i])
        star.filtered = True

scalestars = np.concatenate(scalebins).ravel()

print('Calculating uncertainty-scaled chsq...') 
pb = ProgressBar(len(scalestars))

chsqs = []
for star in scalestars:
    chsqs.append(star.chsq(plot=False))
    pb.update()

# SORT IN DESCENDING ORDER
scalestars = np.flip(scalestars[np.argsort(chsqs)], axis=0)
chsqs = np.flip(np.array(sorted(chsqs)), axis=0)

plt.figure(figsize=(9, 5))
plt.hist(chsqs, bins=np.logspace(np.log10(0.03), np.log10(200), 300))
plt.gca().set_xscale('log')
plt.xlabel('Uncertainty-scaled Chi Squared Nu')
plt.ylabel('n')
plt.title(r'Histogram of Chi Squared Nu, After Multiplying Error by Sqrt(Scale Factor)')
plt.tight_layout()
plt.show()

for i in range(20):
#for i in range(int(len(scalestars)/2), int(len(scalestars)/2)+20):
    star = scalestars[i]
    chsq = chsqs[i]
    mjds = star.data['mjd'] - star.data['mjd'][0]
    fig = plt.figure(figsize=[12,6])
    ax = fig.add_subplot(111)
    ax.errorbar(mjds, star.data['evry_mag'], yerr=star.data['evry_err'], ms=2, fmt='ko', elinewidth=.5, alpha=0.5)
    ax.set_xlabel('Days Since First Observation (Days)')
    ax.set_ylabel('Evryscope Magnitude')
    ax.text(0.1, 0.9, "Reduced Chi Squared: {:.4f}".format(chsq), transform=ax.transAxes)
    plt.title('Mean Magnitude = {:.2f} Mag    Plot Index = {}'.format(np.mean(star.data['evry_mag']), i))
    plt.tight_layout()    
    pickle.dump(ax, open('./interactive_figs/matplotlib_fig_{}.p'.format(i), 'wb'))
#        plt.show()

