import pickle
import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy
from astropy.utils.console import ProgressBar

starbins = pickle.load(open('starbins.p', 'rb'))
scalefactors = pickle.load(open('scalefactors.p', 'rb'))
scalebins = deepcopy(starbins)

for i in range(len(scalebins)):
    for star in scalebins[i]:
        star.data['evry_err'] = star.data['evry_err']*np.sqrt(scalefactors[i])

print('Plotting and saving figures...')
pb = ProgressBar(len(starbins))

for j in range(len(starbins)):
    star = starbins[j][-1]
    scale = scalebins[j][-1]
    mag = np.mean(star.data['evry_mag'])
    scalefactor = scalefactors[j]
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=[10, 10])
    ax1.errorbar(star.data['mjd'], star.data['evry_mag'], yerr=star.data['evry_err'], fmt='ko', ms=2, elinewidth=.5, alpha=0.25)
    ax2.errorbar(scale.data['mjd'], scale.data['evry_mag'], yerr=scale.data['evry_err'], fmt='ko', ms=2, elinewidth=.5, alpha=0.25)
    ax1.set_title('Star with Unscaled Error, EvryMag = {:.2f}'.format(mag))
    ax2.set_title('Same Star with Scaled Error, Scale Factor = {:.4f}'.format(scalefactor))
    ax1.set_xlabel('MJD')
    ax2.set_xlabel('MJD')
    ax1.set_ylabel('EvryMag')
    ax2.set_ylabel('EvryMag')
    plt.subplots_adjust(top=0.961, bottom=0.065, left=0.084, right=0.982, hspace=0.218, wspace=0.2)
#    plt.savefig('./scalefactor_figs/mag{:.2f}_scale{:.4f}.pdf'.format(mag, scalefactor), dpi=300)
    plt.show()
    plt.close()
    pb.update()
