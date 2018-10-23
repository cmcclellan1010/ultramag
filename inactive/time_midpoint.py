from star import Star
from glob import glob
import numpy as np
import matplotlib.pyplot as plt
import argparse
from astropy.utils.console import ProgressBar
import warnings

warnings.filterwarnings('ignore')


parser = argparse.ArgumentParser()
parser.add_argument('nstars', metavar='nstars', type=int, help='number of stars to process')
parser.add_argument('-s', '--scale', action='store_true')
parser.add_argument('-n', '--normalize', action='store_true')

args = parser.parse_args()
nstars = int(args.nstars)
scale = bool(args.scale)
normalize = bool(args.normalize)

colors = plt.rcParams['axes.prop_cycle'].by_key()['color'][:7]

filenames = glob('/media/connor/Lore/ultramag/steve_5899_v2/star_0*.txt')
stars = []

print('\nLoading stars...')
pb = ProgressBar(len(filenames[:nstars]))
for i, filename in enumerate(filenames[:nstars]):
    star = Star(filename)
    if star.data is not None:
        star.filter()
        stars.append(star)
    pb.update()

fig, ax1 = plt.subplots(1, 1, figsize=[10,10])

print('\nSplitting data into nights and plotting...')
pb = ProgressBar(len(stars))

for i, star in enumerate(stars):
    nights = star.split_nights()
    c = colors[i % len(colors)]

    for night in nights:
        try:
            time_midpoint = 0.5*(night['mjd'][0] + night['mjd'][-1])
            if normalize:
                mags = night['evry_mag']/np.mean(night['evry_mag'])
                errs = None
            else:
                mags = night['evry_mag']
                errs = night['evry_err']
            if scale:
                times = np.abs(night['mjd'] - time_midpoint)/np.max(np.abs(night['mjd'] - time_midpoint))
            else:
                times = np.abs(night['mjd'] - time_midpoint)
            
            ax1.errorbar(times, mags, yerr=errs, fmt='o', ecolor=c, mec=c, mfc=c, ms=3, capsize=.5, elinewidth=0.5, alpha=0.5)
        except IndexError:
            continue
    pb.update()

ax1.set_title('Magnitude v. Time of Night For {} Star(s)'.format(len(stars)))
if scale:
    ax1.set_xlabel('Time Difference to Middle of Night (|MJD - Midpoint|) (Scaled)')
else:
    ax1.set_xlabel('Time Difference to Middle of Night In Days (|MJD - Midpoint|)')
if normalize:
    ax1.set_ylabel('Normalized Evryscope Magnitude')
else:
    ax1.set_ylabel('Evryscope Magnitude')
plt.tight_layout()
plt.show()


