from star import Star
from glob import glob
import numpy as np
import matplotlib.pyplot as plt
import argparse
from astropy.utils.console import ProgressBar
import warnings

warnings.filterwarnings('ignore')

def rms(x):
    return (np.absolute(np.mean(x**2) - (np.mean(x))**2))**0.5


parser = argparse.ArgumentParser()
parser.add_argument('nstars', metavar='nstars', type=int, help='number of stars to process')
parser.add_argument('nbins', metavar='nbins', type=int, help='number of bins')
parser.add_argument('-s', '--scale', action='store_true')

args = parser.parse_args()
nstars = int(args.nstars)
nbins = int(args.nbins)
scale = bool(args.scale)

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

all_times = []
all_mags = []

for i, star in enumerate(stars):
    nights = star.split_nights()
    c = colors[i % len(colors)]

    for night in nights:
        try:
            time_midpoint = 0.5*(night['mjd'][0] + night['mjd'][-1])
            mags = night['evry_mag']/np.mean(night['evry_mag'])
            if scale:
                times = np.abs(night['mjd'] - time_midpoint)/np.max(np.abs(night['mjd'] - time_midpoint))
            else:
                times = np.abs(night['mjd'] - time_midpoint)

            all_times.extend(list(times))
            all_mags.extend(list(mags))
            
        except IndexError:
            continue
    pb.update()

sorted_times, sorted_mags = zip(*sorted(zip(all_times, all_mags)))
sorted_times, sorted_mags = np.array(sorted_times), np.array(sorted_mags)
mask = ~np.isnan(sorted_times)

sorted_times = sorted_times[mask]
sorted_mags = sorted_mags[mask]

duration = sorted_times[-1] - sorted_times[0]
bin_width = duration/nbins

bin_rmss = []

for i in range(nbins):
    #get data points in this bin
    bin_mags = sorted_mags[np.where(np.logical_and(sorted_times >= bin_width*i, sorted_times <= bin_width*(i+1)))[0]]
    bin_rmss.append(rms(bin_mags))

ax1.plot(bin_width*range(len(bin_rmss)), bin_rmss)
plt.show()

