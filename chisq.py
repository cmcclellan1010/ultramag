from star import *
from glob import glob
from astropy.utils.console import ProgressBar
import math
import pickle

filenames = glob('/media/connor/Lore/ultramag/steve_5899_v2/star_0*.txt')

stars = []
avg_mags = []

print('Loading stars...')
pb = ProgressBar(len(filenames))

for filename in filenames:
    star = Star(filename)
    if star.data is not None:
        stars.append(star)
        avg_mags.append(np.mean(star.data['evry_mag']))
    pb.update()

stars = np.array(stars)[np.argsort(avg_mags)]
# Now that they're sorted, split them into equal numbered bins and find chsq for each in each bin

nbins = 10
bin_width = int(np.ceil(0.5*len(stars)/nbins)*2)
bin_cutoffs = np.array(range(nbins))*bin_width
bin_cutoffs = np.insert(bin_cutoffs, -1, int(0.5*(bin_cutoffs[-2]+bin_cutoffs[-1])))
bin_cutoffs = np.append(bin_cutoffs, int(0.5*(len(stars)+bin_cutoffs[-1])))
bin_cutoffs = np.append(bin_cutoffs, -1)

bins = []
magnitude_ranges = []
for i in range(len(bin_cutoffs)-1): # 2 bonus bins for bright stars
        bins.append(stars[bin_cutoffs[i]:bin_cutoffs[i+1]])
        some_avg_mags = sorted(avg_mags)[bin_cutoffs[i]:bin_cutoffs[i+1]]
        magnitude_ranges.append((np.min(some_avg_mags), np.max(some_avg_mags)))

print('Binning stars and calculating chsq...')
pb = ProgressBar(len(stars))
chsq_bins = []
for somebin in bins:
    chsq_subbin = []
    for eachstar in somebin:
#        chsq = eachstar.chsq(plot=True)
        chsq = eachstar.chsq(plot=False)
        chsq_subbin.append(chsq)
        pb.update()
    chsq_bins.append(chsq_subbin)

means = [np.mean(subbin) for subbin in chsq_bins]
medians = [np.median(subbin) for subbin in chsq_bins] # MEDIAN IS BETTER

nonvariable_chsqbins = []
nonvariable_starbins = []

# Throw out the 20% most variable stars in each bin (highest chisquared)
for i in range(len(chsq_bins)):

    chsqbin = chsq_bins[i]
    starbin = bins[i]


    n_kept = int(np.ceil(0.8*len(chsqbin)))
    starbin = np.array(starbin)[np.argsort(chsqbin)]
    starbin = starbin[0:n_kept]
    
    chsqbin = sorted(chsqbin)[0:n_kept]

    nonvariable_starbins.append(starbin)
    nonvariable_chsqbins.append(chsqbin)

# Make histograms of the 80% least variable stars, to find a central X^2 that fits most of them

def chsqhist(chsqbin, starbin, plot=True):

    mags = [np.mean(starbin[i].data['evry_mag']) for i in range(len(starbin))]
    maxmag = np.max(mags)
    minmag = np.min(mags)

    mean = np.mean(chsqbin)
    median = np.median(chsqbin)

    print("Mean: {}    Median: {}".format(mean, median))

    if plot:
        plt.figure(figsize=[9, 9])
        plt.hist(chsqbin)
        plt.xlabel('Reduced Chi Squared')
        plt.ylabel('n')
        plt.title('Reduced Chi Squared Histogram for Nonvariable Stars ({:.2f} to {:.2f} Mag Bin)'.format(minmag, maxmag))
        plt.tight_layout()
        plt.show()
    
    return mean, median

means = []
for i in range(len(nonvariable_chsqbins)):
    mean, median = chsqhist(nonvariable_chsqbins[i], nonvariable_starbins[i], plot=True)
    means.append(mean)

pickle.dump(means, open('scalefactors.p', 'wb'))
pickle.dump(bins, open('starbins.p', 'wb'))

# Notes:
# Use mean as scale factor (CHECK)
# May want to use finer bins at the bright end (CHECK)
# split two brightest bins in half to calculate scale factor (CHECK)
