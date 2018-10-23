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

# Now we have all the stars properly scaled and sorted by decreasing 
# variability. Time to prepare the time series data.

print('\nExporting victims...')
victims = [1, 18, 313]
for i in victims:
    star = scalestars[i]
    #tar.export(path='./PRESTO/')
    print("\nIndex: {}      Star ID: {}".format(i, star.id))
    pb.update()
