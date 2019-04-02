import pickle
from copy import deepcopy
import numpy as np
from astropy.utils.console import ProgressBar
import matplotlib.pyplot as plt
from astropy.table import Table
import os.path

starbins = pickle.load(open('/media/connor/Lore/ultramag/starbins.p', 'rb'))
scalefactors = pickle.load(open('/media/connor/Lore/ultramag/scalefactors.p', 'rb'))
scalebins = deepcopy(starbins)

for i in range(len(scalebins)):
    for star in scalebins[i]:
        star.data['evry_err'] = star.data['evry_err']*np.sqrt(scalefactors[i])
        star.filtered = True

scalestars = np.concatenate(scalebins).ravel()

if os.path.isfile('/media/connor/Lore/ultramag/PRESTO/id_chsq_table.dat'):
    chsqs = np.array(list(Table.read('/media/connor/Lore/ultramag/PRESTO/id_chsq_table.dat', format='ascii')['CHISQ']))
    scalestars = np.flip(scalestars[np.argsort(chsqs)], axis=0)
    starids = [scalestars[i].id for i in range(len(scalestars))]
else:
    print('Calculating uncertainty-scaled chsq...') 
    pb = ProgressBar(len(scalestars))

    chsqs = []
    for star in scalestars:
        chsqs.append(star.chsq(plot=False))
        pb.update()
    
    scalestars = np.flip(scalestars[np.argsort(chsqs)], axis=0)
    starids = [scalestars[i].id for i in range(len(scalestars))]
    chsqs = np.flip(np.array(sorted(chsqs)), axis=0)

    t = Table([starids, chsqs], names=['STAR_ID', 'CHISQ'])
    t.write('/media/connor/Lore/ultramag/PRESTO/id_chsq_table.dat', format='ascii', overwrite=True)

# Now we have all the stars properly scaled and sorted by decreasing 
# variability. Time to prepare the time series data.

print('\nExporting stars...')
pb = ProgressBar(len(scalestars))
for i in range(len(scalestars)):
    star = scalestars[i]
    if os.path.isfile('/media/connor/Lore/ultramag/PRESTO/'+str(star.id)+'.dat'):
        pb.update()    
    else:
        star.export(path='/media/connor/Lore/ultramag/PRESTO/')
        pb.update()
