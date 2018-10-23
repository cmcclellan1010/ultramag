
from glob import glob
from star import Star
import matplotlib.pyplot as plt
  
star = Star('/media/connor/Lore/ultramag/steve_5899_v2/star_059185383.txt')

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=[10, 10])
ax1.errorbar(star.data['mjd'], star.data['evry_mag'], yerr=star.data['evry_err'], fmt='o', ms=2)
star.filter()
ax2.errorbar(star.data['mjd'], star.data['evry_mag'], yerr=star.data['evry_err'], fmt='o', ms=2)
ax1.set_title('Unfiltered Magnitude v. MJD')
ax2.set_title('Flags and 10% Highest Error Thrown Out')
ax1.set_xlabel('MJD')
ax1.set_ylabel('EvryScope Magnitude')
ax2.set_xlabel('MJD')
ax2.set_ylabel('EvryScope Magnitude')
plt.subplots_adjust(top=0.926,
bottom=0.121,
left=0.17,
right=0.971,
hspace=0.511,
wspace=0.2)
plt.show()

