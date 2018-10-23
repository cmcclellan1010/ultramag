import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
import numpy as np
import os
from func import Star, performance_star_ids, performance_data

# cwd = "/home/connor/programming/PycharmProjects/presto_projects/"
# os.chdir(cwd)

# SAVE DATA -- DO ONLY ONCE
# mag = np.zeros((7155, len(performance_star_ids)))
# err = np.zeros((7155, len(performance_star_ids)))
# for i, star_id in enumerate(performance_star_ids):
#     star = Star(star_id)
#     mag[:, i][:len(star.Evry_mag)] = star.Evry_mag
#     err[:, i][:len(star.err)] = star.err
#
# np.save('mag', mag)
# np.save('err', err)

mag = np.load('mag.npy')
err = np.load('err.npy')

mag_avgs = np.zeros(np.shape(mag)[1])
for i in range(np.shape(mag)[1]):
    index = np.where(mag[:, i] != 0)
    mag_avgs[i] = np.average(mag[:, i][index])
# print mag_avgs


# Check mins and maxes of arrays
# mags = []
# maxes = []
# for i, star_id in enumerate(performance_star_ids):
#     star = Star(star_id)
#     mags.append(np.average(star.Evry_mag))
#     maxes.append(np.max(star.Evry_mag))
# print mags
# print max(mags)
# print min(mags)
# print np.average(mags)
#
# print maxes
# print max(maxes)
# print min(maxes)
# print np.average(maxes)


index_13to14 = np.where(np.logical_and(13 < mag_avgs, mag_avgs < 14))[0]
index_12to13 = np.where(np.logical_and(12 < mag_avgs, mag_avgs < 13))[0]

# mag_13to14 = mag[:, index_13to14][np.where(mag[:, index_13to14] != 0)]
# mag_12to13 = mag[:, index_12to13][np.where(mag[:, index_12to13] != 0)]

err_13to14 = err[:, index_13to14][np.where(err[:, index_13to14] != 0)]
err_12to13 = err[:, index_12to13][np.where(err[:, index_12to13] != 0)]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(17,11))

ax1.hist(err_12to13, 100)
ax2.hist(err_13to14, 100)
plt.suptitle("Error in Evryscope Magnitudes Of Stars In Healpix 05899")
ax1.set_title("For Stars With Average Magnitudes Between 12 and 13")
ax1.set_ylabel("N")
ax1.set_xlabel("Evryscope Mag Err")
ax2.set_ylabel("N")
ax2.set_xlabel("Evryscope Mag Err")
ax2.set_title("For Stars With Average Magnitudes Between 13 and 14")
ax1.grid(True)
ax2.grid(True)
plt.show()

flat_data = np.ndarray.flatten(mag[np.where(mag != 0)])
flat_err = np.ndarray.flatten(err[np.where(err != 0)])

fig, (ax3, ax4) = plt.subplots(1, 2, figsize=(17, 11))
n, bin_locations, patches = ax3.hist(flat_data, 150)
ax4.hist(flat_err, 150)
bincenters = 0.5*(bin_locations[1:]+bin_locations[:-1])

g_init = models.Gaussian1D(amplitude=100000, mean=np.mean(flat_data), stddev=np.std(flat_data))
fit_g = fitting.LevMarLSQFitter()
g = fit_g(g_init, bincenters, n)

mean1, stddev1 = g.mean[0], g.stddev[0]
print "Mean: ", mean1
print "Stddev: ", stddev1

# ax3.plot(bincenters, g(bincenters), 'r--', linewidth=2)
ax3.set_xlabel('Evryscope Magnitude')
ax3.set_ylabel('N')
ax3.set_xlim([8.5, 14.5])
ax3.grid(True)
plt.suptitle("Histogram Of All Datapoints For Stars In Healpix 05899")
ax3.set_title('Evryscope Magnitudes')
ax3.text(9, 67000, "Mean: %.2f mag" % mean1)
ax3.text(9, 65000, "Stddev: %.2f mag" % stddev1)
ax4.set_title("Error In Evryscope Magnitudes")
ax4.set_xlabel("Error In Evryscope Magnitudes")
ax4.set_ylabel('N')
ax4.grid(True)
plt.show()
print mean1, stddev1