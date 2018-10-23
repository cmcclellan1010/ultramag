from func import *
from matplotlib import pyplot as plt


def spreadplot(id):
    star = Star(id)

    plt.figure(figsize=(12, 12))
    plt.errorbar(star.mjd, star.Evry_mag, yerr=star.err, fmt='o', ms=3)
    plt.xlabel('MJD')
    plt.ylabel("Evryscope Magnitude")
    plt.title("Star "+str(bright_id)+" (Avg mag: {})".format(np.average(star.Evry_mag)))
    plt.show()


bright_id = performance_star_ids[0]
spreadplot(bright_id)

# other = Star(58764897)
#
# plt.figure(figsize=(12,12))
# plt.errorbar(other.mjd, other.Evry_mag, yerr=other.err, fmt='o', ms=3)
# plt.xlabel('MJD')
# plt.ylabel("Evryscope Magnitude")
# plt.title("Star "+str(other.id))
# plt.show()
