from func import Star
import os
import numpy as np
from matplotlib import pyplot as plt
from astropy.modeling import models, fitting

cwd = "/home/connor/programming/PycharmProjects/presto_projects/"
os.chdir(cwd)

performance_star_ids = np.loadtxt("./healpix_05899/performance.txt", dtype=np.int64, delimiter=' ', usecols=1)
performance_data = np.loadtxt("./healpix_05899/performance.txt", dtype=np.float64, delimiter=' ', usecols=(3, 4, 5))


def hist_delta_t(data, ax, bins=100, dt_range=[0, 5]):
    # Convert to minutes
    data_array = np.array(data, dtype=np.float64)
    data_array_mins = data_array*24.*60.
    data_list_mins = list(data_array_mins)
    n, bins, patches = ax.hist(data_list_mins, bins, range=dt_range)


# Collect and organize data by star - do only once

hist_data = []
fig1, ax1 = plt.subplots(1, 1, figsize=(12, 12))
fig2, ax2 = plt.subplots(1, 1, figsize=(12, 12))


for k in range(len(performance_star_ids)):
    star = Star(performance_star_ids[k])
    dt = star.find_delta_t()

    # For histogram
    hist_data.extend(list(dt))

    # For time-dt plot
    sumtime = []
    for i in range(len(dt)):
        sum = np.sum(dt[0:i+1])
        sumtime.append(sum)
    ax1.plot(sumtime, dt)

ax1.set_xlabel("Days since T=0")
ax1.set_ylabel("Delta T (days)")
ax1.set_title("Delta T Between Consecutive Data Points v. Time")

ax2.set_xlabel("Delta T (minutes)")
ax2.set_ylabel("Number of Occurrences")
ax2.set_title("Delta T Between Consecutive Data Points in Evryscope Healpix #05899")

print "Mean: ", np.mean(hist_data)*24.*60., "(min)"
print "Median: ", np.median(hist_data)*24.*60., "(min)"
print "Standard deviation: ", np.std(hist_data)*24.*60., "(min)"

hist_delta_t(hist_data, ax2, bins=200, dt_range=[0, 5])
plt.show()
