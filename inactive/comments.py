import numpy as np

performance_star_ids = np.loadtxt("./healpix_05899/performance.txt", dtype=np.int64, delimiter=' ', usecols=1)
performance_data = np.loadtxt("./healpix_05899/performance.txt", dtype=np.float64, delimiter=' ', usecols=(3, 4, 5))


class Star(object):

    def __init__(self, star_id):
        self.filename = "./healpix_05899/star_0"+str(star_id)+".txt"
        self.main_data = np.loadtxt(self.filename, dtype=np.float64, delimiter=' ', usecols=(0, 1, 2))
        self.row_index = np.where(performance_star_ids == star_id)[0][0]
        self.APASS_gmag = performance_data[self.row_index][0]
        self.RA = performance_data[self.row_index][1]
        self.DEC = performance_data[self.row_index][2]

        sort = self.main_data[:, 0].argsort()
        self.mjd = self.main_data[:, 0][sort]
        self.Evry_mag = self.main_data[:, 1][sort]
        self.err = self.main_data[:, 2][sort]

    def find_delta_t(self):
        sorted = np.sort(self.mjd)
        print (sorted[-1]-sorted[0]), " day range of data"
        delta_t_list = []
        for i in range(len(sorted)-1):
            diff = np.absolute(sorted[i] - sorted[i+1])
            delta_t_list.append(diff)
        return np.array(delta_t_list)


def mk_timeseries(star):
    t_baseline = 2.0/(60.*24.)  # 2-minute sampling baseline, expressed in MJD units
    exp_baseline = 2.0/(60.*24.)  # exposure time baseline, MJD units (for right now = same as t_baseline)
    duration = np.max(star.mjd) - np.min(star.mjd)  # total duration of the time series
    
    float_bins = duration/t_baseline  # number of time elements needed to span duration (non-integer)
    power_round = np.round(np.log10(float_bins)/np.log10(2)+0.5).astype(int)  # float_bins expressed as a power of 2, then rounded up to next integer
    bins = 2**power_round  # length of time series to cover duration and still be a power of 2 (for FFT efficiency)

    print bins
    t0 = np.min(star.mjd)                   # start time of the time series
    t = (np.arange(bins)*t_baseline)+t0     # start times of each resampled bin
    rebinned_flux = np.zeros(bins)          # flux time series; initially all zeros
    Evry_flux = 10.0**(-0.4*star.Evry_mag)    # convert mags to flux
    indices = np.arange(bins)

    # now, we want to re-sample the MJD, flux pairs onto the evenly-sampled time series
    # tct=0

    for i in range(len(star.mjd)):  # loop over the input MJD, flux pairs
        if i == 1000*np.round(i/1000.0):
            print i, len(star.mjd)
        j = np.where(t <= star.mjd[i])[0][-1]  # Closest resampled bin start time to the left of original time stamp
        k = np.where(t <= (star.mjd[i]+exp_baseline))[0][-1]  # Closest resampled bin start time to the left of original time stamp plus the baseline exposure time

        # 3 cases are possible
        # CASE 1 -> all the flux is in one bin
        if j == k:
            rebinned_flux[j] += Evry_flux[i]

        # CASE 2 -> flux is in two time bins
        if k == (j+1):
            frac1 = (t[j]+t_baseline-star.mjd[i])/exp_baseline  # fraction of the exposure in first bin
            rebinned_flux[j] = rebinned_flux[j]+frac1*Evry_flux[i]  # add to the first bin
            frac2 = 1.0-frac1  # rest of the flux
            rebinned_flux[k] = rebinned_flux[k]+frac2*Evry_flux[i]  # add to the second bin

        # CASE 3 -> flux in three time bins
        if k == (j+2):
            frac1 = (t[j]+t_baseline-star.mjd[i])/exp_baseline  # fraction of the exposure in first bin
            rebinned_flux[j] = rebinned_flux[j]+frac1*Evry_flux[i]  # add to the first bin
            frac2 = (t_baseline/exp_baseline)  # fraction of the exposure in the second bin
            rebinned_flux[j+1] = rebinned_flux[j+1]+frac2*Evry_flux[i]  # add to the second bin
            frac3 = 1.0-(frac1+frac2)  # fraction of the exposure in the third bin
            rebinned_flux[k] = rebinned_flux[k]+frac3*Evry_flux[i]  # add to the third bin

    return t, rebinned_flux

