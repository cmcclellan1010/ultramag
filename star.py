import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord as sc
from astropy.table import Table
from matplotlib import pyplot as plt
from scipy.stats import chisquare
import os
import warnings

def dev():
    star = Star('/media/connor/Lore/ultramag/steve_5899_v2/star_058764897')
    return star

class Star(object):
    '''
    A class to handle individual stars from Evryscope data.
    '''

    def __init__(self, starfile):
        self.filename = os.path.abspath(starfile)
        self.datapath = os.path.dirname(os.path.abspath(starfile))+'/'
        self.id = int(starfile.split('star_0')[1].split('.txt')[0])
        self.filtered = False
        self.get_data()


    def get_data(self):

        performance_star_ids = np.loadtxt(self.datapath+"performance.txt", dtype=np.int64, delimiter=' ', usecols=1)
        performance_data = np.loadtxt(self.datapath+"performance.txt", dtype=np.float64, delimiter=' ', usecols=(3, 4, 5))

        self.data = np.loadtxt(self.filename, dtype=np.float64, delimiter=' ', usecols=(0, 1, 2))
        try:
            row_index = np.where(performance_star_ids == self.id)[0][0]
        except IndexError:
            warnings.warn('No performance data found for this star. Skipping.')
            self.data = None
            return
        self.apass_gmag = performance_data[row_index][0]
        self.ra = performance_data[row_index][1]
        self.dec = performance_data[row_index][2]

        s = sc(ra=self.ra*u.degree, dec=self.dec*u.degree).to_string('hmsdms')
        self.ra_hms = s[0:2]+':'+s[3:5]+':'+s[6:13]
        self.dec_dms = s[15:18]+':'+s[19:21]+':'+s[22:29]

        sort = self.data[:, 0].argsort()
        mjd = self.data[:, 0][sort]
        mag = self.data[:, 1][sort]
        err = self.data[:, 2][sort]
        self.data = Table([mjd, mag, err], names=['mjd', 'mag', 'err'], masked=True)


    def find_delta_t(self):
        sorted = np.sort(self.data['mjd'])
        print (sorted[-1]-sorted[0]), " day range of data"
        delta_t_list = []
        for i in range(len(sorted)-1):
            diff = np.absolute(sorted[i] - sorted[i+1])
            delta_t_list.append(diff)
        return np.array(delta_t_list)


    def filter(self, tolerance=0.1):
        if self.filtered == True:
            return
        else:
            self.data.remove_rows(self.data['err'] > 1)
            self.data.remove_rows(self.data['err'] < 0)
            self.data.sort('err')
            self.data.reverse()
            n_to_kill = int(np.round(0.1*len(self.data)))
            self.data.remove_rows(range(n_to_kill))
            self.data.sort('mjd')
            self.filtered = True


    def rebin(self):
  
        if self.filtered is False:
            self.filter()
        
        t_baseline = 2.0/(60.*24.)
        exp_baseline = 2.0/(60.*24.)
        duration = np.max(self.data['mjd']) - np.min(self.data['mjd'])

        float_bins = duration/t_baseline
        power_round = np.round(np.log10(float_bins)/np.log10(2)+0.5).astype(int)
        bins = 2**power_round

        t0 = np.min(self.data['mjd'])
        t = (np.arange(bins)*t_baseline)+t0
        rebinned_flux = np.zeros(bins)
        rebinned_err = np.zeros(bins)
        Evry_flux = 10.0**(-0.4*self.data['mag'])
        Evry_flux_err = np.absolute(Evry_flux*self.data['err']*np.log(10)/(-2.5))

        for i in range(len(self.data['mjd'])):
            j = np.where(t <= self.data['mjd'][i])[0][-1]
            k = np.where(t <= (self.data['mjd'][i]+exp_baseline))[0][-1]

            if j == k:
                rebinned_flux[j] += Evry_flux[i]
                rebinned_err[j] += Evry_flux_err[i]
            if k == (j+1):
                frac1 = (t[j]+t_baseline-self.data['mjd'][i])/exp_baseline
                rebinned_flux[j] += frac1*Evry_flux[i]
                rebinned_err[j] += frac1*Evry_flux_err[i]
                frac2 = 1.0-frac1
                rebinned_flux[k] += frac2*Evry_flux[i]
                rebinned_err[k] += frac2*Evry_flux_err[i]

            if k == (j+2):
                frac1 = (t[j]+t_baseline-self.data['mjd'][i])/exp_baseline
                rebinned_flux[j] += frac1*Evry_flux[i]
                rebinned_err[j] += frac1*Evry_flux_err[i]
                frac2 = (t_baseline/exp_baseline)
                rebinned_flux[j+1] += frac2*Evry_flux[i]
                rebinned_err[j+1] += frac2*Evry_flux_err[i]
                frac3 = 1.0-(frac1+frac2)
                rebinned_flux[k] += frac3*Evry_flux[i]
                rebinned_err[k] += frac3*Evry_flux_err[i]

        t -= t0
        self.flux = rebinned_flux
        self.err = rebinned_err
        return rebinned_flux, rebinned_err

#    def chisqnu_flux(self, plot=True):
#        
#        try:
#            flux_array = self.flux
#            err_array = self.err
#        except AttributeError:
#            flux_array, err_array = self.rebin()
#        
#        nonzero = np.where(flux_array != 0.)
#        flux_array = flux_array[nonzero] * 1e7
#        err_array = err_array[nonzero] * 1e7
#        x = np.linspace(0, len(flux_array) - 1, len(flux_array))

#        z = np.polyfit(x, flux_array, 1)
#        chsq = np.sum(((flux_array - np.polyval(z, x)) ** 2.)/err_array**2.)/(len(flux_array)-1)
#        p = np.poly1d(z)

#        if plot:
#            fig = plt.figure()
#            ax = fig.add_subplot(111)
#            ax.errorbar(x, flux_array, yerr=err_array, ms=2, fmt='ko', elinewidth=.5, alpha=0.5)
#            ax.plot(x, p(x), 'b--')
#            ax.set_xlabel('Index')
#            ax.set_ylabel('Flux (1e-7)')
#            ax.text(0.1, 0.9, "Reduced Chi Squared: %.4f" % chsq, transform=ax.transAxes)
#            plt.show()

#        return chsq

    def chsq(self, plot=False):
        if self.filtered is False:
            self.filter()

        x = np.linspace(0, len(self.data) - 1, len(self.data))

        z = np.polyfit(x, self.data['mag'], 1)
        chsq = np.sum(((self.data['mag'] - np.polyval(z, x)) ** 2.)/self.data['err']**2.)/(len(self.data)-1)
        p = np.poly1d(z)

        if plot:
            fig = plt.figure(figsize=[10,4])
            ax = fig.add_subplot(111)
            ax.errorbar(x, self.data['mag'], yerr=self.data['err'], ms=2, fmt='ko', elinewidth=.5, alpha=0.5)
            ax.plot(x, p(x), 'b--')
            ax.set_xlabel('Index')
            ax.set_ylabel('Evryscope Magnitude')
            ax.text(0.1, 0.9, "Reduced Chi Squared: {:.4f}".format(chsq), transform=ax.transAxes)
            plt.title('Magnitude = {:.2f} Mag'.format(np.mean(self.data['mag'])))
            plt.show()

        return chsq

    def split_nights(self):
        """
        Split a star's time series data by night of observation.

        Parameters
        ----------
        None
        
        Returns
        -------
        List of astropy.table.Table objects
        """
        mjds = np.array(list(self.data['mjd']))
        gaps = mjds[1:] - mjds[:-1] # Time between consecutive data points
        gap_indices = np.where(gaps > 0.01)[0]

        nights = [self.data[:gap_indices[0]]]
        for i in range(1, len(gap_indices)-1):
            nights.append(self.data[gap_indices[i]+1:gap_indices[i+1]])
        nights.append(self.data[(gap_indices[-1]+1):])
        return nights

    def export(self, path='', rescale=True):
        if rescale:
            
            try:
                data = self.flux
            except AttributeError:
                self.filter()
                data = self.rebin()[0]

            data = data*10e6
            zero_loc = np.where(data == 0.)
            nonzero_loc = np.where(data != 0.)
            median = np.median(data[nonzero_loc])
            data[zero_loc] = median
            data.astype(np.float32).tofile(path+str(self.id)+'.dat')
        else:
            self.rebin()[0].astype(np.float32).tofile(path + str(self.id) + '.dat')

        descriptors = np.array([' Data file name without suffix          =  ',
                                ' Telescope used                         =  ',
                                ' Instrument used                        =  ',
                                ' Object being observed                  =  ',
                                ' J2000 Right Ascension (hh:mm:ss.ssss)  =  ',
                                ' J2000 Declination     (dd:mm:ss.ssss)  =  ',
                                ' Data observed by                       =  ',
                                ' Epoch of observation (MJD)             =  ',
                                ' Barycentered?           (1 yes, 0 no)  =  ',
                                ' Number of bins in the time series      =  ',
                                ' Width of each time series bin (sec)    =  ',
                                ' Any breaks in the data? (1 yes, 0 no)  =  ',
                                ' Type of observation (EM band)          =  ',
                                ' Photometric filter used                =  ',
                                ' Field-of-view diameter (arcsec)        =  ',
                                ' Central wavelength (nm)                =  ',
                                ' Bandpass (nm)                          =  ',
                                ' Data analyzed by                       =  ',
                                ' Any additional notes:',
                                '   none'])

        values = np.array([str(self.id), 'Evryscope', 'unset', str(self.id), self.ra_hms, self.dec_dms, 'unset',
                           str(np.min(self.data['mjd'])), '0', '524288', '120', '0', 'Optical', 'Other', '180.00', '500.0',
                           '400.0', 'unset', '', ''], dtype=str)

        inf = np.core.defchararray.add(descriptors, values)
        np.savetxt(path+str(self.id)+'.inf', inf, fmt='%s')

if __name__ == '__main__':
    pass

