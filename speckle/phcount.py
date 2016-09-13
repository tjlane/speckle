
import numpy as np
from scipy.optimize import curve_fit

def multi_gaussian(x, *params):
    y = np.zeros_like(x)
    for i in range(0, len(params), 3):
        ctr = params[i]
        amp = params[i+1]
        wid = params[i+2]
        y = y + amp * np.exp( -((x - ctr)/wid)**2)
    return y


def fit_photon_hist(counts, adu_per_photon=150.0, max_photons=None, plot=False):
    """
    Fit a photon histogram that can be used for photon counting.

    Fits a variable number of Gaussians to an automatically generated 
    photon histogram, and uses that fit to generate bins that can be used
    for photon counting. Parameters good guess for ePix100.

    Parameters
    ----------
    counts : list of ints
        List of droplet counts, in ADUs

    adu_per_photon : float
        A guess of the photon gain on the camera

    plot : bool or str
        If False, do nothing; if True, show a plot of results; if str, save
        a plot of results to that file path

    Returns
    -------
    bins_to_use : list
        A list of the bins to use in photon counting. See np.bincount.
    
    gaussian_parameters : list
        A length 3N list of parameters used in the Gaussian fit. The parameters
        are the center/amplitude/width of each of N Gaussians, in that order.
    """

    counts = np.array(counts)

    # step 1 : generate a histogram
    bins = np.arange(0, counts.max()+1)
    bcnt = np.bincount( np.digitize(counts, bins) )

    # step 2 : curvefit
    if max_photons:
        p_max = max_photons
    else:
        p_max = int(np.ceil(counts.max() / adu_per_photon)) # max photons
    print 'Attempting to fit %d photons...' % p_max

    param_guess = []
    for i in range(p_max):
        param_guess.extend([adu_per_photon*i, 
                            bcnt[int(adu_per_photon)/2:].max()/2.0, 
                            adu_per_photon / 3.0])

    popt, pcov = curve_fit(multi_gaussian, bins, bcnt, p0=param_guess)
    fit = multi_gaussian(bins, *popt)

    # step 3: determine the bins to use
    # (we'll just take the mid points between the peaks)
    centers = np.array([popt[3*i] for i in range(p_max)])
    bins_to_use = (centers[:-1] + centers[1:]) / 2.0


    if plot:

        from matplotlib import pyplot as plt
        plt.plot(bins, bcnt, lw=2)
        plt.plot(bins, fit , 'r-')
        plt.vlines(centers, fit.max()-fit.max()/8., 
                   fit.max(), color='red', lw=2)
        plt.vlines(bins_to_use, 0.0, fit.max())
        plt.xlabel('Droplet ADU')
        plt.ylabel('Observations')

        if type(plot) == str:
            plt.savefig(plot)

        plt.show()


    return bins_to_use, popt


