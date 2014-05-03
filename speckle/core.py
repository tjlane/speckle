
"""
This file contains the main functions for computing quantities of interest for
speckle analysis from CSPAD images.
"""

import os

import numpy as np
from scipy import optimize
from scipy.signal import fftconvolve
from scipy.special import gamma, gammaln
from scipy.special import psi as digamma

try:
    import psana
except ImportError as e:
    print ('Could not import psana, proceeding')


def init_psana(run, expt='', config_file=None):
    """
    Convience function to obtain the datastream and epics objects from psana
    """
    
    if config_file == None:
        config_file = os.path.join("/reg/d/psdm/xcs/%s/res/cfg/%s.cfg" (expt, expt))
    if not os.path.exists(config_file):
        raise IOError('Requested psana configuration file does not exist:'
                      ' %s' % config_file)
    
    psana.setConfigFile(config_file)
    psana.setOption('psana.l3t-accept-only',0)
    print "Loading psana config file:    %s" % config_file
    
    ds = psana.DataSource('exp=%s:run=%d' % (expt, run))
    epics = ds.env().epicsStore()
    
    return ds, epics


def autocorrelate_image(cspad_image):
    """
    Autocorrelate an image to determine the distribution of speckle sizes.
    
    Parameters
    ----------
    cspad_image : np.ndarray
        The image, (32, 185, 388) format

    Returns
    -------
    acf : np.ndarray
        A (`window`, `window`) shape array containing the autocorrelation
        of the image. Not normalized.
        
    Notes
    -----
    Tests indicate current implementation will run at ~2.4s/image = 0.4 Hz.
    """
    
    if not cspad_image.shape[1:] == (185, 388):
        raise ValueError('`cspad_image` incorrect shape. Expected (X, 185, '
                         '388), got %s' % str(cspad_image.shape))
    
    acf = np.zeros((369,775))
    
    for two_by_one in cspad_image:
        x = two_by_one - two_by_one.mean()
        acf += fftconvolve(x, x[::-1,::-1])
        
    acf /= float( cspad_image.shape[0] )
        
    return acf
    
    
def speckle_profile(autocorrelation_image, resolution=1.0):
    """
    Compute the speckle profile from an autocorrelation.
    
    Parameters
    ----------
    autocorrelation_image : np.ndarray
        A (`window`, `window`) shape array containing the autocorrelation
        of the image. Not normalized.
        
    resolution : float
        The resolution of the profile, in pixel units.
        
    Returns
    -------
    profile : np.ndarray
        The radial profile of the speckle autocorrelation.
    """
    
    if not autocorrelation_image.shape == (369, 775):
        raise ValueError()
    
    c = np.array([184, 387])[:,None,None]
    r = np.sqrt( np.sum( np.square(np.mgrid[:369,:775] - c), axis=0 ) )
    
    bins = int(r.max() / resolution)
    
    values, edges = np.histogram(r.flatten(), bins=bins, 
                                 weights=autocorrelation_image.flatten() / np.square(r.flatten() + 1e-100))
    profile = np.zeros(( len(values), 2 ))
    profile[:,0] = edges[:-1] + 0.5 * (edges[1] - edges[0])
    profile[:,1] = values
    
    return profile

    
def ADU_to_photons(cspad_image, cuts=None):
    """
    Given `cuts` that demarkate photon bins in ADU units, transform 
    `cspad_image` from ADU units to photon counts.
    """
    if cuts == None:
         cuts = np.arange(21.0, 1000.0, 37.5)
    return np.digitize(cspad_image.flatten(), cuts).reshape(cspad_image.shape)

    
def fit_negative_binomial(samples, method='ml', limit=1e-16):
    """
    Estimate the parameters of a negative binomial distribution, using either
    a maximum-likelihood fit or an analytic estimate appropriate in the low-
    photon count limit.
    
    Parameters
    ----------
    samples : np.ndarray, int
        An array of the photon counts for each pixel.
    
    method : str, {"ml", "expansion", "lsq"}
        Which method to use to estimate the contrast.
        
    Returns
    -------
    contrast : float
        The contrast of the samples.
        
    sigma_contrast : float
        The first moment of the parameter estimation function.
        
    References
    ----------
    ..[1] PRL 109, 185502 (2012)
    """
    
    k = samples.flatten()
    N = float( len(k) )
    k_bar = np.mean(samples)
    
    if method == 'ml': # use maximium likelihood estimation
        
        def logL_prime(contrast):
            M = 1.0 / contrast
            t1 = -N * (np.log(k_bar/M + 1.0) + digamma(M))
            t2 = np.sum( (k_bar - k)/(k_bar + M) + digamma(k + M) )
            return t1 + t2
       
        try: 
            contrast = optimize.brentq(logL_prime, limit, 1.0)
        except ValueError as e:
            print e
            #raise ValueError('log-likelihood function has no maximum given'
            #                 ' the empirical example provided. Please samp'
            #                 'le additional points and try again.')
            contrast = limit
            sigma_contrast = 1.0
            return contrast, sigma_contrast
        
        def logL_dbl_prime(contrast):
            M = 1.0 / (contrast + 1e-100)
            t1 = np.sum( (np.square(k_bar) - k*M) / (M * np.square(k_bar + M)) )
            t2 = - N * digamma(M)
            t3 = np.sum( digamma(k + M) )
            return t1 + t2 + t3
            
        sigma_contrast = logL_dbl_prime(contrast)
        if sigma_contrast < 0.0:
            raise RuntimeError('Maximum likelihood optimization found a local '
                               'minimum instead of maximum! sigma = %s' % sigma_contrast) 
    
    elif method == 'lsq': # use least-squares fit

        empirical_pmf = np.bincount(samples)
        k_range = np.arange(len(empirical_pmf))

        err = lambda contrast : negative_binomial_pmf(k_range, k_bar, contrast) - empirical_pmf
        
        c0 = 0.5
        contrast, success = optimize.leastsq(err, c0)
        if success:
            contrast = contrast[0]
            sigma_contrast = success
        else:
            raise RuntimeError('least-squares fit did not converge.')
    
    elif method == 'expansion': # use low-order expansion
        # directly from the SI of the paper in the doc string
        p1 = np.sum( k == 1 ) / N
        p2 = np.sum( k == 2 ) / N
        contrast = (2.0 * p2 * (1.0 - p1) / np.square(p1)) - 1.0
        
        # this is not quite what they recommend, but it's close...
        # what they recommend is a bit confusing to me atm --TJL
        sigma_contrast = np.power(2.0 * (1.0 + contrast) / N, 0.5) / k_bar
        
        
    else:
        raise ValueError('`method` must be one of {"ml", "ls", "expansion"}')
    
    return contrast, sigma_contrast
    
    
def negative_binomial_pmf(k_range, k_bar, contrast):
    """
    Evaluate the negative binomial probablility mass function.
    
    Parameters
    ----------
    k_range : ndarray, int
        The integer values in the domain of the PMF to evaluate.
        
    k_bar : float
        The mean count density. Greater than zero.
        
    contrast : float
        The contrast parameter, in [0.0, 1.0).
        
    Returns
    -------
    pmf : np.ndarray, float
        The normalized pmf for each point in `k_range`.
    """

    M = 1.0 / contrast
    norm = np.exp(gammaln(k_range + M) - gammaln(M) - gammaln(k_range+1))
    f1 = np.power(1.0 + M/k_bar, -k_range)
    f2 = np.power(1.0 + k_bar/M, -M)
    
    return norm * f1 * f2



def negative_binomial_ef(k_range, k_bar, contrast):
    """
    Evaluate the negative binomial probablility error function.
    
    Parameters
    ----------
    k_range : ndarray, int
        The integer values in the domain of the PMF to evaluate.
        
    k_bar : float
        The mean count density. Greater than zero.
        
    contrast : float
        The contrast parameter, in [0.0, 1.0).
        
    Returns
    -------
    pmf : np.ndarray, float
        The normalized pmf for each point in `k_range`.
    """
    
    M = 1.0 / contrast
    norm = np.exp(gammaln(k_range + M) - gammaln(M) - gammaln(k_range+1))
    f1 = np.power(1.0 + M/k_bar, -k_range)
    f2 = np.power(1.0 + k_bar/M, -M)
    
    return norm * f1 * f2

    
    
