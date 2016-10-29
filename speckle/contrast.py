
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



def fit_negative_binomial(samples, method='ml', limit=1e-4):
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

    limit : float
        If method == 'ml', then this sets the lower limit for which contrast
        can be evaluated. Note that the solution can be sensitive to this
        value.
        
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
    h = np.bincount(samples)
    return fit_negative_binomial_from_hist(h, method=method, limit=limit)
    
    # copied from 347fe0e --TJL

    # k = samples.flatten()
    # N = float( len(k) )
    # k_bar = np.mean(samples)
    #
    # if method == 'ml': # use maximium likelihood estimation
    #
    #     def logL_prime(contrast):
    #         M = 1.0 / contrast
    #         t1 = -N * (np.log(k_bar/M + 1.0) + digamma(M))
    #         t2 = np.sum( (k_bar - k)/(k_bar + M) + digamma(k + M) )
    #         return t1 + t2
    #
    #     try:
    #         contrast = optimize.brentq(logL_prime, limit, 1.0)
    #     except ValueError as e:
    #         print e
    #         #raise ValueError('log-likelihood function has no maximum given'
    #         #                 ' the empirical example provided. Please samp'
    #         #                 'le additional points and try again.')
    #         contrast = limit
    #         sigma_contrast = 1.0
    #         return contrast, sigma_contrast
    #
    #     def logL_dbl_prime(contrast):
    #         M = 1.0 / (contrast + 1e-100)
    #         t1 = np.sum( (np.square(k_bar) - k*M) / (M * np.square(k_bar + M)) )
    #         t2 = - N * digamma(M)
    #         t3 = np.sum( digamma(k + M) )
    #         return t1 + t2 + t3
    #
    #     sigma_contrast = logL_dbl_prime(contrast)
    #     if sigma_contrast < 0.0:
    #         raise RuntimeError('Maximum likelihood optimization found a local '
    #                            'minimum instead of maximum! sigma = %s' % sigma_contrast)
    #
    #
    # elif method == 'ml-nozero': # use maximium likelihood estimation
    #
    #     def logL_prime(contrast):
    #
    #         M = 1.0 / contrast
    #         p = 1.0 / (k_bar / M + 1.0)
    #
    #         # this is the non-zeros correction to the norm const (d log N / dr)
    #         pr = np.power(1.0 - p, M)
    #         t0 = pr * np.log(1.0 - p) / (pr + 1.0)
    #
    #         t1 = -N * (np.log(k_bar/M + 1.0) + digamma(M))
    #         t2 = np.sum( (k_bar - k)/(k_bar + M) + digamma(k + M) )
    #
    #         return t0 + t1 + t2
    #
    #     try:
    #         contrast = optimize.brentq(logL_prime, limit, 1.0)
    #     except ValueError as e:
    #         print e
    #         #raise ValueError('log-likelihood function has no maximum given'
    #         #                 ' the empirical example provided. Please samp'
    #         #                 'le additional points and try again.')
    #         contrast = limit
    #         sigma_contrast = 1.0
    #         return contrast, sigma_contrast
    #
    #     # def logL_dbl_prime(contrast):
    #     #     M = 1.0 / (contrast + 1e-100)
    #     #     t1 = np.sum( (np.square(k_bar) - k*M) / (M * np.square(k_bar + M)) )
    #     #     t2 = - N * digamma(M)
    #     #     t3 = np.sum( digamma(k + M) )
    #     #     return t1 + t2 + t3
    #
    #     #sigma_contrast = logL_dbl_prime(contrast)
    #     sigma_contrast = 10.0
    #     if sigma_contrast < 0.0:
    #         raise RuntimeError('Maximum likelihood optimization found a local '
    #                            'minimum instead of maximum! sigma = %s' % sigma_contrast)
    #
    #
    # elif method == 'lsq': # use least-squares fit
    #
    #     cut = 0
    #     fit_p = False
    #
    #     empirical_pmf = np.bincount(samples)
    #     k_range = np.arange( len(empirical_pmf) )
    #
    #     def _nb_pmf(M, p):
    #         t1 = np.exp(gammaln(k_range + M) - gammaln(M) - gammaln(k_range + 1))
    #         t2 = np.power(p, k_range)
    #         t3 = np.power(1.0 - p, M)
    #         return t1 * t2 * t3
    #
    #     def err(args):
    #
    #         if fit_p:
    #             M, p, s = args
    #         else:
    #             M, s = args
    #             p = 1.0 / (k_bar/M + 1.0)
    #
    #         if (p < 0.0) or (p > 1.0) or (M < 1.0):
    #             return np.inf
    #
    #         obs = s * empirical_pmf[cut:] / np.sum(empirical_pmf[cut:])
    #         return np.sum(np.square(_nb_pmf(M, p)[cut:] - obs))
    #
    #     if fit_p:
    #         args0 = (2.0, 0.5, 10.0)
    #     else:
    #         args0 = (2.0, 10.0)
    #
    #     res = optimize.minimize(err, args0, method='Powell')
    #
    #     print res.x
    #     contrast = 1.0/res.x[0]
    #     #contrast = 1.0/res.x
    #
    #     sigma_contrast = 1.0
    #
    # elif method == 'expansion': # use low-order expansion
    #     # directly from the SI of the paper in the doc string
    #     p1 = np.sum( k == 1 ) / N
    #     p2 = np.sum( k == 2 ) / N
    #     contrast = (2.0 * p2 * (1.0 - p1) / np.square(p1)) - 1.0
    #
    #     # this is not quite what they recommend, but it's close...
    #     # what they recommend is a bit confusing to me atm --TJL
    #     sigma_contrast = np.power(2.0 * (1.0 + contrast) / N, 0.5) / k_bar
    #
    # elif method == 'felix':
    #     p0 = np.sum( k == 0 ) / N
    #     p1 = np.sum( k == 1 ) / N
    #     contrast = p0 / p1 - 1.0 / k_bar
    #     sigma_contrast = 1.0
    #
    #
    # else:
    #     raise ValueError('`method` must be one of {"ml", "ls", "expansion"}')
    #
    # return contrast, sigma_contrast
    




def fit_negative_binomial_from_hist(empirical_pmf, method='ml', limit=1e-4):
    """
    Just like fit_negative_binomial, but uses the argument `empirical_pmf`,
    simply the empirical distribution of counts. I.e. empirical_pmf[n] should yield
    an integer counting how many samples there were with `n` counts.

    Parameters
    ----------
    samples : np.ndarray, int
        An array of the photon counts for each pixel.
    
    method : str, {"ml", "expansion", "lsq"}
        Which method to use to estimate the contrast.

    limit : float
        If method == 'ml', then this sets the lower limit for which contrast
        can be evaluated. Note that the solution can be sensitive to this
        value.    

    Returns
    -------
    contrast : float
        The contrast of the samples.
        
    sigma_contrast : float
        The first moment of the parameter estimation function

    See Also
    --------
    fit_negative_binomial : function
    """
    
    N = float( empirical_pmf.sum() )
    n = np.arange( len(empirical_pmf) )
    k_bar = float(np.sum(n * empirical_pmf)) / N 
    if k_bar == 0.0:
        return 0.0, 0.0
    pmf = empirical_pmf # just a shorthand
    
    if method == 'ml': # use maximium likelihood estimation

        def logL_prime(M):
            t1 = np.sum([ pmf[n] * digamma(n + M) for n in range(len(pmf)) ])
            t2 = -N * digamma(M)
            t3 = N * np.log(M / (M + k_bar))
            return t1 + t2 + t3
       
        try: 
            M_hat = optimize.brentq(logL_prime, 1.0, 1000.0, maxiter=10000)
            contrast = 1.0 / M_hat
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
            n = np.arange( len(empirical_pmf) )
            t1 = np.sum( pmf * ((np.square(k_bar) - n*M) / (M * np.square(k_bar + M))) )
            t2 = - N * digamma(M)
            t3 = np.sum( pmf * digamma(n + M) )
            return t1 + t2 + t3
            
        sigma_contrast = logL_dbl_prime(contrast)
        if sigma_contrast < 0.0:
            raise RuntimeError('Maximum likelihood optimization found a local '
                               'minimum instead of maximum! sigma = %s' % sigma_contrast) 

    elif method == 'ml_nodiv':

        def neglogL(M):
            r = M / k_bar
            t1 = gammaln(pmf + M) - gammaln(pmf + 1.0) - gammaln(M)
            t2 = - pmf * np.log(1.0 + r)
            t3 = - M *   np.log(1.0 + 1.0/r)
            return -1.0 * np.sum(t1 + t2 + t3)

        M_hat = optimize.brent(neglogL, brack=(1.0,1.0/limit))
        contrast = 1.0 / M_hat
        sigma_contrast = 0.0
            
    
    elif method == 'lsq': # use least-squares fit

        k_range = len(empirical_pmf)
        norm_pmf = empirical_pmf.astype(np.float) / (np.sum(empirical_pmf) + 1.0e-8)

        def lsq_error(args):
            c, scale = args
            return negative_binomial_pmf(k_range, k_bar, c) - scale*empirical_pmf
        
        def brute_error(c):
            a = norm_pmf
            b = negative_binomial_pmf(k_range, k_bar, c)
            return np.sum(np.abs(a - b))

        #contrast = optimize.golden(error, brack=(limit, 0.99999))

        #a0 = (0.5, 1.0/np.sum(empirical_pmf))
        #res, success = optimize.leastsq(error, a0)
        #contrast = res[0]


        contrast = optimize.brute(brute_error, [(limit, 1.0)], Ns=100)

        sigma_contrast = 0.0


    elif method == 'expansion': # use low-order expansion
        # directly from the SI of the paper in the doc string
        p1 = float(empirical_pmf[1]) / N
        p2 = float(empirical_pmf[2]) / N
        contrast = (2.0 * p2 * (1.0 - p1) / np.square(p1)) - 1.0
        
        # this is not quite what they recommend, but it's close...
        # what they recommend is a bit confusing to me atm --TJL
        #sigma_contrast = np.power(2.0 * (1.0 + contrast) / N, 0.5) / k_bar
        sigma_contrast = 0.0
        
        
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

    if type(k_range) == int:
        k_range = np.arange(k_range)
    elif type(k_range) == np.ndarray:
        pass
    else:
        raise ValueError('invalid type for k_range: %s' % type(k_range))

    M = 1.0 / contrast
    norm = np.exp(gammaln(k_range + M) - gammaln(M) - gammaln(k_range+1))
    f1 = np.power(1.0 + M/k_bar, -k_range)
    f2 = np.power(1.0 + k_bar/M, -M)
    
    return norm * f1 * f2


def negative_binomial_samples(k_bar, contrast, size=1):
    M = 1.0 / contrast
    p = 1.0 / (k_bar/M + 1.0)
    samples = np.random.negative_binomial(M, p, size=size)
    return samples


def kbar_plot(samples, filename=None, show=True, k_range=5, k_bars=np.logspace(-3, 0, 100)):
    """
    Generate a "k_bar" plot, which shows the statistics of a series of samples
    along side the expected bounds of a negative binomial distribution.
    
    Parameters
    ----------
    samples : ndarray
        A shape (N, ...) array of N samples. Can have any other number of dimensions,
        typically 1 or 2 (2d image).
        
    filename : str or None
        Pass a string to write the figure to disk.
        
    show : bool
        Show the plot!
    """

    from matplotlib import pyplot as plt
    
    colors = ['black', 'blue', 'red', 'green', 'teal', 'orange']
    n_colors = len(colors)

    min_pmf = np.zeros((len(k_bars), k_range))
    max_pmf = np.zeros((len(k_bars), k_range))

    for i,k_bar in enumerate(k_bars):
        min_pmf[i,:] = negative_binomial_pmf(k_range, k_bar, 1e-9)
        max_pmf[i,:] = negative_binomial_pmf(k_range, k_bar, 1.0)

    plt.figure()

    for sample in samples:
        k_bar = np.mean(sample)
        p = np.bincount(sample.flatten(), minlength=5) / float(np.product(sample.shape))
        for i in range(k_range):
            plt.plot(k_bar, p[i], '.', color=colors[i%n_colors])

    for i in range(k_range):
        plt.plot(k_bars, min_pmf[:,i], color=colors[i%n_colors])
        plt.plot(k_bars, max_pmf[:,i], color=colors[i%n_colors])
        if i > 0:
            plt.text(k_bars[13*i], max_pmf[13*i,i], 'k=%d'%i, color=colors[i%n_colors],
                     bbox=dict(facecolor='white', alpha=0.7, 
                               edgecolor=colors[i%n_colors], boxstyle='round,pad=.3'))


    plt.xscale('log')
    plt.yscale('log')

    plt.xlim([1e-3, 1])
    plt.ylim([1e-6, 1e-1])

    if filename:
        plt.savefig(filename)

    if show:
        plt.show()

    return
