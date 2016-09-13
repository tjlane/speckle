
"""
"""

import os

import numpy as np
from scipy.signal import fftconvolve

try:
    import psana
except ImportError as e:
    print ('Could not import psana, proceeding')


def autocorrelate_image(image):
    """
    Autocorrelate an image to determine the distribution of speckle sizes.
    
    Parameters
    ----------
    image : np.ndarray
        The image or stack of images.

    Returns
    -------
    acf : np.ndarray
        A (`window`, `window`) shape array containing the autocorrelation
        of the image. Not normalized.
    """
    
    if len(image.shape) == 3:
        n_images = image.shape[0]
        img_shp  = np.array(image.shape[1:])
    elif len(image.shape) == 2:
        n_images = 1
        img_shp  = np.array(image.shape)
        image = image.reshape(1, img_shp[0], img_shp[1])
    else:
        raise TypeError('`image` is not a valid shape (must be 2d or 3d)')

    acf = np.zeros(img_shp * 2 - 1)
    
    for i in range(n_images):
        img = image[i]
        x = img - img.mean()
        acf += fftconvolve(x, x[::-1,::-1])
        
    acf /= float( n_images )
        
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

    shp = np.array(autocorrelation_image.shape)
    
    c = np.array(shp - 1)[:,None,None]

    g = np.mgrid[:c[0]*2-1,:c[1]*2-1]
    r = np.sqrt( np.sum( np.square(np.mgrid[:369,:775] - c), axis=0 ) )
    
    bins = int(r.max() / resolution)
    
    values, edges = np.histogram(r.flatten(),
                                 bins=bins,
                                 weights=autocorrelation_image.flatten())
    rvalues, redges = np.histogram(r.flatten(), bins=bins)
    profile = np.zeros(( len(values), 2 ))
    profile[:,0] = edges[:-1] + 0.5 * (edges[1] - edges[0])
    profile[:,1] = values / rvalues
    
    return profile

    
