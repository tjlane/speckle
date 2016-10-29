
"""
Droplet Algorithms
"""

import numpy as np
import scipy.ndimage.measurements as smt 
import scipy.ndimage.morphology as smf
from scipy.ndimage import filters

try:
    from skbeam.core.accumulators.droplet import dropletfind, dropletanal
    _SKBEAM = True
except ImportError as e:
    _SKBEAM = False
    
    
    
def photonize(img, adu_per_photon, remainder_required=0.9, remainder_min=0.5):
    """
    
    """
    
    img[img < 0.0] = 0.0
    full_photons = (img // adu_per_photon).astype(np.int)
    remainder    = img % adu_per_photon
    
    neighbour_max = filters.maximum_filter(remainder, 
                                           footprint=np.array([[0,1,0],
                                                               [1,0,1],
                                                               [0,1,0]]))
    neighbour_max[neighbour_max < remainder_min] = 0.0 # filter small values
    
    neighbour_max_wc = filters.maximum_filter(remainder, 
                                              footprint=np.array([[0,1,0],
                                                                  [1,1,1],
                                                                  [0,1,0]]))
    local_maxima = (remainder == neighbour_max_wc)
    split_photons = ((remainder + neighbour_max) > remainder_required) *\
                    local_maxima
    
    
    photon_img = full_photons + split_photons
    
    return photon_img


def dropletize(img, threshold=10.0, dilate=1, return_all=False):
    """
    A simple and very effective threshold-based droplet algorithm.

    Works only in the case where droplets form disjoint regions (sparse case).
    The algorithm thresholds the image and looks for connected regions above
    the threshold. Each connected region becomes a droplet.

    Parameters
    ----------
    img : np.ndarray
        The two-D image to search for droplets in.

    threhsold : float
        The threshold for the image. Should be between 0 and 1/2 a photon
        in detector gain units.

    dilate : int
        Optionally extend the droplet regions this amount around the thresholded
        region. This lets you be sure you capture all intesity if using
        an aggresive threshold.

    return_coms : bool
        Whether or not to return the droplet positions (centers of mass).

    Returns
    -------
    adus : list of floats
        The summed droplet intensities in detector gain units (ADUs).

    coms : np.ndarray
        The x,y positions of each droplet found.
    """

    bimg = (img > threshold)
    if dilate > 0:
        bimg = smf.binary_dilation(bimg, iterations=dilate)
    limg, numlabels = smt.label(bimg)

    adus = smt.sum(img, labels=limg, index=np.arange(2,numlabels))

    if return_all:
        coms = np.array(smt.center_of_mass(img, labels=limg,
                                        index=np.arange(2,numlabels))) 
        size = smt.sum(np.ones_like(img), labels=limg, index=np.arange(2,numlabels))
        return adus, coms, size

    else:
        return adus


def skbeam_dropletize(img, threshold=10.0):
    """
    Scikit-beam's droplet algorithm. This is originally from Mark Sutton.
    """
    if _SKBEAM is False:
        raise ImportError('You need the droplet branch of skbeam installed')
    bimg = (img > threshold).astype(np.int)
    npeaks, limg = dropletfind(bimg) 
    npix, xcen, ycen, adus, idlist = dropletanal(img.astype(np.int), 
                                                 limg, npeaks) 
    return adus, (xcen, ycen)

