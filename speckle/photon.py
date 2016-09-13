
"""
Photon Counting
* fit a histogram to get photon gain
* convert an image into photons given a gain

photon gains will be represented by an array, so that e.g.

>>> gain = [0, 10, 25, 100]

would mean that between 0 and 10 ADU, you get zero photons, 10 to 25 is one,
photon, etc etc. 
"""


def fit_photon_gain(histogram_x, histogram_y, plot=False):
    """
    Fit a photon histogram.

    Parameters
    ----------
    histogram_x : np.ndarray
        The x-axis of a photon histogram, specifically the bin centers, in ADU

    histogram_y : np.ndarray
        The number of photons observed in that bin

    Returns
    -------
    gain : list
        A gain array
    """

    raise NotImplementedError()

    return gain


def bin_image(image, gain):

    return photon_counts

