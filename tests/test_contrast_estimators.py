#!/usr/bin/env python

"""
"""

import numpy as np
from speckle import contrast as sc

from matplotlib import pyplot as plt
from matplotlib import cm

# ----------------- parameters ---------------
sample_size = 100

k_range = 15
k_bar_min = -2

contrasts = np.linspace(0.1, 1.0, 100)
k_bars    = np.logspace(k_bar_min,  1,   100)
# --------------------------------------------

methods = ["ml", "lsq", "expansion"]

shp = (len(contrasts), len(k_bars))


f, axes = plt.subplots(len(methods), 3)
extent  = [k_bar_min, 1, np.max(contrasts), np.min(contrasts)]


for im,method in enumerate(methods):
    print method

    c_hat = np.zeros(shp)
    sigma = np.zeros(shp)
    err   = np.zeros(shp)

    for i,contrast in enumerate(contrasts):
        for j,k_bar in enumerate(k_bars):

            nb_pmf = sc.negative_binomial_pmf(k_range, k_bar, contrast)
            nb_pmf /= np.sum(nb_pmf) + 1e-8
            nb_pmf *= sample_size
            nb_pmf = nb_pmf.astype(np.int)

            print nb_pmf

            c_hat[i,j], sigma[i,j] = sc.fit_negative_binomial_from_hist(nb_pmf,
                                                                        method=method,
                                                                        limit=1e-4)
            err[i,j] = np.abs(c_hat[i,j] - contrast) / contrast

    axes[im,0].imshow(c_hat, interpolation='nearest', extent=extent, aspect='auto', vmin=0, vmax=2)
    axes[im,1].imshow(sigma, interpolation='nearest', extent=extent, aspect='auto')
    axes[im,2].imshow(err,   interpolation='nearest', extent=extent, aspect='auto', 
                      cmap=cm.hot_r, vmin=0, vmax=1)

    for ix in range(3):
        axes[im,ix].set_xlabel(r'$\log_{10} \bar{k}$')
        axes[im,ix].set_ylabel(r'$\beta$')

    print c_hat


axes[0,0].set_title(r'$\hat{\beta}$')
axes[0,1].set_title(r'$\hat{\beta} / \sigma_\beta$')
axes[0,2].set_title(r'Error')



plt.tight_layout()
plt.show()




