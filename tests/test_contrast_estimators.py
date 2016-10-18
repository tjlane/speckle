#!/usr/bin/env python

"""
"""

import numpy as np
from speckle import contrast as sc

from matplotlib import pyplot as plt
from matplotlib import cm
from matplotlib.ticker import FormatStrFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable

# ----------------- parameters ---------------
sample_size = 1000
no_zero_samples = False

k_range = 15
k_bar_min = -2

contrasts = np.linspace(0.01, 0.99, 99) # 99
k_bars    = np.logspace(k_bar_min,  1, 100) # 100
# --------------------------------------------

methods = ["ml", "lsq", "expansion"]

shp = (len(contrasts), len(k_bars))
extent  = [k_bar_min, 1, np.max(contrasts), np.min(contrasts)]

fig, big_axes = plt.subplots( figsize=(10.0, 10.0) , nrows=len(methods), ncols=1, sharey=True)
fig.suptitle('Contrast Estimator Performance (%d samples)' % sample_size)

for im,method in enumerate(methods):

    print method
    big_ax = big_axes[im]
    big_ax.set_title(method + '\n', fontsize=16)
    big_ax.tick_params(labelcolor=(1.,1.,1., 0.0), top='off', bottom='off', left='off', right='off')
    big_ax._frameon = False

    c_hat = np.zeros(shp)
    sigma = np.zeros(shp)
    err   = np.zeros(shp)

    for i,contrast in enumerate(contrasts):
        for j,k_bar in enumerate(k_bars):

            if no_zero_samples:
                ss = sample_size * 100
                samples = []
                while len(samples) < sample_size:
                    ss *= 2
                    samples = sc.negative_binomial_samples(k_bar, contrast, size=ss)
                    samples = samples[ samples > 0 ]
                    samples[:sample_size]
            else:
                samples = sc.negative_binomial_samples(k_bar, contrast, size=sample_size)

            c_hat[i,j], sigma[i,j] = sc.fit_negative_binomial(samples, method=method)
            err[i,j] = np.abs(c_hat[i,j] - contrast) / contrast

    axes = [fig.add_subplot(len(methods), 3, i + 3*im + 1) for i in range(3)]
    im1 = axes[0].imshow(c_hat, interpolation='nearest', extent=extent, aspect='auto', vmin=0, vmax=1)
    axes[0].set_title(r'$\hat{\beta}$')

    im2 = axes[1].imshow(sigma, interpolation='nearest', extent=extent, aspect='auto', vmin=0, vmax=5)
    axes[1].set_title(r'$\hat{\beta} / \sigma_\beta$')

    im3 = axes[2].imshow(err,   interpolation='nearest', extent=extent, aspect='auto', 
                        cmap=cm.hot_r, vmin=0, vmax=1)
    axes[2].set_title(r'Error (%)')

    ims = [im1, im2, im3]
    for iax,ax in enumerate(axes):
        ax.set_xlabel(r'$\log_{10} \bar{k}$')
        ax.set_ylabel(r'$\beta$')
        ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))

        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="10%", pad=0.05)
        plt.colorbar(ims[iax], cax=cax, format="%.1f")

    print c_hat

fig.set_facecolor('w')
plt.tight_layout()
plt.show()




