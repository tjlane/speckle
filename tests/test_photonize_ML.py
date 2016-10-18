
#!/usr/bin/env python

from speckle.phcount import fit_photon_hist
from speckle.droplet import dropletize
from speckle.contrast import *

import psana
import numpy as np
from matplotlib import pyplot as plt

SHOTS = 100
ADU_PER_PHOTON = 145.0

ds = psana.DataSource('exp=xcsm9816:run=88:smd')
epix = psana.Detector('epix_1')
ipm6 = psana.Detector('XCS-IPM-gon')
cspad = psana.Detector('cspad2x2_diff')
mask = np.load('/reg/d/psdm/xcs/xcsm9816/results/sellberg/masks/mask-droplet_avg_run43_epix1_1000_shots+borders.npy')
mask = np.logical_not(mask)

summed = 0
photon_hist = np.zeros(6, dtype=np.int)
cspad_intensity = []

for ie,evt in enumerate(ds.events()):

    print ie, summed

    img  = epix.calib(evt)
    if img is not None:
        img *= mask
    else:
        continue
    pimg = epix.photons(evt, adu_per_photon=ADU_PER_PHOTON) * mask
    ph = np.bincount(pimg.flatten())
    photon_hist[:len(ph)] += ph

    summed += 1
    if summed > SHOTS: break


contrast, err = fit_negative_binomial_from_hist(photon_hist, limit=1e-4)
k_range = len(photon_hist)
k_bar = np.sum(np.arange(k_range) * photon_hist) / float(np.sum(photon_hist))
print contrast, err, k_range, k_bar

plt.figure()
plt.bar(np.arange(k_range), photon_hist / float(np.sum(photon_hist)), 0.8, alpha=0.8)
plt.plot(negative_binomial_pmf(k_range, k_bar, contrast), lw=2, color='r')
plt.yscale('log')
plt.xlabel('photons / pixel')
plt.ylabel('Observations')
plt.show()

