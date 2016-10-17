
#!/usr/bin/env python

from speckle.phcount import fit_photon_hist
from speckle.droplet import dropletize
from speckle.contrast import fit_negative_binomial

import psana
import numpy as np
from matplotlib import pyplot as plt

SHOTS = 50
ADU_PER_PHOTON = 147.0

ds = psana.DataSource('exp=xcsm9816:run=355:smd')
epix = psana.Detector('epix_1')
ipm6 = psana.Detector('XCS-IPM-gon')
mask = np.load('/reg/d/psdm/xcs/xcsm9816/results/sellberg/masks/mask-droplet_avg_run43_epix1_1000_shots+borders.npy')

adus = []

for ie,evt in enumerate(ds.events()):
    print ie

    if ipm6.channel(evt) is not None:
        if ipm6.channel(evt)[0] < 1.0:
            continue

    img  = epix.calib(evt)
    if img is not None:
        img *= mask
    else:
        continue
    pimg = epix.photons(evt, adu_per_photon=90.0) * mask
    adus.extend( dropletize(img) )

    if ie > SHOTS: break

adus = np.array(adus).flatten()
bins = np.arange(0.5, 10.5) * ADU_PER_PHOTON
photon_counts = np.digitize(adus, bins)

plt.figure()

plt.hist((adus / ADU_PER_PHOTON) + 0.5, 100, histtype='step', lw=2, normed=True)
plt.hist(photon_counts, np.arange(10), histtype='step', lw=2, normed=True)
plt.hist(pimg[pimg>0].flatten(), np.arange(10), histtype='step', lw=2, normed=True)
plt.yscale('log')

plt.legend(['Raw ADUs', 'droplets', 'photonize'])

plt.show()

