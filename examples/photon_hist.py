#!/usr/bin/env python

from speckle.phcount import fit_photon_hist
from speckle.droplet import dropletize
from speckle.contrast import fit_negative_binomial

import psana
import numpy as np
from matplotlib import pyplot as plt

RUN    = 'exp=xcsm9816:run=400:smd'
SHOTS  = 1000
DILATE = 1
max_photons = 4

ds = psana.DataSource(RUN)
epix = psana.Detector('epix_1')

adus  = []
coms  = []
sizes = []

for ie, evt in enumerate(ds.events()):
    print ie
    img = epix.calib(evt)
    if img is not None:
        ret = dropletize(img, return_all=True, dilate=DILATE)
        adus.extend(ret[0])
        coms.extend(ret[1])
        sizes.extend(ret[2])

    if ie > SHOTS: break

try:
    bins, diagnostics = fit_photon_hist(adus, plot=True, max_photons=max_photons)
except:
    bins = np.arange(-0.5, max_photons+0.5) * 147.0

photon_counts = np.digitize(adus, bins)
print 'counts:', np.bincount(photon_counts)
print 'bins:', bins

plt.figure(figsize=(9,5), facecolor='white')
plt.suptitle(RUN + ' | %d shots | dilation=%d' % (SHOTS, DILATE)) 

plt.subplot(121)
plt.plot(adus, sizes, '.')
plt.xlabel('droplet intensity (ADU)')
plt.ylabel('droplet size (px)')

plt.subplot(122)
plt.plot(photon_counts, sizes, '.')
plt.xlabel('droplet intensity (photons)')
plt.ylabel('droplet size (px)')


plt.show()


contrast, err = fit_negative_binomial(photon_counts, method='ml', limit=1e-4)
print 'contrast estimate: %f +/- %f' % (contrast, err)


