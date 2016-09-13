#!/usr/bin/env python

from speckle.phcount import fit_photon_hist
from speckle.droplet import dropletize
from speckle.contrast import fit_negative_binomial

import psana
import numpy as np
from matplotlib import pyplot as plt

RUN   = 'exp=xcs01116:run=120'
#RUN   = 'exp=xpptut15:run=263'
SHOTS = 100


ds = psana.DataSource(RUN)
epix = psana.Detector('epix100a_ladm')

adus = []

for ie, evt in enumerate(ds.events()):
    print ie
    img = epix.calib(evt)
    if img is not None:
        adus.extend( dropletize(img) )
    if ie > SHOTS: break


bins, diagnostics = fit_photon_hist(adus, plot=True, max_photons=4)

photon_counts = np.digitize(adus, bins)
print 'counts:', np.bincount(photon_counts)
print 'bins:', bins
print diagnostics

contrast, err = fit_negative_binomial(photon_counts, method='ml', limit=1e-4)
print 'contrast estimate: %f +/- %f' % (contrast, err)


