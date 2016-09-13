#!/usr/bin/env python

from speckle.phcount import fit_photon_hist
from speckle.droplet import dropletize

import psana
from matplotlib import pyplot as plt

RUN   = 'exp=xcs01116:run=120'
#RUN   = 'exp=xpptut15:run=263'
SHOTS = 20


ds = psana.DataSource(RUN)
epix = psana.Detector('epix100a_ladm')

adus = []

for ie, evt in enumerate(ds.events()):
    print ie
    img = epix.calib(evt)
    if img is not None:
        adus.extend( dropletize(img) )
    if ie > SHOTS: break


fit_photon_hist(adus, plot=True, max_photons=4)

