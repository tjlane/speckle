#!/usr/bin/env python


import psana
import numpy as np
from matplotlib import pyplot as plt

from speckle import acf

#RUN   = 'exp=xcs01116:run=120'
#RUN   = 'exp=xpptut15:run=263'
RUN = 'exp=xcsm9816:run=288:smd'
SHOTS = 20


ds = psana.DataSource(RUN)
epix = psana.Detector('epix_1')

imgs = []
for ie, evt in enumerate(ds.events()):
    print ie
    img = epix.calib(evt)
    if img is not None:
        imgs.append(img)

    if ie > SHOTS: break

a = acf.autocorrelate_image(np.array(img), normalize='xpcs')
p = acf.speckle_profile(a)

plt.figure()
plt.subplot(121)
plt.imshow(np.log(np.abs(a)), interpolation='nearest')
plt.subplot(122)
plt.plot(p[:,0], p[:,1])
plt.show()

