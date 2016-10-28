
import sys
import numpy as np
import matplotlib.pyplot as plt
from speckle import acf


#test_image = np.random.randn(32, 185, 388)
test_image = np.zeros((32, 185, 388))

ng = 1000
sigma = 150.0

for i in range(ng):
    p = np.random.randint(0, 32)
    x = np.random.randint(0, 388)
    y = np.random.randint(0, 185)
    
    c = np.array([y, x])[:,None,None]
    r2 = np.sum( np.square(np.mgrid[:185,:388] - c), axis=0 )
    gaussian_image = np.exp( - r2 / sigma )
    
    
    test_image[p,:,:] += gaussian_image

cf = acf.autocorrelate_image(test_image, normalize=True)
#cf[184,387] = 0.0


profile = acf.speckle_profile(cf)

plt.figure(figsize=(15,4))
plt.subplot(131)
plt.imshow(test_image[0], interpolation='nearest')
plt.subplot(132)
plt.imshow(cf, interpolation='nearest')
plt.subplot(133)
plt.plot(profile[:,0], profile[:,1])
plt.show()

