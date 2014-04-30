
import sys
import numpy as np
import matplotlib.pyplot as plt

sys.path.append('../speckle')
import core

test_image = np.random.randn(32, 185, 388)

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

acf = core.autocorrelate_image(test_image)
acf[184,387] = 0.0


plt.figure()
plt.imshow(acf, interpolation='nearest')
plt.show()


profile = core.speckle_profile(acf)
plt.figure()
plt.plot(profile[:,0], profile[:,1])
plt.show()

