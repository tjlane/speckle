
import sys
import numpy as np
import matplotlib.pyplot as plt

sys.path.append('../speckle')
import core

samples = np.random.negative_binomial(10, .9, size=1e5)

x = np.arange(samples.min(), samples.max())
beta_hat, err = core.fit_negative_binomial(samples, method='expansion')
curve = core.negative_binomial_pdf(x, samples.mean(), beta_hat)
print x, beta_hat, curve

plt.figure()
plt.hist(samples, histtype='step', normed=True)
plt.plot(x, curve, color='r')
plt.show()