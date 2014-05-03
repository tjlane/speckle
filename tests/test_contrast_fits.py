
import sys
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt

sys.path.append('../speckle')
import core

plt.ion()
fig = plt.figure()

k_bar       = 5.0
contrast    = 0.1
sample_size = 1e4


#class negative_binomial(stats.rv_discrete):
#    def _pmf(self, n, param):
#        return 

def get_randos(k_bar, contrast, samples):
    x = np.arange(100)
    pmf = core.negative_binomial_pmf(x, k_bar, contrast)
    vals = [x, pmf]
    nbinom = stats.rv_discrete(name='nbinom', values=vals)
    return nbinom.rvs(size=samples)


print "k_bar\tk_hat\tbeta\tb_hat\tp\tr\tmu"
print "-----\t"*7

for contrast in np.linspace(0.1, 0.9, 101):

    Bk = k_bar * contrast
    p = Bk / (Bk + 1.0)
    r = 1.0 / contrast
    mu = (p*r) / (1.-p)

    #samples = np.random.negative_binomial(r, p, size=sample_size)
    #samples = stats.nbinom.rvs(r, p, size=sample_size)
    samples = get_randos(k_bar, contrast, sample_size)

    x = np.arange(samples.min(), samples.max())
    beta_hat, err = core.fit_negative_binomial(samples, method='lsq')
    curve = core.negative_binomial_pmf(x, samples.mean(), beta_hat)

    print "%.2f\t"*7 % (k_bar, samples.mean(), contrast, beta_hat, p, r, mu)

    plt.cla()
    plt.hist(samples, histtype='step', normed=True)
    plt.plot(x, curve, color='r')
    fig.canvas.draw()


