
import sys
import argparse
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt

sys.path.append('../speckle')
import core

k_bar       = 1.0
contrast    = 0.5
sample_size = 2*185*388

# ---- parse arguments
parser = argparse.ArgumentParser(description='CONTRAST TESTER')

parser.add_argument('-k', '--k_bar', type=float,
                    help='Average photon count rate [photon/pixel]',
                    default=0.5)
parser.add_argument('-c', '--contrast', type=float,
                    help='Contrast',
                    default=0.5)
parser.add_argument('-n', '--sample_size', type=float,
                    help='Sample size [pixels]',
                    default=1e4)

args = parser.parse_args()

plt.ion()
fig = plt.figure()

#class negative_binomial(stats.rv_discrete):
#    def _pmf(self, n, param):
#        return 

def get_randos(k_bar, contrast, samples):
    x = np.arange(100)
    pmf = core.negative_binomial_pmf(x, k_bar, contrast)
    vals = [x, pmf]
    nbinom = stats.rv_discrete(name='nbinom', values=vals)
    return nbinom.rvs(size=samples)

contrast_error = []
many_photon_events = []
contrasts = np.linspace(0.1, 0.9, 101)
#k_bars = np.power(10, np.linspace(-5, 1, 101))
sample_size = 1e5

#print "k_bar\t\tk_hat\t\tbeta\tb_hat\terr\tp\tr\tmu\tN"
print "k_bar\tk_hat\tbeta\tb_hat\terr\tp\tr\tmu\tN"
#print "-----\t\t"*2 + "-----\t"*7
print "-----\t"*9

for contrast in contrasts:
#for k_bar in k_bars:

    Bk = k_bar * contrast
    p = Bk / (Bk + 1.0)
    r = 1.0 / contrast
    mu = (p*r) / (1.-p)

    #samples = np.random.negative_binomial(r, p, size=sample_size)
    #samples = stats.nbinom.rvs(r, p, size=sample_size)
    samples = get_randos(k_bar, contrast, sample_size)

    x = np.arange(samples.min(), samples.max())
    beta_hat, sigma = core.fit_negative_binomial(samples, method='ml')
    curve = core.negative_binomial_pmf(x, samples.mean(), beta_hat)

    many_photon_events.append(len(samples[samples > 1]))
    err = np.abs(contrast - beta_hat)
    contrast_error.append(err)

    #print "%.2f\t"*8 % (k_bar, samples.mean(), contrast, beta_hat, err, p, r, mu)
    #print ("%f\t"*2 + "%.2f\t"*6 + "%.0f") % (k_bar, samples.mean(), contrast, beta_hat, err, p, r, mu, sample_size)
    print ("%.2f\t"*8 + "%.0f") % (k_bar, samples.mean(), contrast, beta_hat, err, p, r, mu, sample_size)

    plt.cla()
    plt.hist(samples, histtype='step', normed=True)
    plt.plot(x, curve, color='r')
    plt.yscale('symlog')
    fig.canvas.draw()

plt.ioff()
plt.figure()

plt.subplot(121)
#plt.plot(k_bars, contrast_error, 'bo')
plt.plot(sample_sizes, contrast_error, 'bo')
plt.xscale('symlog')
#plt.xlabel("Mean count rate [photons/pixel]")
plt.xlabel("Number of pixels")
plt.ylabel("Error in extracted contrast")
#plt.title("%d randomly drawn samples from a" % (sample_size))
plt.title("negative binomial distribution")

plt.subplot(122)
#plt.plot(k_bars, many_photon_events, 'ro')
plt.plot(sample_sizes, many_photon_events, 'ro')
plt.xscale('symlog')
#plt.xlabel("Mean count rate [photons/pixel]")
plt.xlabel("Number of pixels")
plt.ylabel("Many-photon events (k > 1)")
#plt.title("negative binomial distribution with contrast = %.1f" % (contrast))
plt.title("with contrast = %.1f and count rate = %.2f photons/pixel" % (contrast, k_bar))

plt.show()

#raw_input('Dont close me yet...')
