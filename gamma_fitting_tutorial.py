import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
import pickle 

with open('Aedes_prior_gamma_fits.pickle', 'rb') as handle:
    hyperparam = pickle.load(handle)



gamma = stats.gamma

h1= hyperparam['EFD']['tau'][0]
h2= hyperparam['EFD']['tau'][1]

a, loc, scale = h1, 0, h2
size = 2000
y = np.random.gamma(a,scale, size=size)
x = np.linspace(0, y.max(), 100)

# plot the histogram
plt.hist(y, normed=True, bins=30)

plt.show()