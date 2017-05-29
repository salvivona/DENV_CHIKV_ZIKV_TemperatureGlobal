
import numpy as np
import math

#importing modules containing functions

import temp_functions_all as tfa
import mcmc_utils_all as mua


import matplotlib.pyplot as plt

R0 = np.loadtxt('R0_Aegypti_DENV.csv', delimiter=',')


Temps= np.arange(0,50,0.1)
R0_mean = np.mean(R0, axis=1)
print(Temps.shape)
print(R0_mean.shape)
R0_mean_temp= np.stack((Temps,R0_mean), axis=-1)
print(R0_mean_temp.shape)
plt.plot(R0_mean_temp[:,0], R0_mean_temp[:,1])
R0_mean_temp.tolist()

plt.show()
