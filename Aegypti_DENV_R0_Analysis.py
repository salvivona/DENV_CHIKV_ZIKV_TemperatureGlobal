

import pandas as pd 
import numpy as np
from pymc3 import Model, Normal, Gamma, summary, Uniform, traceplot, Slice, Metropolis, sample, find_MAP
import math
import scipy as sp
from scipy import optimize

#importing modules containing functions

import temp_functions_all as tfa
import mcmc_utils_all as mua


#importing plotting functions

import matplotlib.pyplot as plt
import pickle




figure_count= 0 

with open('Aegypti_IPFC_samps.pickle', 'rb') as handle:
    samps = pickle.load(handle)

with open('lifespan_samps_AEG_ALB.pickle', 'rb') as handle:
    lfsamps = pickle.load(handle)






a = np.asmatrix(samps['GCR'])
PDR= np.asmatrix(samps['PDR'])
MDR= np.asmatrix(samps['MDR'])
EFD= np.asmatrix(samps['EFD'])
pEA= np.asmatrix(samps['pEA'])
b= np.asmatrix(samps['b'])
c= np.asmatrix(samps['c'])
lf= np.asmatrix(lfsamps['lf_Aeg'])


#assuming all of the variablesa the same dim
sample_len= a.shape

R0= np.zeros((sample_len[0],sample_len[1]))


for i in range(0,sample_len[0]):
	R0[i,:]=mua.myR0(np.squeeze(np.asarray(a[i,:])), np.squeeze(np.asarray(b[i,:])), np.squeeze(np.asarray(c[i,:])), np.squeeze(np.asarray(PDR[i,:])), np.squeeze(np.asarray(MDR[i,:])), np.squeeze(np.asarray(EFD[i,:])), np.squeeze(np.asarray(pEA[i,:])), np.squeeze(np.asarray(lf[i,:])))
	
np.savetxt('R0_Aegypti_DENV.csv', R0, delimiter=',')

'''
Temps= np.arange(0,50,0.1)
R0_mean = np.mean(R0, axis=1)
print(Temps.shape)
print(R0_mean.shape)
R0_mean_temp= np.stack((Temps,R0_mean), axis=-1)
print(R0_mean_temp.shape)
plt.plot(R0_mean_temp[:,0], R0_mean_temp[:,1])
R0_mean_temp.tolist()

plt.show()'''