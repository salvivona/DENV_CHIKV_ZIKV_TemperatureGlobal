

import pandas as pd 
import numpy as np
from pymc3 import Model, Normal, Gamma, summary, Uniform, traceplot, Slice
import math
import scipy as sp
np.random.seed(123)


data_all= pd.read_csv('aegyptiDENVmodelTempData_2016-03-30.csv')


# Exclude the Focks & Barrera 2006 data because they're from a model
data = data_all.query('ref!= "Focks_Barrera_2006_Research&TrainingTropicalDis_Geneva_Paper"')

data = data.query('trait_name== "GCR"')

#GET THE VALUES
Y = data['trait']
Y=Y.as_matrix()
Y= Y.flatten()

#Get the temperatures
T= data['T']
T=T.as_matrix()
T= T.flatten()


size= len(Y)
#TRANSFORM FROM DATAFRAME TO NP ARRAY FORMAT


# Specifing the parameters that control the MCMC (these will be used throughout the code). 

basic_model = Model()


with basic_model: 
	#priors for unknown model parameters
	c = Gamma('c',10,1)
	Tm= Uniform('Tm',31,45)
	T0= Uniform('T0',5,24)
	tau= Gamma('tau',0.0001, 0.0001)

	mu_temp= c*T*((T-T0)*(T0<T))*np.sqrt((Tm-T)*(Tm>T))
	mu= 0*(mu_temp<0) + mu_temp*(mu_temp>0)

	Y_obs = Normal('Y_obs',mu=mu, sd=tau, observed= Y)


from pymc3 import Metropolis, sample, find_MAP
from scipy import optimize
trace_copy= {}
with basic_model:  

    # obtain starting values via MAP
    start = find_MAP(fmin=optimize.fmin_powell)

    #instantiate sampler
 

    # draw 5000 posterior samples

    trace= sample(100, step= Metropolis(), start=start)
    trace_copy= trace

thin_factor=2

print(trace['c'][0:9])
trace= trace[:][0::thin_factor]
print(trace['c'][0:9])

#summary(trace)
#traceplot(trace); 













'''

import numpy as np
import matplotlib.pyplot as plt
from pymc3 import Model, Normal, HalfNormal, summary

# Initialize random number generator
np.random.seed(123)

# True parameter values
alpha, sigma = 1, 1
beta = [1, 2.5]

# Size of dataset
size = 10

# Predictor variable
X1 = np.random.randn(size)
X2 = np.random.randn(size) * 0.2

# Simulate outcome variable
Y = alpha + beta[0]*X1 + beta[1]*X2 + np.random.randn(size)*sigma


basic_model = Model()


with basic_model:
    # Priors for unknown model parameters
	alpha = Normal('alpha', mu=0, sd=10)
	beta = Normal('beta', mu=0, sd=10, shape=2)
	sigma = HalfNormal('sigma', sd=1)

	mu = alpha + beta[0]*X1 + beta[1]*X2
	Y_obs = Normal('Y_obs', mu=mu, sd=sigma, observed=Y)

from pymc3 import NUTS, sample
from scipy import optimize

with basic_model:
	trace = sample()

print(trace['alpha'][-100:])

summary(trace)
# Code to add widgets will go here...
'''
