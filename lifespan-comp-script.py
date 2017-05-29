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
figure_count= 0 



#############################################
################# lifespan of Aegypti #######################
#############################################

data_all= pd.read_csv('aegyptiDENVmodelTempData_2016-03-30.csv')


sample_size= 5000

data_days = data_all.query('trait_name=="p/days"')
data_1mu = data_all.query('trait_name=="1/mu"')

Y_days=data_days['trait']
Y_days= Y_days.as_matrix()
Y_days= Y_days.flatten()
Y_days= np.reciprocal(Y_days)

T_days= data_days['T']
T_days=T_days.as_matrix()
T_days= T_days.flatten()


Y_1mu= data_1mu['trait']
Y_1mu= Y_1mu.as_matrix()
Y_1mu= Y_1mu.flatten()

T_1mu= data_1mu['T']
T_1mu=T_1mu.as_matrix()
T_1mu= T_1mu.flatten()

Y= np.append(Y_days,Y_1mu)
T= np.append(T_days,T_1mu)


size= len(Y)
#TRANSFORM FROM DATAFRAME TO NP ARRAY FORMAT


# Specifing the parameters that control the MCMC (these will be used throughout the code). 

basic_model_lf_Aeg = Model()


with basic_model_lf_Aeg: 
	#priors for unknown model parameters
	c = Gamma('c',1,1)
	Tm= Uniform('Tm',25,45)
	T0= Uniform('T0',0,24)
	tau= Gamma('tau',0.0001, 0.0001)


	mu= -c*((T-T0)*(T0<T))*((T-Tm)*(Tm>T))

	Y_obs = Normal('Y_obs',mu=mu, sd=tau, observed= Y)


from pymc3 import Metropolis, sample, find_MAP
from scipy import optimize

with basic_model_lf_Aeg:  

    # obtain starting values via MAP
    start = find_MAP(fmin=optimize.fmin_powell)

    # draw 5000 posterior samples

    trace= sample(sample_size, step= Metropolis(), start=start)
  



#thin the samples by selecting every 5 samples
thin_factor=5

#summary(trace)
#traceplot(trace); 


#PLOTTING THE HISTOGRAM

figure_count= mua.create_2x2_histograms(trace, figure_count)


#Create the Brier Function
Temps= np.arange(0,50,0.1)
lf_Aeg_samps= mua.make_sims_temp_resp("quad",trace, Temps, thin_factor )
samps={}
samps['lf_Aeg']=lf_Aeg_samps

#FIND THE QUANTILES
q= mua.temp_sim_quants(lf_Aeg_samps,Temps)


#PLOTTING THE DATA

figure_count= mua.add_sim_lines(Temps, lf_Aeg_samps, q,Y , T, figure_count)


#############################################
################# lifespan of Alb #######################
#############################################


data_all= pd.read_csv('albopictusCHIKVmodelTempData_2016-03-26.csv')
data_all= data_all.query('trait2=="sugar-fed"')

sample_size= 5000

data_prop = data_all.query('trait_name=="prop.dead"')
data_1mu = data_all.query('trait_name=="1/mu"')

Y_prop=data_prop['trait']
Y_prop= Y_prop.as_matrix()
Y_prop= Y_prop.flatten()
Y_prop= -np.log(1-Y_prop)
Y_prop= np.reciprocal(Y_prop)


T_prop= data_prop['T']
T_prop=T_prop.as_matrix()
T_prop= T_prop.flatten()

Y_1mu= data_1mu['trait']
Y_1mu= Y_1mu.as_matrix()
Y_1mu= Y_1mu.flatten()

T_1mu= data_1mu['T']
T_1mu=T_1mu.as_matrix()
T_1mu= T_1mu.flatten()

Y= np.append(Y_prop,Y_1mu)
T= np.append(T_prop,T_1mu)



size= len(Y)
#TRANSFORM FROM DATAFRAME TO NP ARRAY FORMAT


# Specifing the parameters that control the MCMC (these will be used throughout the code). 

basic_model_lf_Alb = Model()


with basic_model_lf_Alb: 
	#priors for unknown model parameters
	c = Gamma('c',1,1)
	Tm= Uniform('Tm',25,45)
	T0= Uniform('T0',0,24)
	tau= Gamma('tau',0.0001, 0.0001)


	mu= -c*((T-T0)*(T0<T))*((T-Tm)*(Tm>T))

	Y_obs = Normal('Y_obs',mu=mu, sd=tau, observed= Y)


from pymc3 import Metropolis, sample, find_MAP
from scipy import optimize

with basic_model_lf_Alb:  

    # obtain starting values via MAP
    start = find_MAP(fmin=optimize.fmin_powell)

    # draw 5000 posterior samples

    trace= sample(sample_size, step= Metropolis(), start=start)
  



#thin the samples by selecting every 5 samples
thin_factor=5

#summary(trace)
#traceplot(trace); 


#PLOTTING THE HISTOGRAM

figure_count= mua.create_2x2_histograms(trace, figure_count)


#Create the Brier Function
Temps= np.arange(0,50,0.1)
lf_Alb_samps= mua.make_sims_temp_resp("quad",trace, Temps, thin_factor )
samps['lf_Alb']=lf_Alb_samps


#FIND THE QUANTILES
q= mua.temp_sim_quants(lf_Alb_samps,Temps)


#PLOTTING THE DATA

figure_count= mua.add_sim_lines(Temps, lf_Alb_samps, q,Y , T, figure_count)


import pickle
with open('lifespan_samps_AEG_ALB.pickle', 'wb') as handle:
    pickle.dump(samps, handle, protocol=pickle.HIGHEST_PROTOCOL)
