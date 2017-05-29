

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


data_all= pd.read_csv('aegyptiDENVmodelTempData_2016-03-30.csv')

# Exclude the Focks & Barrera 2006 data because they're from a model
data_all = data_all.query('ref!= "Focks_Barrera_2006_Research&TrainingTropicalDis_Geneva_Paper"')


# PYMC parameters

sample_size= 5000


#############################################
################# GCR #######################
#############################################


#pEA, b, 
data = data_all.query('trait_name== "GCR"')

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

basic_model_GCR = Model()


with basic_model_GCR: 
	#priors for unknown model parameters
	c = Gamma('c',1,10)
	Tm= Uniform('Tm',25,45)
	T0= Uniform('T0',0,24)
	tau= Gamma('tau',0.0001, 0.0001)

	mu_temp= c*T*((T-T0)*(T0<T))*np.sqrt((Tm-T)*(Tm>T))
	mu= 0*(mu_temp<0) + mu_temp*(mu_temp>0)

	Y_obs = Normal('Y_obs',mu=mu, sd=tau, observed= Y)


from pymc3 import Metropolis, sample, find_MAP
from scipy import optimize

with basic_model_GCR:  

    # obtain starting values via MAP
    start = find_MAP(fmin=optimize.fmin_powell)

    # draw 5000 posterior samples

    trace= sample(sample_size, step= Metropolis(), start=start)
  



#thin the samples by selecting every 5 samples
thin_factor=5

#summary(trace)
#traceplot(trace); 


#PLOTTING THE HISTOGRAM

#figure_count= mua.create_2x2_histograms(trace, figure_count)


#Create the Brier Function
Temps= np.arange(0,50,0.1)
a_samps= mua.make_sims_temp_resp("briere",trace, Temps, thin_factor )

samps={}
samps['GCR']=a_samps
samps['T']=Temps

q= mua.temp_sim_quants(a_samps,Temps)


#PLOTTING THE DATA

#figure_count= mua.add_sim_lines(Temps, a_samps, q,Y , T, figure_count)



#CHECK COUNT:*
#############################################
################# EFD #######################
#############################################


#pEA, b, 
data = data_all.query('trait_name== "EFD"')

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

basic_model_EFD = Model()


with basic_model_EFD: 
	#priors for unknown model parameters
	c = Gamma('c',9,1000)
	Tm= Uniform('Tm',25,45)
	T0= Uniform('T0',0,24)
	tau= Gamma('tau',10, 1)

	mu_temp= c*T*((T-T0)*(T0<=T))*np.sqrt((Tm-T)*(Tm>=T))
	mu= 0*(mu_temp<0) + mu_temp*(mu_temp>0)

	Y_obs = Normal('Y_obs',mu=mu, sd=tau, observed= Y)




with basic_model_EFD:  

    # obtain starting values via MAP
    start = find_MAP(fmin=optimize.fmin_powell)

    # draw 5000 posterior samples

    trace= sample(sample_size, step= Metropolis(), start=start)
  



#thin the samples by selecting every 5 samples
thin_factor=5

#summary(trace)
#traceplot(trace); 


#PLOTTING THE HISTOGRAM

#figure_count= mua.create_2x2_histograms(trace, figure_count)


#Create the Brier Function
Temps= np.arange(0,50,0.1)
EFD_samps= mua.make_sims_temp_resp("briere",trace, Temps, thin_factor )
samps['EFD']=EFD_samps


q= mua.temp_sim_quants(EFD_samps,Temps)


#PLOTTING THE DATA

#figure_count= mua.add_sim_lines(Temps, EFD_samps, q,Y , T, figure_count)


q= mua.temp_sim_quants(EFD_samps,Temps)


#PLOTTING THE DATA

#figure_count= mua.add_sim_lines(Temps, EFD_samps, q,Y , T, figure_count)



#CHECK COUNT:*
#############################################
################# b #######################
#############################################

data = data_all.query('trait_name== "b"')
data= data.query('ref!="Lambrects_et_al_2011PNAS"')

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


basic_model_b= Model()

with basic_model_b:
	#priors for unknown model parameters
	c = Gamma('c',7,10000)
	Tm= Uniform('Tm',25,45)
	T0= Uniform('T0',0,24)
	tau= Gamma('tau',10, 1)

	mu_temp= c*T*((T-T0)*(T0<T))*np.sqrt((Tm-T)*(Tm>=T))
	mu= 1*(mu_temp>=1) + mu_temp*(mu_temp>0)*(mu_temp<1) + 0*(mu_temp<=0)

	Y_obs = Normal('Y_obs',mu=mu, sd=tau, observed= Y)


# Specifing the parameters that control the MCMC (these will be used throughout the code). 

with basic_model_b:  
    # obtain starting values via MAP
    start = find_MAP(fmin=optimize.fmin_powell)
    # draw 5000 posterior samples
    trace= sample(sample_size, step= Metropolis(), start=start)

#thin the samples by selecting every 5 samples
thin_factor=5

#PLOTTING THE HISTOGRAM

#figure_count= mua.create_2x2_histograms(trace, figure_count)


#Create the Brier Function
Temps= np.arange(0,50,0.1)
b_samps= mua.make_sims_temp_resp("briere_trunc",trace, Temps, thin_factor )
samps['b']=b_samps

#calculate the quants
q= mua.temp_sim_quants(b_samps,Temps)


#PLOTTING THE DATA

#figure_count= mua.add_sim_lines(Temps, b_samps, q,Y , T, figure_count)



#CHECK COUNT:*
#############################################
################# c #######################
#############################################



data = data_all.query('trait_name== "c"')
data= data.query('ref!="Lambrects_et_al_2011PNAS"')
data= data.query('ref!="Alto&Bettinardi_2013_AJTMH"')

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


basic_model_c= Model()

with basic_model_c:
	#priors for unknown model parameters
	c = Gamma('c',1,10)
	Tm= Uniform('Tm',25,45)
	T0= Uniform('T0',0,24)
	tau= Gamma('tau',0.0001, 0.0001)

	mu_temp= c*T*((T-T0)*(T0<T))*np.sqrt((Tm-T)*(Tm>=T))
	mu= 1*(mu_temp>1) + mu_temp*(mu_temp>0)*(mu_temp<1) + 0*(mu_temp<0)

	Y_obs = Normal('Y_obs',mu=mu, sd=tau, observed= Y)


# Specifing the parameters that control the MCMC (these will be used throughout the code). 

with basic_model_c:  
    # obtain starting values via MAP
    start = find_MAP(fmin=optimize.fmin_powell)
    # draw 5000 posterior samples
    trace= sample(sample_size, step= Metropolis(), start=start)

#thin the samples by selecting every 5 samples
thin_factor=5

#PLOTTING THE HISTOGRAM

#figure_count= mua.create_2x2_histograms(trace, figure_count)


#Create the Brier Function
Temps= np.arange(0,50,0.1)
c_samps= mua.make_sims_temp_resp("briere_trunc",trace, Temps, thin_factor )
samps['c']=c_samps



#calculate the quants
q= mua.temp_sim_quants(c_samps,Temps)


#PLOTTING THE DATA

#figure_count= mua.add_sim_lines(Temps, c_samps, q,Y , T, figure_count)



#CHECK COUNT:*
#############################################
################# MDR #######################
#############################################


#pEA, b, 
data = data_all.query('trait_name== "MDR"')

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

basic_model_MDR = Model()


with basic_model_MDR: 
	#priors for unknown model parameters
	c = Gamma('c',1,10)
	Tm= Uniform('Tm',25,45)
	T0= Uniform('T0',0,24)
	tau= Gamma('tau',0.0001, 0.0001)

	mu_temp= c*T*((T-T0)*(T0<T))*np.sqrt((Tm-T)*(Tm>T))
	mu= 0*(mu_temp<0) + mu_temp*(mu_temp>0)

	Y_obs = Normal('Y_obs',mu=mu, sd=tau, observed= Y)


from pymc3 import Metropolis, sample, find_MAP
from scipy import optimize

with basic_model_MDR:  

    # obtain starting values via MAP
    start = find_MAP(fmin=optimize.fmin_powell)

    # draw 5000 posterior samples

    trace= sample(sample_size, step= Metropolis(), start=start)
  



#thin the samples by selecting every 5 samples
thin_factor=5

#summary(trace)
#traceplot(trace); 


#PLOTTING THE HISTOGRAM

#figure_count= mua.create_2x2_histograms(trace, figure_count)


#Create the Brier Function
Temps= np.arange(0,50,0.1)
MDR_samps= mua.make_sims_temp_resp("briere",trace, Temps, thin_factor )
samps['MDR']=MDR_samps


q= mua.temp_sim_quants(MDR_samps,Temps)


#PLOTTING THE DATA

#figure_count= mua.add_sim_lines(Temps, MDR_samps, q,Y , T, figure_count)

#CHECK COUNT:*
#############################################
################# pEA #######################
#############################################


#pEA, b, 
data = data_all.query('trait_name== "pEA"')

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

basic_model_pEA = Model()


with basic_model_pEA: 
	#priors for unknown model parameters
	c = Gamma('c',1,1)
	Tm= Uniform('Tm',25,45)
	T0= Uniform('T0',0,24)
	tau= Gamma('tau',0.0001, 0.0001)


	mu= -c*((T-T0)*(T0<T))*((T-Tm)*(Tm>T))

	Y_obs = Normal('Y_obs',mu=mu, sd=tau, observed= Y)


from pymc3 import Metropolis, sample, find_MAP
from scipy import optimize

with basic_model_pEA:  

    # obtain starting values via MAP
    start = find_MAP(fmin=optimize.fmin_powell)

    # draw 5000 posterior samples

    trace= sample(sample_size, step= Metropolis(), start=start)
  



#thin the samples by selecting every 5 samples
thin_factor=5

#summary(trace)
#traceplot(trace); 


#PLOTTING THE HISTOGRAM

#figure_count= mua.create_2x2_histograms(trace, figure_count)


#Create the Brier Function
Temps= np.arange(0,50,0.1)
pEA_samps= mua.make_sims_temp_resp("quad_trunc",trace, Temps, thin_factor )

#FIND THE QUANTILES
q= mua.temp_sim_quants(pEA_samps,Temps)
samps['pEA']=pEA_samps


#PLOTTING THE DATA

#figure_count= mua.add_sim_lines(Temps, pEA_samps, q,Y , T, figure_count)



#CHECK COUNT:*
#############################################
################# PDR #######################
#############################################


#pEA, b, 
data1 = data_all.query('trait_name== "PDR"')

#GET THE VALUES
Y1 = data1['trait']
Y1=Y1.as_matrix()
Y1= Y1.flatten()

#Get the temperatures
T1= data1['T']
T1=T1.as_matrix()
T1= T1.flatten()


data2 = data_all.query('trait_name== "EIP"')

#GET THE VALUES
Y2 = data2['trait']
Y2=Y2.as_matrix()
Y2= Y2.flatten()
Y2= np.reciprocal(Y2)

#Get the temperatures
T2= data2['T']
T2=T2.as_matrix()
T2= T2.flatten()

#COMBINE
Y= np.append(Y1,Y2)
T= np.append(T1,T2)



size= len(Y)

#TRANSFORM FROM DATAFRAME TO NP ARRAY FORMAT


# Specifing the parameters that control the MCMC (these will be used throughout the code). 

basic_model_PDR = Model()


with basic_model_PDR: 
	#priors for unknown model parameters
	c = Gamma('c',1,10)
	Tm= Uniform('Tm',25,45)
	T0= Uniform('T0',0,24)
	tau= Gamma('tau',0.0001, 0.0001)

	mu_temp= c*T*((T-T0)*(T0<T))*np.sqrt((Tm-T)*(Tm>T))
	mu= 0*(mu_temp<0) + mu_temp*(mu_temp>0)

	Y_obs = Normal('Y_obs',mu=mu, sd=tau, observed= Y)


from pymc3 import Metropolis, sample, find_MAP
from scipy import optimize

with basic_model_PDR:  

    # obtain starting values via MAP
    start = find_MAP(fmin=optimize.fmin_powell)

    # draw 5000 posterior samples

    trace= sample(sample_size, step= Metropolis(), start=start)
  



#thin the samples by selecting every 5 samples
thin_factor=5

#summary(trace)
#traceplot(trace); 


#PLOTTING THE HISTOGRAM

#figure_count= mua.create_2x2_histograms(trace, figure_count)


#Create the Brier Function
Temps= np.arange(0,50,0.1)
PDR_samps= mua.make_sims_temp_resp("briere",trace, Temps, thin_factor )


#FIND THE QUANTILES
q= mua.temp_sim_quants(PDR_samps,Temps)
samps['PDR']=PDR_samps


#PLOTTING THE DATA

#figure_count= mua.add_sim_lines(Temps, PDR_samps, q,Y , T, figure_count)


import pickle




aeg_samps= samps
with open('Aegypti_IPFC_samps.pickle', 'wb') as handle:
    pickle.dump(aeg_samps, handle, protocol=pickle.HIGHEST_PROTOCOL)
