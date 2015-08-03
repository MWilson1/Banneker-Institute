#!/usr/bin/env python

import numpy as np
import math
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from numpy.linalg import inv
#from matplotlib import rc
import copy


'''
Eventual Goals: 
Create generalized script for performing the matrix multiplications to find A_vec
Create generalized script for performing Lev-Mar Statistics
Create generalized script for performing Monte Carlo and Bootstrapping Statistics

Give user the ability to choose which ones to use 

'''    


def poisson_fit(init_pars_g=[1,1,1,1], par_sigmas=1, trials = 1000, x_vals=np.array([1,3,4,5,6,7,8,10,11,12], dtype=float), Data=np.array([1,2,4,6,3,3,1,1,2,1], dtype=float), plot=False):

    # x_vals : pixels
    # Data   : photons
    
    figs = 1

    if type(par_sigmas) ==float or type(par_sigmas)==int: par_sigmas = np.ones(len(init_pars_g))*par_sigmas

    if type(init_pars_g) == list: init_pars_g = np.array(init_pars_g)
    
    mu = gauss_calc(x_vals, init_pars_g)

    '''
    from Wilson_Lev_Mar_Stat import LM_fitter
    
    # Because of poisson, I won't use LM fitting.

    # Levenberg-Marquardt to find the precise parameters, and thus the mu
    
    model_calc = copy.copy(gauss_calc)
    Arguments_model = x_vals,
    lnlikely_calc = copy.copy(lnlikely_poisson)
    Arguments_lnlikely = x_vals,
    
    #execfile('Lev-Mar_Stat.py')    
      # from Lev-Mar_Stat, run the DEF LM_fitter
      #eps_data = 1

      
    #LM_pars, new_mu = LM_fitter(init_pars_g, Arguments_model, model_calc, Data, 1, mu, Arguments_lnlikely, lnlikely_calc, resid=False)
    LM_pars, new_mu = LM_fitter(init_pars_g, Arguments_model, model_calc, Data, 1, mu, Arguments_lnlikely, lnlikely_calc)

    print(LM_pars)
     
      # plot the new_mu with the original data
    '''
    
    '''
    # Monte Carlo and Bootstrapping  to find the uncertainties

    execfile('MC_Bootstrapping.py')
    
    pars_tracker, mc_mu = MC_walking(model_calc, *Arguments_model, likely_calc, *Arguments_likely, trials, pars_g, pars_g)        
        
     # create histogram:  Amount of trials for parameter vs. parameter value

     # with histogram find CDF and best fitting parameter and confidence interval        
    '''
    
    '''
    try:
        del(MC_walking)
    except:
        pass
    '''
    
    # execute MC and Bootstrapping script
    from Wilson_MCMC_Boot import MC_walking
    
    #model_calc = copy.copy(gauss_calc)
    Arguments_model = x_vals, init_pars_g    # the comma is necessary to keep everything inside of a tuple. e.g. args = x_vals,
    all_priors_index = 1
    lnlikely_calc = copy.copy(lnlikely_poisson)
    Arguments_lnlikely = mu, x_vals
    mu_index = 0
    
    #pars_tracker, pars_accept_counter, pars_trial_tracker, pars_trial_counter = MC_walking(model_calc, LM_pars, Arguments_model, lnlikely_calc, new_mu, Arguments_lnlikely, trials, par_sigmas)
    pars_tracker, pars_accept_counter, pars_trial_tracker, pars_trial_counter = MC_walking(gauss_calc, all_priors_index, Arguments_model, lnlikely_calc, mu_index, Arguments_lnlikely, trials, par_sigmas)
    

    if plot:
        # plot initial data
        plt.figure(figs)
        figs += 1
        plt.clf()
        plt.scatter(x_vals, Data, marker='+')
        plt.ylabel('# of Photons')
        plt.xlabel('Pixel position')
        plt.show() # in most cases unnecessary
        
        # plot MC parameter changes
        for ind in xrange(len(pars_tracker)):
            plt.figure(figs)
            figs+= 1
            plt.clf()
            plt.scatter(pars_trial_tracker[ind], pars_tracker[ind], marker='+')
            plt.ylabel(r'Values of a$_'+str(ind)+'$ used to step')
            plt.xlabel('Trial at which step was taken')
            plt.show() # in most cases unnecessary
            
    return




def lnlikely_poisson(mu, x_vals):
    '''
    probs = []
    [probs.append([]) for i in xrange(len(x_vals))]
    probs = [(mu**x_vals_i*np.exp(-mu))/math.factorial(x_vals_i) for x_vals_i in x_vals]
    
    lnLikelihood = np.log(np.prod(probs))
    #print(np.prod(probs))  # probs are zeros
    #print(lnLikelihood)
    '''
    lnLikelihoods = np.array([])
    for i in xrange(len(x_vals)):
        #print(mu[i])
        try:
            evil_fact_term = np.log(math.factorial(x_vals[i]))  
        except AttributeError:
            evil_fact_term = x_vals[i]*np.log(x_vals[i]) - x_vals[i] # Stirling's approximation
        
        lnLikelihoods =  np.concatenate((lnLikelihoods,  [ x_vals[i]*np.log(mu[i]) - mu[i] - evil_fact_term ]  ))
    return np.sum(lnLikelihoods)



def gauss_calc(x_vals, pars):

    A = pars[0]
    B = pars[1]
    C = pars[2]
    E = pars[3]
    
    return A*np.exp(-(1.0/2.0)*((B-x_vals)/C)**2)+E


    
'''
def mw_cdf(N, probab, a_coeff, figs, plot=False):
    max_a = np.sum(probab)
    area_a = np.ones(len(probab))
    for el in xrange(len(probab)):
        area_a[el] = np.sum(probab[0:el+1])


    c_d_f = area_a/max_a
    interp = interp1d(c_d_f, N)
    a_best = interp(0.5)
    a_limits = interp(np.array([0.5 - 0.683/2.0, 0.5 + 0.683/2.0]))

    decim = [math.trunc(np.abs(np.log10(a_best - a_limits[0])))+2, math.trunc(np.abs(np.log10(a_limits[1] - a_best)))+2]

    uncertainties = np.array([round(a_best - a_limits[0], decim[0]), round(a_limits[1] - a_best, decim[1])])

    if plot:
        plt.figure(figs)
        figs += 1
        plt.clf()
        plt.scatter(N, c_d_f, marker='+')
        plt.plot((a_best, a_best), (( c_d_f.max(), 0)), 'g')
        plt.errorbar(a_best, 0.5, xerr=[[uncertainties[0]], [uncertainties[1]]], fmt='^', color='red')
        plt.ylabel('CDF ')
        plt.xlabel(r'a$_'+str(a_coeff)+'$ values')

        plt.title(r'Result: a$_'+str(a_coeff)+' = '+str(round(a_best, np.max(decim)))+'_{-'+str(uncertainties[1])+'}^{+'+str(uncertainties[0])+'}$')

    return figs
'''


'''
def poisson(Data=np.array([74, 15, 11, 4, 3, 1], dtype=float), N=np.arange(6, dtype=float), coeffs=3, plot=False):
    figs = 1

    power_mat = np.resize(np.arange(coeffs)*1.0, (len(N), coeffs)).T
    
    A_vec = mwpolyfit(N, Data, power_mat)

    if A_vec.shape[1] == 1:
        pars_g = np.array(A_vec.T)[0]
    else:
        pars_g = np.array(A_vec)[0]
        
    mu = poly_calc(N, pars_g)
    
    Likelihood, probs = likely_poisson(mu, N)

    # Monte Carlo and Bootstrapping
    
     # Use likelihood to determine of parameter should be changed

    
    if plot:
        plt.figure(figs)
        figs += 1
        plt.clf()
        plt.scatter(N, Data)
        plt.plot(N, mu)
        plt.title(r'Data and $\mu$')
    
        
    # create histogram:  Amount of trials for parameter vs. parameter value

    # with histogram find CDF and best fitting parameter and confidence interval

    #figs = mw_cdf(N, probs[ind], ind, figs, plot)
        

    return

'''
