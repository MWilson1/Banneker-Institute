#!/usr/bin/env python

'''
Author: Maurice Wilson

Concept
---------
Find confidence interval of parameters given using Markov Chain Monte Carlo and Bootstrapping Statistics. 
In theory, this can also find the best fitting parameters themselves if they are not given.  
However, Levenberg-Marquardt Statistics is much more highly recommended for determining the best fitting parameters
(and then use MC Stats to evaluate uncertainties).
The statistical method employed applies the concept of the "Random Walk."
'''


import numpy as np
import math

def MC_walking(model_calc, pars_g, Arguments_model, lnlikely_calc, mu, Arguments_lnlikely, trials, par_sigmas):
    '''Determine best fitting parameters and check Likelihoods around them to find each of their confidence intervals

    Parameters
    ----------
    model_calc : function
                Your model used to fit the data. E.g. model_calc = copy.copy(your_function_name)
    
    pars_g : np.ndarray or list
                Parameters used as priors in the "model_calc" function

    Arguments_model : tuple
                Every argument that the "model_calc" function needs, 
                BESIDES "pars_g" should be here. E.g. Arguments_model = x_values, 
                The comma will yield a 1-element tuple.
    
    lnlikely_calc : function
                Your calculation of the ln(Likelihood) 
                
    mu : any type that is necessary for your "lnlikely_calc" function
        Mu is the result of your model.
    
    Arguments_lnlikely : tuple
                Every argument that the "lnlikely_calc" function needs, 
                BESIDES "mu" should be here
    
    trials : integer or float
                Number of times a ln(likelihood) would be found to evaluate newfound fitting parameters   
        
    par_sigmas : integer, float, list or np.ndarray
                Acts as the most common change/difference from a parameter that will be 
                attempted when the MC walk attempts to step
    
    
    Return
    ------
    pars_tracker : np.ndarray
                 The values of the parameters as they change throughout the trials
    
    pars_accept_counter : np.ndarray
                  Number of times that the parameter actually did accept the step
    
    pars_trial_counter : np.ndarray
                  Number of times that the parameter attempted or proceeded to step

    '''
    # Make a few Error Messages
    

    # Prepare counters and trackers for Random Walk
    pars_trial_counter = np.zeros(len(pars_g))
    pars_accept_counter = np.zeros(len(pars_g))
    
    pars_tracker = []
    pars_trial_tracker = []
    for parameter_num in xrange(len(pars_g)):
        pars_tracker.append([])
        pars_trial_tracker.append([])
        
    # The Random Walk
    for trial in xrange(trials):
        
        # calculate ln(likelihood)
        lnL_old = lnlikely_calc(mu, *Arguments_lnlikely)
        
        # choose which parameter to adjust
        pars_g_ind = math.trunc( np.random.uniform()*len(pars_g) )
        
        # normal distribution around parameter will be chosen
        delta_par = np.random.randn()*par_sigmas[pars_g_ind]

        pars_g[pars_g_ind] = pars_g[pars_g_ind] + delta_par
        
        # calculate new ln(Likelihood) with adjusted parameter
          # calculate the new mu 
        mu_stepped = model_calc(pars_g, *Arguments_model)
        
        lnL_new = lnlikely_calc(mu_stepped, *Arguments_lnlikely)
            
        pars_trial_counter[pars_g_ind] += 1
        
        if lnL_new > lnL_old: # then take the step and use another mu
            pars_tracker[pars_g_ind].append(pars_g[pars_g_ind])
            mu = mu_stepped
            pars_accept_counter[pars_g_ind] += 1
            pars_trial_tracker[pars_g_ind].append(trial)
            
        elif np.random.uniform() < np.exp(lnL_new)/np.exp(lnL_old):  # this must go after the lnL_new > lnL_old condition
            # take the step and use another mu
            pars_tracker[pars_g_ind].append(pars_g[pars_g_ind])
            mu = mu_stepped
            pars_accept_counter[pars_g_ind] += 1
            pars_trial_tracker[pars_g_ind].append(trial)
            
        else: # do not take the step
            # mu stays the same
            pass
        

    pars_accept_counter = np.array(pars_accept_counter, dtype=float)
    pars_trial_counter = np.array(pars_trial_counter, dtype=float)
    
    # check if the percentage of acceptances for each parameter is above 40%
    perc_check = pars_accept_counter/pars_trial_counter

    if len(perc_check[perc_check  < 0.4]) > 0 or len(perc_check) == 0 or len(np.where(np.isnan(perc_check) == True)[0]) > 0:
        print('\nBad choice of parameter-sigmas \n')
        
    print('\nAcceptance rates : '+str(perc_check)+'\n')
    print('\nNumber of accepted steps : '+str(pars_accept_counter)+'\n')
    print('\nNumber of attempted and accepted steps : '+str(pars_trial_counter)+'\n')
    
    return pars_tracker, pars_accept_counter, pars_trial_tracker, pars_trial_counter
