#!/usr/bin/env python

'''
Author: Maurice Wilson

Concept
---------
Find confidence interval of parameters given using Markov Chain Monte Carlo and Bootstrapping Statistics. 
In theory, this can also find the best fitting parameters themselves if they are not given.  
However, Levenberg-Marquardt Statistics is much more highly recommended for determining the best fitting parameters
(and then use MC Stats to evaluate uncertainties).

The statistical method employed applies the concept of the "Random Walk" or "Drunk Walk."
'''


import numpy as np
import math

def MC_walking(model_calc, all_priors_index, Arguments_model, lnlikely_calc, mu_index, Arguments_lnlikely, trials, par_sigmas):
    '''Determine best fitting parameters and check Likelihoods around those parameters to find each parameter's confidence interval

    Parameters
    ----------
    model_calc : function
                Your model used to fit the data. 
    
    all_priors_index : int
                This is the index for where the prior parameters are within your "Arguments_model" tuple.  Priors are the parameters chosen preceding the calculation of your model (which occurs in the "model_calc" function).

    Arguments_model : tuple
                Every argument that the "model_calc" function needs including your prior parameters.  E.g. Arguments_model = x_values, 
                The comma will yield a 1-element tuple. 
    
    lnlikely_calc : function
                Your calculation of the ln(Likelihood) 
                
    mu_index : int
                This is the index for where Mu is within your "Arguments_lnlikely." Mu is the result(s) of your model. 
    
    Arguments_lnlikely : tuple
                Every argument that the "lnlikely_calc" function needs including the "Mu"
    
    trials : int or float
                Number of times a ln(likelihood) would be found to evaluate newfound fitting parameters   
        
    par_sigmas : int, float, list or np.ndarray
                Acts as the most common change/difference/displacement from a parameter that will be 
                attempted when the MC walk attempts to step.  
    
    
    Return
    ------
    pars_tracker : list
                 The values of the parameters as they change throughout the trials
    
    pars_accept_counter : np.ndarray
                  Number of times that the parameter actually did accept the step
    
    pars_trial_counter : np.ndarray
                  Number of times that the parameter attempted or proceeded to step

    pars_trial_tracker : list
                  The trial numbers at which the parameter was chosen for a step attempt
    

    Example
    -------
    >>>gauss_calc = lambda ...
    >>>Arguments_model = x_vals, initial_parameters
    >>>all_priors_index = 1
    >>>lnlikely_calc = copy.copy(lnlikely_poisson)
    >>>Arguments_lnlikely = mu, x_vals
    >>>mu_index = 0
    >>>pars_tracker, pars_accept_counter, pars_trial_tracker, pars_trial_counter = MC_walking(gauss_calc, all_priors_index, Arguments_model, lnlikely_calc, mu_index, Arguments_lnlikely, trials, par_sigmas)
    '''
    
    # Make a few Error Messages


    priors = Arguments_model[all_priors_index]
    
    if type(par_sigmas)==int or type(par_sigmas)==float: par_sigmas=par_sigmas*np.ones(len(priors))   

    mu = Arguments_lnlikely[mu_index]
    
    # Prepare counters and trackers for Random Walk
    pars_trial_counter = np.zeros(len(priors))
    pars_accept_counter = np.zeros(len(priors))
    
    pars_tracker = []
    pars_trial_tracker = []
    for parameter_num in xrange(len(priors)):
        pars_tracker.append([])
        pars_trial_tracker.append([])
        
    # The Random Walk
    for trial in xrange(trials):
        
        # calculate ln(likelihood)
        lnL_old = lnlikely_calc(*Arguments_lnlikely)
        
        # choose which parameter to adjust
        priors_ind = math.trunc( np.random.uniform()*len(priors) )
        
        # normal distribution around parameter will be chosen
        delta_par = np.random.randn()*par_sigmas[priors_ind]

        priors[priors_ind] = priors[priors_ind] + delta_par
        
        # calculate new ln(Likelihood) with adjusted parameter
          # calculate the new mu

        Arguments_model = list(Arguments_model)        
        Arguments_model[all_priors_index] = priors

        Arguments_model = tuple(Arguments_model)
        mu_stepped = model_calc(*Arguments_model)

        Arguments_lnlikely = list(Arguments_lnlikely)
        Arguments_lnlikely[mu_index] = mu_stepped

        Arguments_lnlikely = tuple(Arguments_lnlikely)
        lnL_new = lnlikely_calc(*Arguments_lnlikely)
            
        pars_trial_counter[priors_ind] += 1
        
        if lnL_new > lnL_old: # then take the step and use another mu
            pars_tracker[priors_ind].append(priors[priors_ind])
            mu = mu_stepped
            pars_accept_counter[priors_ind] += 1
            pars_trial_tracker[priors_ind].append(trial)
            
        elif np.random.uniform() < np.exp(lnL_new)/np.exp(lnL_old):  # this must go after the lnL_new > lnL_old condition
            # take the step and use another mu
            pars_tracker[priors_ind].append(priors[priors_ind])
            mu = mu_stepped
            pars_accept_counter[priors_ind] += 1
            pars_trial_tracker[priors_ind].append(trial)
            
        else: # do not take the step
            # mu stays the same
            pass

        Arguments_lnlikely = list(Arguments_lnlikely)
        Arguments_lnlikely[mu_index] = mu
        Arguments_lnlikely = tuple(Arguments_lnlikely)

        
    pars_accept_counter = np.array(pars_accept_counter, dtype=float)
    pars_trial_counter = np.array(pars_trial_counter, dtype=float)
    
    # check if the percentage of acceptances for each parameter is above 40%
    perc_check = pars_accept_counter/pars_trial_counter

    if len(perc_check[perc_check  < 0.4]) > 0 or len(perc_check) == 0 or len(np.where(np.isnan(perc_check) == True)[0]) > 0:
        print('\nBad choice of parameter-sigmas, also known as "step size" \n')
        
    print('\nAcceptance rates : '+str(perc_check)+'\n')
    print('\nNumber of accepted steps : '+str(pars_accept_counter)+'\n')
    print('\nNumber of attempted and accepted steps : '+str(pars_trial_counter)+'\n')
    
    return pars_tracker, pars_accept_counter, pars_trial_tracker, pars_trial_counter
