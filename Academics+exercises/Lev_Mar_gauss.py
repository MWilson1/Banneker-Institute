#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
#%matplotlib inline    #       type this if you are in QT Console of ipython
from scipy.optimize import leastsq
#import sys

plt.ion()

'''

     This code should be ran under the "mcmc_walk" script.  But it is powerful enough to be ran under this script alone as long as you call out two of the definitions.
'''

#print(sys.argv)

#figs = 1
'''
Ag = 20
Bg = 7
Cg = 5.
Eg = 9
'''

#eps_data = sigmas


def poly_calc(x_vals, pars):
    
    x_mat = np.resize(x_vals*1.0,(len(pars),len(x_vals)))
    power_mat = np.resize(np.arange(len(pars))*1.0, (len(x_vals), len(pars))).T
    a_mat = np.resize(pars*1.0,(len(x_vals),len(pars))).T
    mats= x_mat**power_mat*a_mat
    
    poly_fit = np.sum(mats, axis=0)
        
    return poly_fit



def gauss_gen(figs, x_vals, pars, data_sigmas=1, plot=False):  # only use this definition if you are running ONLY this script!  Because the model requires precisely 4 parameters
    #x_vals = np.linspace(0,10, 300)
    '''
    A = 10.
    B=5.
    C=2.
    E=4.
    '''
    #sig = 0.25

    #sigmas = np.ones(len(x_vals))*sig

    ##eps_data = sig

    #figs = 1

    A, B, C, E = pars

    noise = data_sigmas*np.random.randn(len(x_vals))  
    Data = A*np.exp(-(1.0/2.0)*((B-x_vals)/C)**2)+E + noise

    model_gauss = A*np.exp(-(1.0/2.0)*((B-x_vals)/C)**2)+E

    #sim_data_gauss = True

    if plot:
        plt.figure(figs)
        figs +=1
        plt.clf()
        plt.scatter(x_vals, Data, marker='o', label='noise added')
        plt.plot(x_vals, model_gauss, label='without noise')
        plt.legend()
        plt.title('Data of Gaussian profile with and without noise')
        plt.xlabel('X values')
        plt.ylabel('Data points')
        plt.show() # in most cases unnecessary
        
    return Data, model_gauss, figs


   
def residual(vars, x_vals, sim_data_gauss, Data, eps_data):
    
    # vars is the pars_g
    
    if sim_data_gauss == True: 
        amp = vars[0]
        centroid = vars[1]
        width = vars[2]
        offset = vars[3]

        new_model = amp*np.exp(-(1.0/2.0)*((centroid-x_vals)/width)**2)+offset # just the 'Data' without the noise
        
    elif sim_data_gauss == False:
        
        pars_g = np.array([])  # this is a redundant part
        for par_g_ind in xrange(len(vars)):
            pars_g = np.concatenate((pars_g, np.array([vars[par_g_ind]])))            
        
        new_model = poly_calc(x_vals, pars_g)    
    
    return (Data-new_model)/eps_data



#def mwgaussfit(x_vals, Data, sigmas, plot=False):
def mw_LMfit(figs, x_vals, Data, pars_g, sim_data_gauss, data_sigmas, plot=False):
    
    #vars = [A_g, B_g, C_g, E_g]
    vars = pars_g
    
    eps_data = data_sigmas
    
    params = leastsq(residual, vars, args=(x_vals, sim_data_gauss, Data, eps_data))[0]
    print('\nLev-Marq best fitting parameters: '+str(params)+'\n')

    if sim_data_gauss == True:
        new_D = params[0]*np.exp(-(1.0/2.0)*((params[1]-x_vals)/params[2])**2)+params[3]
        
    elif sim_data_gauss == False:
        new_D = poly_calc(x_vals, params)

    if plot:
        plt.figure(figs)
        figs += 1
        plt.clf()
        plt.plot(x_vals, new_D, 'r')
        plt.scatter(x_vals, Data)
        plt.xlabel('X values')
        plt.ylabel('Measurement')
        plt.title('Fitting immediately after LM stats')
        plt.show() # unnecessary for the typical ipython
        
    return new_D, params, figs


'''
out = mwgaussfit(x_vals, Data, eps_data)
new_D = out[0][0]*np.exp(-(1.0/2.0)*((out[0][1]-x_vals)/out[0][2])**2)+out[0][3]

plt.plot(x_vals, new_D, 'r')
plt.scatter(x_vals, Data)
'''
