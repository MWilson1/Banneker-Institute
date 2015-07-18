#!/usr/bin/env python

'''
@author: Maurice Wilson

Objective
---------
Given a data set, find the best fitting parameters with Levenberg-Marquardt statistics, 
then find the confidence interval of each parameter with MC and Bootstrapping statistical methods.
Chi^2 calculations will govern the Lev-Mar method.  
#This is applicable because the PDF of the data set abides by a Gaussian distribution (due to random noise generator).
'''


import numpy as np
import matplotlib.pyplot as plt
#import numpy.random as rand
#from pylab import *
#import numpy.random as rand  # this must come after pylab!
#from matplotlib import rc
from scipy.interpolate import interp1d
#from numpy.linalg import inv
#from mpl_toolkits.mplot3d import Axes3D
import math
#from scipy.stats import linregress

plt.ion()  
#plt.clf()
#plt.close('all')


#figs = 1


execfile('Lev_Mar_gauss.py')  # you can now generate data and run any definitions from there


########### good choice of parameters for function

#walk_mcmc(x_vals = np.arange(20), data_sigmas = 2, pars = np.array([10, 500, 2000, 1000],dtype=float), pars_g = np.array([90, 1, 2, 50], dtype = float), par_sigmas=.07, trials = 1000, sim_data_gauss = True, plot=True)
#walk_mcmc(x_vals = np.arange(20), data_sigmas = 2, pars = np.array([10, 500, 2000],dtype=float), pars_g = np.array([90, 200, 2], dtype = float), par_sigmas=np.array([.00009, 0.00009, 0.0009], dtype=float), trials = 2000, sim_data_gauss = False, plot=True)
############
def walk_mcmc(x_vals = np.arange(20), data_sigmas = 3, pars = np.array([10, 500, 2000],dtype=float), pars_g = np.array([90, 1, 2], dtype = float), par_sigmas=np.array([.00009, 0.00009, 0.0009], dtype=float), trials = 1000, sim_data_gauss = False, plot=True):

    if sim_data_gauss == True and (len(pars) != 4 or len(pars_g) !=4):
        print('\nOoops!\nIf you want to simulate data with a Gaussian profile, you MUST use 4 parameters.\n4 numbers for "pars" and 4 numbers for "pars_g"')
        return
    
    figs = 1

    
    #x_vals = np.array([1, 2, 3, 4, 5], dtype=float)
    #x_vals = np.array([-1, -0.5, 0, 0.5, 1], dtype=float)*10.0


    '''
    d_mean = np.array([100])
    pars = np.array([d_mean]) # only 1 param for now
    '''

    #pars = np.array([10, 5, 2, 4],dtype=float)

    #pars_g = np.array([9, 7, 3.5, 2.8], dtype = float)   # [A_g, B_g, C_g, E_g]

    
    pars_trial_counter = np.zeros(len(pars))
    pars_n_accept = np.zeros(len(pars))
    
    pars_tracker = []
    pars_trial_tracker = []
    for parameter_num in xrange(len(pars_g)):
        pars_tracker.append([])
        pars_trial_tracker.append([])
    
    #sig = 2
    #sigmas = np.ones(len(x_vals))*sig
    #sigmas = np.array([0.5, 0.3, 1.1, 0.15, 2.0])
    


    # DEF : Generate phony data
     # 5 points, randn() for the noise
     # for now, only use 1 parameter (D_mean)

    def poly_calc(x_vals, pars):
    
        x_mat = np.resize(x_vals*1.0,(len(pars),len(x_vals)))
        power_mat = np.resize(np.arange(len(pars))*1.0, (len(x_vals), len(pars))).T
        a_mat = np.resize(pars*1.0,(len(x_vals),len(pars))).T
        mats= x_mat**power_mat*a_mat
    
        poly_fit = np.sum(mats, axis=0)
        
        return poly_fit
 

    def poly_gen(figs, x_vals, pars, data_sigmas=1, plot=False):
        
        model_poly = poly_calc(x_vals, pars)

        noise = data_sigmas*np.random.randn()
        y_vals = model_poly + noise   # Data
            
        if plot:
            fig1 = plt.figure(figs)
            figs+=1
            plt.clf()
            plt.scatter(x_vals, y_vals)
            #plt.plot(x_vals, mu_vals)
            plt.plot(x_vals, model_poly)
            plt.title('Simulated Data with and without noise')
            #plt.plot(x_vals, np.poly1d(np.fliplr([a_values])[0])(x_vals))

            #plt.show() # unnecessary for the typical ipython
            
        return y_vals, model_poly, figs




    #DEF : calculate a ln(likelihood)
    def ln_likely(Data, mu, data_sigmas):
        
        constant =  -(1.0/2.0)*np.sum( np.log(data_sigmas*np.sqrt(2*np.pi)) )
        # or constant = -(1.0/2.0)*np.sum( np.log(2*np.pi*(data_sigmas**2)) )
        
        chisq = (1.0/2.0)*np.sum( ((Data - mu)/data_sigmas)**2)
    
        ln_likelihood = constant - chisq

        return ln_likelihood


    #DEF : Prepare histogram values for CDF
    def histog(pars_tracker, param, figs, binning=50, plot=False):
        param_col = pars_tracker[param]
        hist_vals, bin_edges = np.histogram(param_col, bins=binning)
        
        x_hist_vals = np.array([np.average(bin_edges[el: el+2]) for el in xrange(len(bin_edges)-1)]) # list comprehension
        # it's freaking amazing that this "el+2" works!
        '''
        if plot:
        plt.figure(figs)
        plt.clf()
        figs+=1
        plt.hist(pars_tracker[param], bins=binning)
        plt.xlabel(r'a$_'+str(param)+'$ sampling')
        plt.ylabel(r'Amount of trials for a$_'+str(param)+'$')
        '''
          
        # Call the mw_cdf function
        figs = mw_cdf(x_hist_vals, hist_vals, param, figs, plot)

        return figs



    def mw_cdf(x_hist_vals, hist_vals, a_coeff, figs, plot=False):
        max_a = np.sum(hist_vals)
        area_a = np.ones(len(hist_vals))
        for el in xrange(len(hist_vals)):
            area_a[el] = np.sum(hist_vals[0:el+1])


        c_d_f = area_a/max_a
        interp = interp1d(c_d_f, x_hist_vals)
        a_best = interp(0.5)
        a_limits = interp(np.array([0.5 - 0.683/2.0, 0.5 + 0.683/2.0]))

        decim = [math.trunc(np.abs(np.log10(a_best - a_limits[0])))+2, math.trunc(np.abs(np.log10(a_limits[1] - a_best)))+2]

        uncertainties = np.array([round(a_best - a_limits[0], decim[0]), round(a_limits[1] - a_best, decim[1])])

        if plot:
            plt.figure(figs)
            figs += 1
            plt.clf()
            plt.scatter(x_hist_vals, c_d_f, marker='+')
            plt.plot((a_best, a_best), (( c_d_f.max(), 0)), 'g')
            plt.errorbar(a_best, 0.5, xerr=[[uncertainties[0]], [uncertainties[1]]], fmt='^', color='red')
            plt.ylabel('CDF ')
            plt.xlabel(r'a$_'+str(a_coeff)+'$ values')

            plt.title(r'Result: a$_'+str(a_coeff)+' = '+str(round(a_best, np.max(decim)))+'_{-'+str(uncertainties[1])+'}^{+'+str(uncertainties[0])+'}$')

        return figs


    # DEF : mcmc
     # try a D_mean value (the direction)
     # randn()*sigma can determine the direction of the new guess/parameter     

    # CALL likelihood calculator function
 
    # if ln(L_new) > ln(L_old) then take step
     # record new param
     # add to number of acceptances
 
    # elseif random.uniform() > np.exp[ln(L_new) - ln(L_old)] then take step
     # record new param
     # add to number of acceptances

    # else reject step, and record the old parameter

    # N_acceptances should be 40% of the total trials


    #trials = 1000

    #def mw_mcmc(figs, x_vals, trials=1000, par_sigmas=0.5, plot=False):  # this only works for one parameter.  change "mu" to "d_mean" and/or other params to accomodate for multiple pa
    if type(par_sigmas) ==float or type(par_sigmas)==int: par_sigmas = np.ones(len(pars_g))*par_sigmas

    if type(data_sigmas) ==float or type(data_sigmas)==int: data_sigmas = np.ones(len(x_vals))*data_sigmas
             
     # CALL the poly_gen function for Data
    if sim_data_gauss ==False:
         initial_data = poly_gen(figs, x_vals, pars, data_sigmas, plot)
         Data = initial_data[0]
         figs = initial_data[2]
         
       # if Gauss data was desired, call the gauss_gen function
    else:
        Data, model_gauss, figs = gauss_gen(figs, x_vals, pars, data_sigmas, plot)
    
    
     # CALL the LM fitter script to get initial guesses
    LM_info = mw_LMfit(figs, x_vals, Data, pars_g, sim_data_gauss, data_sigmas, plot)

        
    #start_d_mean = 110 #Data[0] #np.average(Data) 

    #mu = start_d_mean

    mu = LM_info[0]

    pars_g_LM = LM_info[1]
    
    figs = LM_info[2]
    
    # LOOP for the random walk
    for trial in xrange(trials):
        lnL_old = ln_likely(Data, mu, data_sigmas)
        
        # calcuate L_new
        
          # choose which parameter to adjust
        pars_g_ind = math.trunc( np.random.uniform()*len(pars_g_LM) )
        
        delta_par = np.random.randn()*par_sigmas[pars_g_ind]

        pars_g_LM[pars_g_ind] = pars_g_LM[pars_g_ind] + delta_par

          # calculate the new mu
        if sim_data_gauss == False:
            mu_stepped = poly_calc(x_vals, pars_g_LM)
            
        elif sim_data_gauss == True:
            mu_stepped = pars_g_LM[0]*np.exp(-(1.0/2.0)*((pars_g_LM[1]-x_vals)/pars_g_LM[2])**2)+pars_g_LM[3]
            
        lnL_new = ln_likely(Data, mu_stepped, data_sigmas)
            
        #print(lnL_new/lnL_old)
        #print(lnL_new, lnL_old)

        pars_trial_counter[pars_g_ind] += 1
        
        
        if lnL_new > lnL_old:
            pars_tracker[pars_g_ind].append(pars_g_LM[pars_g_ind])
            #pars_trial_counter[pars_g_ind] += 1
            mu = mu_stepped
            pars_n_accept[pars_g_ind] += 1
            pars_trial_tracker[pars_g_ind].append(trial)
            
        elif np.random.uniform() < np.exp(lnL_new)/np.exp(lnL_old):  # have the " > " would make more sense if you truly wanted to find the precise value.  But " < " means that you are giving it a better chance at traveling AROUND the best value to get a feel and eventually determine the uncertainties
            pars_tracker[pars_g_ind].append(pars_g_LM[pars_g_ind])
            #pars_trial_counter[pars_g_ind] += 1
            mu = mu_stepped
            pars_n_accept[pars_g_ind] += 1
            pars_trial_tracker[pars_g_ind].append(trial)
            
        else:
            # d_mean/mu stays the same
            pass

    pars_n_accept = np.array(pars_n_accept, dtype=float)
    pars_trial_counter = np.array(pars_trial_counter, dtype=float)
    
    perc_check = pars_n_accept/pars_trial_counter
    if len(perc_check[perc_check  < 0.4]) > 0 or len(perc_check) == 0 or len(np.where(np.isnan(perc_check) == True)[0]) > 0:
        print('\nBad choice of "sig"\n No Plot #'+str(figs)+' for you!')
        plot = False
        
    print('Acceptance rates : '+str(perc_check)+'\n')
    print('Number of acceptances : '+str(pars_n_accept)+'\n')
   
                
    if plot:
        #pars_tracker = np.array(pars_tracker)
        for ind in xrange(len(pars_tracker)):
            '''
            #only plot histogram values AFTER the Burn-in
            check = 0
            stdev_start = np.sqrt(np.sum( (pars_tracker[ind][0] - np.average(pars_tracker[ind][0:]))**2 )/np.float(len(pars_tracker[ind][0:])))
            while check < len(pars_tracker[ind])-1:
                stdev = np.sqrt(np.sum( (pars_tracker[ind][check] - np.average(pars_tracker[ind][check:]))**2 )/np.float(len(pars_tracker[ind][check:])))
                #print(stdev)
                if stdev < 0.3*stdev_start:
                    start_ind = check
                    check = len(pars_tracker[ind])                
                check += 1
            '''
                
            # Plot scatter
            plt.figure(figs)
            figs += 1
            plt.clf()
            plt.scatter(pars_trial_tracker[ind], pars_tracker[ind], marker='+')
            #plt.plot((pars_trial_tracker[ind][start_ind], pars_trial_tracker[ind][start_ind]), (np.max(pars_tracker[ind]),np.min(pars_tracker[ind]) ), 'green')
            plt.ylabel(r'Values of a$_'+str(ind)+'$ used to step')
            plt.xlabel('Trial at which step was taken')     
            #plt.show() # unnecessary for the typical ipython

            
            # Plot histograms        
            plt.figure(figs)
            figs += 1
            plt.clf()            
        
            #binning = math.trunc( len(pars_tracker[ind][start_ind:])*0.8 )
            
            #print(len(pars_tracker[ind]),len(pars_tracker[ind][start_ind:]), binning)

            binning = math.trunc( len(pars_tracker[ind])*.3 )
            
            print(len(pars_tracker[ind]), binning)

            #plt.hist(pars_tracker[ind][start_ind:], bins=binning)
            plt.hist(pars_tracker[ind], bins=binning)
            plt.xlabel(r'Values of a$_'+str(ind)+'$ used to step')
            plt.title(r'Histogram of parameters-values found after LM fitting')
            #plt.show() # unnecessary for the typical ipython
            
            #Find the CDF and uncertainties
            figs = histog(pars_tracker, ind, figs, binning, plot)
    
    return #pars_tracker, pars_trial_counter

'''
# DEF : use likelihood and chi^2 to make matrix
 # use matrix to plot contour


# DEF mcmc

 # for now, start at plot point (0,0)

 # choose direction of step
  # parameters[round( random.uniform()*(len(parameters) - 1) )]

 # decide how large the step is
  # use randn() * sigma   # sigma determines how large the step is

'''





'''


        
def polyfitter(x_vals, a_vals, figs, sigmas, plot=False):
    a_1st_set = np.array([1,2])
    
    x_mat = np.resize(x_vals*1.0,(len(a_1st_set),len(x_vals)))
    power_mat = np.resize(np.arange(len(a_1st_set))*1.0, (len(x_vals), len(a_1st_set))).T
    a_mat = np.resize(a_1st_set*1.0,(len(x_vals),len(a_1st_set))).T
    mats= x_mat**power_mat*a_mat
    
    d_vals = np.sum(mats, axis=0) + sigmas*rand.randn(len(x_vals))
    
    #print(d_vals)
    #d_vals = np.poly1d(np.fliplr([a_1st_set])[0])(x_vals) + sigmas*rand.randn(len(x_vals))

    ln_Likel = linefit(d_vals, sigmas, log_L, row)
    Likel = np.exp(ln_Likel - np.max(ln_Likel))
    
    # Marginalization
    x_Likel = np.sum(Likel, axis=0)
    y_Likel = np.sum(Likel, axis=1)
    
'''
'''
    x_mat = np.resize(x_vals*1.0,(len(a_vals),len(x_vals)))
    power_mat = np.resize(np.arange(len(a_vals))*1.0, (len(x_vals), len(a_vals))).T
    a_mat = np.resize(a_vals*1.0,(len(x_vals),len(a_vals))).T
    mats= x_mat**power_mat*a_mat
    y_vals = np.sum(mats, axis=0) + sigma*rand.randn(len(x_vals))
'''
'''
    #x_a = x_mat.dot(a_mat)  # we are not doing matrix multiplication in this exercise    
    #y_vals = np.poly1d(a_vals)(x_vals) + sigma*rand.randn(len(x_vals))
    #coeffs = np.polyfit(x_vals, y_vals, len(a_vals)-1)
    #y_vals = np.poly1d(coeffs)(x_vals) + sigma*rand.randn(len(x_vals))
    #stats=np.linregress(x_vals, y_vals)
    #uncert = rand.random(len(x_vals))*7
    if plot:
        
        fig1 = plt.figure(figs)
        figs+=1
        plt.clf()
        plt.scatter(x_vals, d_vals)
        #plt.plot(x_vals, mu_vals)
        #plt.plot(x_vals, np.poly1d(a_1st_set)(x_vals))
        plt.xlabel(r'x_vals')
        plt.ylabel(r'd_vals')
        plt.title('Simulated Data')
        
        fig2 = plt.figure(figs)
        figs+=1
        plt.clf()
'''
'''
        #ax1d = fig2.add_subplot(111, projection='3d')
        ax = Axes3D(fig2)
        #ax = plt.axes(projection='3d')

        ax3d = fig2.add_subplot(212, projection = '3d')
        con_plot = ax1d.contour(a_vals[0], a_vals[1], Likel) # contourf() would fill the ellipses with color
        ax1d.clabel(con_plot, inline=True, fontsize=10)
        ax1d.set_xlabel(r'a$_0$ values')
        ax1d.set_ylabel(r'a$_1$ values')
        
        
        ax.contour(a_vals[0], a_vals[1], Likel)
        
        ax3d.set_xlabel(r'a$_0$ values')
        ax3d.set_ylabel(r'a$_1$ values')
        ax1d.set_title('2D Contour of Likelihood Matrix')
'''
'''  
        #x_coords, y_coords = meshgrid(a_vals[0], a_vals[1]) # unnecessary
        con_plot = plt.contour(a_vals[0], a_vals[1], Likel) # contourf() would fill the ellipses with color
        plt.clabel(con_plot, inline=True, fontsize=10)
        plt.xlabel(r'a$_0$ values')
        plt.ylabel(r'a$_1$ values')
        plt.title('2D Contour of Likelihood Matrix')
        
         
        fig3 = plt.figure(figs)
        figs+=1
        plt.clf()
        plt.scatter(a_vals[0], x_Likel)
        plt.xlabel(r'a$_0$ values')
        plt.ylabel(r'Sum of Likelihood-columns')
        plt.title(r'Marginalization of 2D Likelihood Matrix for a$_0$') 

        fig4 = plt.figure(figs)
        figs+=1
        plt.clf()
        plt.scatter(a_vals[1], y_Likel)
        plt.xlabel(r'a$_1$ values')
        plt.ylabel(r'Sum of Likelihood-rows')
        plt.title(r'Marginalization of 2D Likelihood Matrix for a$_1$') 

        
        #!plt.scatter(x_vals, d_vals)
        #plt.plot(x_vals, np.poly1d(np.fliplr([a_vals])[0])(x_vals))
        #!plt.plot(x_vals, mu_vals)
        #plt.errorbar(x_vals, y_vals, yerr = uncert, fmt='s', linewidth= 1.5)
        ##plt.show()
    return Likel, x_Likel, y_Likel, figs




def linefit(d_vals, sigmas, log_L, row): # get the ln(Likelihood)
    for j in a_vals[1]:
        row += 1.0
        col = -1.0
        for i in a_vals[0]:
            col += 1.0
            
            a_set = np.array([i, j])

            #print(str(row)+'  '+str(col))
            mu_vals = a_set[0] + a_set[1]*x_vals
            chi_sq = np.sum(((d_vals - mu_vals)/sigmas)**2)
            log_L[row][col] = -(1.0/2.0)*np.sum(np.log(2*np.pi*sigmas**2)) - (1.0/2.0)*chi_sq
            #log_L[row][col] = -(1.0/2.0)*chi_sq  # this works too.  The constant isn't important
            
    return log_L
'''



'''
def mwpolyfit(x_vals, data, power_mat, sigmas=1):
    if type(sigmas) ==float or type(sigmas)==int: sigmas = np.ones(len(x_vals))*sigmas

    # fill in B array
    # raise x-column to respective power and divide by respective sigma
    B_mat = np.resize(x_vals, (len(power_mat), len(x_vals)))          # the "sigmas" variable is missing!
    B_mat = B_mat**power_mat  # element-wise operations
    B_mat_T = B_mat.T
    A_vec = inv( np.mat(B_mat)*np.mat(B_mat_T) )*np.mat(B_mat)*np.mat(data_vec).T  # A-->nx1 vector
    # matrix dot products
    #A_vec = np.mat(D_vec)*np.mat(B_mat_T) * inv( np.mat(B_mat)*np.mat(B_mat_T) ) # A--> 1xn vector
    #'''
'''
    B_mat = np.resize(x_vals, (len(x_vals), len(power_mat)))
    power_mat=power_mat.T
    B_mat = B_mat**power_mat
    B_mat_T = B_mat.T
    A_vec = inv(np.mat(B_mat_T)*np.mat(B_mat))*np.mat(B_mat_T)*np.mat(D_vec).T  # A-->nx1 vector
'''
'''
    return A_vec

'''
