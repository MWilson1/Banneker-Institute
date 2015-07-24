#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
#import numpy.random as rand
#from pylab import *
import numpy.random as rand  # this must come after pylab!
from matplotlib import rc
from scipy.interpolate import interp1d
from numpy.linalg import inv
#from mpl_toolkits.mplot3d import Axes3D
import math

#from scipy.stats import linregress

plt.ion()
#plt.clf()
#plt.close('all')

#x_vals = np.array([1, 2, 3, 4, 5], dtype=float)
x_vals = np.array([-1, -0.5, 0, 0.5, 1], dtype=float)*10.0

'''
a_0 = np.linspace(0.4,1.6,120)
a_1 = np.linspace(0,3,100)
a_vals = np.array([a_0, a_1])
'''

#a_values = np.array([1,2], dtype=float)
a_values=np.array([1, 0, 2], dtype=float)

#sigmas = np.ones(len(x_vals))*0.5
sigmas = np.array([0.5, 0.3, 1.1, 0.15, 2.0])

#log_L = np.ones((len(a_vals[1]),  len(a_vals[0])))  # a0 is y-intercept corresponding to each column (x-axis in log_L matrix)

#row = -1.0

figs = 1


def mwpoly(x_vals, a_values, figs, sigmas=1, plot=False):
    if type(sigmas) ==float or type(sigmas)==int: sigmas = np.ones(len(x_vals))*sigmas

    x_mat = np.resize(x_vals*1.0,(len(a_values),len(x_vals)))
    power_mat = np.resize(np.arange(len(a_values))*1.0, (len(x_vals), len(a_values))).T
    a_mat = np.resize(a_values*1.0,(len(x_vals),len(a_values))).T
    mats= x_mat**power_mat*a_mat

    noise = sigmas*rand.randn(len(x_vals))
    y_vals = np.sum(mats, axis=0) + noise

    model_y_vals = np.sum(mats, axis=0)

    if plot:
        fig1 = plt.figure(figs)
        figs+=1
        plt.clf()
        plt.scatter(x_vals, y_vals)
        #plt.plot(x_vals, mu_vals)
        plt.plot(x_vals, model_y_vals)
        plt.title('Simulated Data with and without noise')
        #plt.plot(x_vals, np.poly1d(np.fliplr([a_values])[0])(x_vals))
        plt.show()
    return y_vals, power_mat, model_y_vals, figs

####################################################################


# Monte Carlo Bootstrap



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
    
    return A_vec




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
        plt.show() #in most cases unnecessary
        
    return figs



D_P_M_fig = mwpoly(x_vals, a_values, figs, sigmas, plot=False)
init_data = D_P_M_fig[0]
power_mat = D_P_M_fig[1]
figs = D_P_M_fig[3]

data_vec = init_data

trials = 1000

A_trials_mat = np.ones((trials, len(a_values))) #blank mat; ones

for run in xrange(trials):
    
    A_vec = mwpolyfit(x_vals, data_vec, power_mat, sigmas)
    
    if A_vec.shape[1] == 1:
        A_trials_mat[run] = np.array(A_vec.T)[0]
    else:
        A_trials_mat[run] = np.array(A_vec)[0]
        
    noise = sigmas*rand.randn(len(x_vals))
    data_vec = init_data + noise


    
def histog(A_trials_mat, figs, binning=50, plot=False):  # run this after execfile( )
    for a_coeff in xrange(len(a_values)):
        a_col = A_trials_mat.T[a_coeff]
        hist_vals, bin_edges = np.histogram(a_col, bins=binning)
        #a_col_sor = np.array(sorted(a_col))  # bin_edges takes care of this
        
        #x_hist_vals = [np.concatenate((x_hist_vals, np.array([np.average(bin_edges[el: el+2])]))) for el in xrange(len(bin_edges)-1)] # doesn't work
        x_hist_vals = np.array([np.average(bin_edges[el: el+2]) for el in xrange(len(bin_edges)-1)]) # list comprehension
        
        if plot:
          plt.figure(figs)
          plt.clf()
          figs+=1
          plt.hist(A_trials_mat.T[a_coeff], bins=binning)
          plt.xlabel(r'a$_'+str(a_coeff)+'$ sampling')
          plt.ylabel(r'Amount of trials for a$_'+str(a_coeff)+'$')
          plt.show() #in most cases unnecessary
          
        # Call the mw_cdf function
        figs = mw_cdf(x_hist_vals, hist_vals, a_coeff, figs, plot)

    return 


#print(A_trials_mat)
#A_mat[run]= np.array(mwpolyfit(x_vals, a_values, figs, sigmas, plot=False).T
#A_vec = [np.array(mwpolyfit(x_vals, a_values, figs, sigmas, plot=False))[0] for trials in xrange(5)]
