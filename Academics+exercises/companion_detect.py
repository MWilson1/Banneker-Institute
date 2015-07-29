#!/usr/bin/env python

'''
Author: Maurice Wilson

Objective
---------
Answer problem 4 in ps3.pdf
'''

import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import numpy as np
import math

plt.ion()


text_f_info = np.genfromtxt('ps3_q4data.txt', names=True, dtype=None)

metallicities = np.array([text_f_info[ind][1] for ind in xrange(len(text_f_info))])


data = np.array([text_f_info[ind][2] for ind in xrange(len(text_f_info))])

metal_1 = metallicities[np.where(data == 1)]

metal_0 = metallicities[np.where(data == 0)]



#alphas: 0, 0.1

#betas: 1, 3


# call this to start the magic!
def p_posterior(metal_0, metal_1, alphas=np.linspace(0,1,100), betas=np.linspace(0,1,120), plot=True): 
    figs = 1
    Posterior_m = np.ones((len(alphas), len(betas)))
    
    row = -1.0
    for a in alphas:
        row += 1.0
        col = -1.0
        for b in betas:
            col += 1.0            
            
            fracts_0 = model_calc(a, b, metal_0)
            fracts_1 = model_calc(a, b, metal_1)
            
            Posterior_m[row][col] = a*b*likel(fracts_0, fracts_1)

    # Do PDF and then CDF

    if plot:
        # Contour Plot
        plt.figure(figs)
        figs+=1
        plt.clf()
        con_plot = plt.contour(betas, alphas, Posterior_m)
        plt.clabel(con_plot, inline=True, fontsize=10)
        plt.xlabel(r'$\alpha$')
        plt.ylabel(r'$\beta$')
        plt.title('2D Contour of Posterior Matrix')

    figs = pdf_2_cdf(alphas, betas, Posterior_m, figs, plot)
    
    return 



def model_calc(alpha, beta, M):
    fract = alpha*10**(beta*M)
    return fract



def likel(fracts_0, fracts_1):
    return np.prod(fracts_1)*np.prod(1-fracts_0)



def pdf_2_cdf(alphas,betas, PDF_Matrix, figs,plot=False):
    param_conglom = [betas, alphas]
    
    for axi in [0,1]:
        
        param_pdf = np.sum(PDF_Matrix, axis=axi) # axi=0 --> beta
        max_a = np.sum(param_pdf)
        area_a = np.ones(len(param_pdf))
        for el in xrange(len(param_pdf)):
            area_a[el] = np.sum(param_pdf[0:el+1]) # even at last element, this will work

        c_d_f = area_a/max_a

        interp = interp1d(c_d_f, param_conglom[axi])
        param_best = interp(0.5)
        
        param_limits = interp(np.array([0.5 - 0.683/2.0, 0.5 + 0.683/2.0]))

        decim = [math.trunc(np.abs(np.log10(param_best - param_limits[0])))+2, math.trunc(np.abs(np.log10(param_limits[1] - param_best)))+2]

        uncertainties = np.array([round(param_best - param_limits[0], decim[0]), round(param_limits[1] - param_best, decim[1])])

        if plot:
            plt.figure(figs)
            figs += 1
            plt.clf()
            plt.scatter(param_conglom[axi], c_d_f, marker='+')
            plt.plot((param_best, param_best), (( c_d_f.max(), 0)), 'g')
            plt.errorbar(param_best, 0.5, xerr=[[uncertainties[0]], [uncertainties[1]]], fmt='^', color='red')
            plt.ylabel('CDF ')
            if axi==0:
                parameter = 'alpha'
            else:
                parameter = 'beta'
            plt.xlabel(r'$\%s$ values' % (parameter))
            best_and_interval = str(round(param_best, np.max(decim)))
            lower_unc = str(uncertainties[1])
            upper_unc = str(uncertainties[0])
            
            plt.title(r'Result: $\%s = %s_{-%s}^{+%s}$' % (parameter, best_and_interval, lower_unc, upper_unc))

    return figs
