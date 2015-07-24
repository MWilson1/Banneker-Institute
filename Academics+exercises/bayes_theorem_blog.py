#!/usr/bin/env python

'''
@author: Maurice Wilson

'''
import matplotlib.pyplot as plt
import numpy as np
plt.ion()

results = np.array([1,1,0,0,0,0,1,1,1,0,0,1,1,1,0,1,1,1,1,1])
#[heads,heads,tails,tails,tails,tails,heads,heads,heads,tails,tails,heads,heads,heads,tails,heads,heads,heads,heads,heads]

#results = np.concatenate((results, np.zeros(100))) # testing
#print(results)


figs = 1
plt.figure(figs)
plt.clf()
#figs+=1
h_vals=np.linspace(0,1,100)

def gauss_pro(centroid=0.25, x_vals=h_vals, width=0.05):
    return np.exp(-(1.0/2.0)*((centroid-x_vals)/width)**2)

def likel(k, n, h_vals):
    return (h_vals**k)*((1-h_vals)**(n-k))

plt.plot(h_vals, gauss_pro(), label='Prior = 0.25')
priors = gauss_pro()
#print(priors)

for d in xrange(len(results)):
    #if d == 4 or d== 9 or d==14 or d==19: #Bad!
    if (d+1)%5 == 0:
        likelies = likel(np.sum(results[0:d]), d+1, h_vals)
        
        priors = priors*likelies
    #if (d+1)%5 == 0:    # correct. but precision
        plt.plot(h_vals, priors/np.max(priors), label='Prior #'+str(d+1)) 
        plt.title('PDF Shifting throughout Iterations')
        plt.xlabel('Fitting Parameter(s)')
        plt.ylabel('Probability')
        plt.show()
        
plt.legend()
