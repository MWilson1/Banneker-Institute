#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from timeit import default_timer

# circle radius = pi*rad^2

rad = 1  # unit circle

# square length = 2*rad

square_l = 2.0*rad

circle_center = np.array([1.0,1.0])

#generated_darts = [square_l * np.random.random_sample((N_pts,2)) for N_pts in [100, 1000, 10000, 1e6]]


'''
generated_darts = Table(names=('100', '1000', '10000', '1e6'))
for N_pts in [100, 1000, 10000, 1e6]:
    generated_darts[str(N_pts)] = square_l * np.random.random_sample((N_pts,2))
    timer = time.time()
'''

def generate_darts(N_pts):
    arr = square_l * np.random.random_sample((2,N_pts)) - 1.0   
    return arr


def in_the_circle(arr):
    N_c_boolean = np.sqrt(arr[0]**2 + arr[1]**2) < 1.0
    N_c_pos = np.array([arr[0][N_c_boolean], arr[1][N_c_boolean]])
    N_inside = len(N_c_pos[0])
    Pi_val = 4*N_inside/N_pts
    return N_c_pos, N_inside, Pi_val
    #return N_inside




times = []
plt.ion()
fig = plt.figure()
x_pt = 0

trials = ['100', '1000', '10000', '1e6']

for N_pts in [100, 1000, 10000, 1e6]: 
    start = default_timer()    
    generated_darts = generate_darts(N_pts)
    N_c_pos_N_in_Pi = in_the_circle(generated_darts)
    ending = default_timer() - start
    ax1 = fig.add_subplot(111)
    ax1.scatter([x_pt], [ending])
    x_pt += 1

plt.xticks(xrange(x_pt+1), trials)
plt.show()


plt.figure()
fig1 = plt.gcf()
generated_darts_100 = generate_darts(100)
circ = plt.Circle((0,0), 1, color='r', fill=False)
plt.scatter(generated_darts_100[0], generated_darts_100[1])
fig1.gca().add_artist(circ)

plt.show()


for itera in xrange(50):
    generated_N = generate_darts(100)
    his_N_c_pos_N_in_Pi = in_the_circle(generated_darts)
    his_N_c_pos_N_in_Pi[1]
    


