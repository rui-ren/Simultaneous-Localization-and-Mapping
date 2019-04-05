# -*- coding: utf-8 -*-
"""
Created on Sat Mar 30 14:36:12 2019

@author: rur4893
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Mar 22 16:54:55 2019

@author: rur4893
"""

import matplotlib.pyplot as plt
from scipy.integrate import quad
import numpy as np
import pandas as pd
import math

# input the basic drilling parameters
m = 0.65                 # flow behavior index
K = 1.85                 # fluid consistency
tau_y = 5.62             # yield stress pa
r_w = 311.2 /2 /1000     # wellbore radius

# define the dimensionless function for this calculation
a = lambda w: (2*m+1)/(m+1)*(2*r_w)/w*tau_y/delta_p
f = lambda x: 2.**((m+1)/m)*x*(((x**(1-m))-1)/(1-m))**(1/m)/(1-a*(x-1))**(1/m)

delta_ps = [4.2*10**6, 3.2*10**6, 2.2*10**6]
# calculation time step
# dynamic list append

Rad_ds = {}
Time_ds = {}

for delta_p in delta_ps:
    Rad_ds['%i' % delta_p] = []
    Time_ds['%i' % delta_p] = []
    
# iterate the pressure
for j, delta_p in enumerate(delta_ps):
    # dimensionless lost circulation parameters
    # Dimensionless_a = [0.0001,0.008,0.001]  The original dimensionless group
    # Dimensionless_a = [0.0001,0.001,0.008]

    #width = [0.0055, 0.00055, 0.000055]
    #width = [0.0008]
    
#    w = lambda a: (2*m+1)/(m+1)*(2*r_w)/a*tau_y/delta_p
#    width = [w(a) for a in Dimensionless_a]
#    a = lambda w: (2*m+1)/(m+1)*(2.*r_w)/w*tau_y/delta_p
 #   Dimensionless_a = [a(w) for w in width ]
    Dimensionless_a = [0.0001,0.0002,0.0003,0.0004,0.0005,0.0006,0.0007,0.0008,0.0009,
                   0.001,0.002,0.003,0.004,0.005,0.005,0.006,0.007,0.008,0.009,0.01,0.02,0.03,0.04]
    
    for k,v in enumerate(Dimensionless_a):

        Radius_d_init = 1.1                            # the initial R_d
        Radius_d_alti = 1 + 1./v                       # the altimate R_d

        Radius_d = np.arange(Radius_d_init, Radius_d_alti)
        total_time_step = len(Radius_d)
        # calculation the time step matrices
        Rad_ds['%i' % delta_p].append(Radius_d);

    # seven matrices
    # use quad integration

    for k,a in enumerate(Dimensionless_a):
        Y = [quad(f,1.01,int_val)[0] for int_val in Rad_ds['%i' % delta_p][k]]
        Time_ds['%i' % delta_p].append(Y)

R_d = pd.DataFrame(Rad_ds)
T_d = pd.DataFrame(Time_ds)

fig, ax = plt.subplots()
ax.xaxis.grid(True, which = 'Major', linestyle='dotted')
ax.yaxis.grid(True, which = 'Major', linestyle='dotted') 

color = ['r', 'b', 'g']
col = 0
for i in delta_ps:         # iterate for delta_ps
    for j in range(len(Dimensionless_a)):     # iterate for dimensionless_a
        ax.plot(np.log10(T_d['%i'%i][j][:]),np.log10(R_d['%i'%i][j][:]**2 - 1), '.', c = color[col])
    col += 1
    
plt.title('Theoretical type curve figure:%i (Mpa)' %(delta_p/1000000.))
plt.xlabel('Dimensionless time (t)')
plt.ylabel('Dimensionless invasion radius $({r^2 - 1})$')
# plt.legend()
ax.annotate('*gp = 2.2 Mpa', xy = (7.8,4))
ax.annotate('*bp = 3.2 Mpa', xy = (11, 6))
ax.annotate('*rp = 4.2 Mpa', xy = (7.5,7.8))
# Good time_step algorithm.
# Radius_d  the fluid invasion radius.
plt.show()


'''
def Hough_Transform(x, y):
    X = np.array(x)
    Y = np.array(y)
    n = len(x)
    d = np.zeros(shape=(9,181))
    # generate matrix for the calculation
    for j in range(n):
        for i in range(181):
            theta = i * math.pi/180
            a = X[j] * math.cos(theta) + Y[j] * math.sin(theta)
            d[j,i] = a
            
    # data visulization for the calculation
    for j in range(n):
        plt.plot(range(181),d[j,:],'r')
    
    ax = plt.axes()
    plt.xlim([0, 180])
    ax.set_xlabel('theta')
    ax.set_ylabel('d')
    ax.annotate('Point',xy=(140,2))
    ax.xaxis.grid(True, which = 'major', linestyle = 'dotted')
    ax.yaxis.grid(True, which = 'major', linestyle = 'dotted')
    b = 'Polar Representation for Line'
    plt.title(b)
'''
data = [[math.log10(T_d['3200000'][0][i]) for i in range(R_d['3200000'][0].shape[0])],
        [math.log10(R_d['3200000'][0][i]**2 - 1) for i in range(R_d['3200000'][0].shape[0])]]



