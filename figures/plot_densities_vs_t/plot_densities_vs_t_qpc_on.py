#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 11 15:24:24 2020

@author: baptiste
"""


import kwant
import tkwant
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erf
import kwant_spectrum
import functools
import pickle
import glob
from scipy.optimize import curve_fit

def plot_densities(file, W = 1, tmin =0, tmax = 500,imin =10, imax = 1400, fac = 1e7, M =1, title ='') :
    times , charge_density, current = pickle.load(open(file, 'rb'))

        
    #print(times)
  #Ef = float(file.split('_')[-2].split('Ef')[-1])
    f, ax = plt.subplots(1,1, figsize = (10,10))
    f.suptitle(title, fontsize = 20)

    L = int(len(charge_density[0])/W)
    X = [n for n in range(len(charge_density[0]))]
    X1 = [n for n in range(L)][imin:imax]
    N = len(charge_density[0])
    lines = []
    T = []
    i = 0
    for c, time in zip(charge_density, times):
        c = c - charge_density[0]
        c_reduced = []
        #c -= c[-1]
        if tmin< time < tmax:
            T.append(time)
            for n in range(L):
                m = [c[n*W+i] for i in range(0,W)]
                c_reduced.append(np.sum(m, axis = 0)/W)
               # c_reduced.append(c[n*W + 1])
            c_reduced = np.array(c_reduced[imin:imax])
            if i%M == 0:
           
                lines.append(ax.plot(X1, c_reduced*fac+ time, color = 'black',
                                 label = 't = {:2.0f}'.format(time))[0])
            #ax.axhline(i*offset, linestyle = '--', alpha = 0.3, color = 'black')
            i+=1

    ax.tick_params(axis= 'both', labelsize = 18)
    ax.set_xlabel('x', fontsize = '20')
    ax.set_ylabel('time $t$', fontsize = '20')

files = glob.glob('*V03*')
plot_densities(files[0], W =1, tmin =1000, tmax = 6200,imax = 1600, fac = 6e6, M =1, title ='') 
plt.savefig('plot_density_vs_t_qpc_on.pdf', format = 'pdf', dpi = 60)
  