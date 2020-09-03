#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  3 11:40:45 2020

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



def plot_densities(file, W = 1,imin=0, imax = 2000) :
    
    
    times , charge_density , c= pickle.load(open(file, 'rb'))
    print(times)
  #Ef = float(file.split('_')[-2].split('Ef')[-1])
    f, ax = plt.subplots(1,1, figsize = (10,10))
    

    
   # tm = 1000

   
    L = int(len(charge_density[0])/W)
    X = [n for n in range(len(charge_density[0]))]
    X1 = [n for n in range(L)][imin:imax]
    N = len(charge_density[0])
    lines = []
    T = []
    i = 0
    ax.plot(X[imin:imax], charge_density[0][imin:imax])

    #ax.legend(lines, [l.get_label() for l in lines])
   # ax.set_yticks([])
    ax.tick_params(axis= 'both', labelsize = 18)
   # print(T)
    #ax.set_yticklabels(T)

    ax.set_xlabel('position $x$', fontsize = '20')
    #ax.set_ylabel('charge density $n(i)$ = $\sum_j n(i,j)$', fontsize = '20')
    ax.set_ylabel('$n(t=0)$ $[1/a]$', fontsize = '20')


files = glob.glob('*2.02*.npy')
print(files)
plot_densities(files[0], W = 1, imin = 1700,imax = 3500)
plt.savefig('plot_density_t0_Ef202.pdf', format ='pdf', dpi = 60)