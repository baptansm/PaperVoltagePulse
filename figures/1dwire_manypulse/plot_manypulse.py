#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  5 21:21:04 2020

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



def plot_densities(file, W = 1, tmin =0, tmax = 500, imax = 2000, QCP = 1, fac = 1e7, M =1, title ='') :
    
    
    times , charge_density , c= pickle.load(open(file, 'rb'))
    print(times)
  #Ef = float(file.split('_')[-2].split('Ef')[-1])
    f, ax = plt.subplots(1,1, figsize = (10,10))
    f.suptitle(title, fontsize = 20)
    
    imin = 10

    
   # tm = 1000

   
    i = 0
    L = int(len(charge_density[0])/W)
    X = [n for n in range(len(charge_density[0]))]
    X1 = [n for n in range(L)][imin:imax]
    N = len(charge_density[0])
    lines = []
    T = []
    i = 0
    
    print('L = ', L)
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
    #ax.legend(lines, [l.get_label() for l in lines])
   # ax.set_yticks([])
    ax.tick_params(axis= 'both', labelsize = 18)
   # print(T)
    #ax.set_yticklabels(T)

    ax.set_xlabel('position $x$', fontsize = '20')
    #ax.set_ylabel('charge density $n(i)$ = $\sum_j n(i,j)$', fontsize = '20')
    ax.set_ylabel('time $t$', fontsize = '20')
    #ax.ticklabel_format(axis ='y', style = 'sci', scilimits= (0,10))
  #  syst = make_system_DC(W =W, L= L, V0 =V0,lx = 5).finalized()
   # plt.savefig('bands.pdf', format = 'pdf', dpi = 150)
  #  bands = kwant.physics.Bands(syst.leads[0])
  #  Velo = velocities(bands, Ef)
  #  lines  =[ ]
  #  colors = ['blue', 'orange', 'lime', 'red', 'purple']
 
   

          #      line_x.append(x)
           #     line_y.append(x/v) 
      #  lines.append(ax.plot(line_x, line_y, label = 'mode {:d}'.format(n), color = colors[n],
       #                      linewidth = 3,linestyle = '--',alpha = 1)[0])
files= glob.glob('d*V03*2.02*.npy')
plot_densities(files[0], W = 1, tmin =800, tmax = 10000000, QCP = 1, imax = 4000, 
               fac = 3e7, M = 6)
