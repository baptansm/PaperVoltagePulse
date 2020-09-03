#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  6 13:48:50 2020

@author: baptiste
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 20 20:56:02 2020

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
#from scipy.fft import fft, fftshift, fftfreq

def get_freq_t0(T,D,t_tar, imin = 0, imax = 150):
    N = np.argmin([abs(t-t_tar) for t in T])
    L = len(D[0])
    X = np.array([n for n in range(L)])[imin:imax]
    c = D[0]
    mfreq =[]
    c = np.array(c[imin:imax])    
    offset = np.mean(c)
    c = c - offset
    print(c)
    for c_i in c:
        if c_i >= 0:
            mfreq.append(1)
        if c_i < 0:
            mfreq.append(-1)
   # mfreq = np.array(mfreq[imin:imax])
    mfreq = np.array(mfreq)
    x = np.linspace(imin, imax, 400)
    
  #  f, ax = plt.subplots(1,1, figsize = (15,14))
   # ax.set_xlabel('position $x$', fontsize = 25)
   # ax.tick_params(axis= 'both', labelsize = 20)
   # ax.set_ylabel('density $n$',  fontsize = 25)
    #ax.plot(X,c)
    #ax.plot(X,mfreq*max(c)*0.8) 
    freq = find_freq_sq(mfreq)
    return 2*np.pi*freq
    #ax.set_ylim(min(c_moy)*
def find_freq_sq(mfreq):
    ilist = []
    for i in range(1, len(mfreq)):
        if mfreq[i-1] < 0 and mfreq[i] >=0:
            ilist.append(i)
    ib = min(ilist)
    ie = max(ilist)
    lambd = abs(ie-ib)/len(ilist)
    print(ilist)
    return 1/lambd

#KWANT SIMULATION    
def make_system_DC(a = 1, t = 1.0,L=50,W = 2, x0 = 20, y0 = -0.4,iP = 5,lx = 4, ly = 0.2,
                U = 5, V0 = 0,
               Vp = 0.005, tau = 100):
    lat = kwant.lattice.square(a=1, norbs = 1)
    syst = kwant.Builder()
    x0 = L/2
    def onsite(site):
        (x,y) = site.pos
        return np.exp(-((x-x0)/lx)**2)*V0 + 4*t

    syst[(lat(x,y) for x in range(0,L) for y in range(0,W))] = onsite
    syst[lat.neighbors()] = -t  
    lead = kwant.Builder(kwant.TranslationalSymmetry((-a,0)))
    lead[(lat(0,y) for y in range(0,W))] = 4*t
    lead[lat.neighbors()] = -t
    syst.attach_lead(lead)
    syst.attach_lead(lead.reversed())
    colormap = plt.cm.plasma 
    
    def region_colors(site):
        (x,y) = site.pos
        maxpot =  0.8
        rap = (onsite(site)-4*t)/maxpot
        return colormap(rap)
    #kwant.plot(syst, site_color= region_colors, fig_size = (20,15), colorbar = 0)
 #   plt.savefig('syst_V0{:d}_W{:d}.pdf'.format(int(V0),W), format = 'pdf', dpi =60)
   # plt.savefig('syst.pdf', format = 'pdf', dpi =100)
    lead_fin = lead.finalized()
    return syst

def velocities(bands, Ef):
    N = len(bands(0))
    Velo = []
    for n in range(N):
        if bands(0)[n] < Ef:
            E0 = bands(np.pi/2)[n]
            A = E0 - bands(0)[n]
            kf = np.arccos((E0- Ef)/A)
            v = A*np.sin(kf)
        else:
            v = np.nan 
        Velo.append(v)
    return Velo

def get_kf(Ef = 2,n = 0, W= 3, L = 5 , V0 = 0,lx = 5):
    syst = make_system_DC(W =W, L= L, V0 =V0,lx = lx).finalized()
   # kwant.plotter.bands(syst.leads[0], show=False)
    bands = kwant.physics.Bands(syst.leads[0])
    N = len(bands(0))
    if bands(0)[n] < Ef:
        E0 = bands(np.pi/2)[n]
        A = E0 - bands(0)[n]
        kf = np.arccos((E0- Ef)/A)
        v = A*np.sin(kf)
    else:
        kf = np.nan 
    return kf

def get_kf_freq(file, t_tar=300,imin = 0, imax = 100, n=0):
    times , charge, c = pickle.load(open(file, 'rb'))
    #print(times)
    Ef = float(file.split('_')[-3].split('Ef')[-1])
   # print('Ef =', Ef)
    Q = get_freq_t0(times, charge,t_tar,imin = imin, imax = imax ) 
    kf = get_kf(W=1,Ef= Ef, n =n)
    return Q, kf, Ef

files = glob.glob('*.npy')

#get_kf_freq(files[0], t_tar=0,imin = 2000, imax = 2150, n=0)
KF = []
Q = []
EF = []
for file in files:
     q, kf, Ef = get_kf_freq(file, t_tar=0,imin = 500, imax = 2150, n=0)
     KF.append(kf)
     EF.append(Ef)
     Q.append(q)
     
f, ax = plt.subplots(1,1,figsize = (10,6))
ax.scatter(KF,Q)
ax.plot([min(KF),max(KF)],[min(KF)*2,max(KF)*2], color = 'k', label = '$q=2k_F$')
ax.legend(fontsize = 20)
ax.set_xlabel('$k_F$', fontsize =  '22')
ax.set_ylabel('$q$', fontsize =  '22')
ax.tick_params(axis= 'both', labelsize = 20)
plt.savefig('q_vs_kf_at_t0.pdf', format ='pdf', dpi = 60)


