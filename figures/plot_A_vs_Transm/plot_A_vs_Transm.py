#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 11 16:14:44 2020

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

def get_freq2(T,D,t_tar, imin = 0, imax = 150):
    
    W= 1
    N = np.argmin([abs(t-t_tar) for t in T])
    f, ax = plt.subplots(1,1, figsize = (15,14))
    ax.set_xlabel('position $x$', fontsize = 25)
    ax.tick_params(axis= 'both', labelsize = 20)
    ax.set_ylabel('density $n$',  fontsize = 25)
    L = int(len(D[0])/W)
    X = np.array([n for n in range(L)])[imin:imax]
    c = D[N] - D[0]
    c_moy = []
    mp = []
    mm = []
    mfreq = []
    for n in range(L-1):
        m = [c[n*W+i] for i in range(0,W)] 
        c_i = np.sum(m, axis = 0)/W
        c_moy.append(np.sum(m, axis = 0)/W)
        if c_i >= 0:
            mfreq.append(1)
        if c_i < 0:
            mfreq.append(-1)
    c_moy = c
    c_moy = np.array(c_moy[imin:imax])
    
    ax.plot(X,c_moy)
    
    #spectrum = np.fft.fft(c_moy)
    #freq = np.fft.fftfreq(len(spectrum), d =1)
    #ax1.plot(freq,abs(spectrum))
    
    mfreq = np.array(mfreq[imin:imax])
    x = np.linspace(imin, imax, 400)
    ax.plot(mfreq*max(c_moy)*0.8)
    
    freq = find_freq_sq(mfreq)
    
        
   # ax.plot(np.cos(2*np.pi*freq*X)*max(c_moy)*0.4)
        
        

    return 2*np.pi*freq
    #ax.set_ylim(min(c_moy)*0.6, max(c_moy)*0.4)
  #  popt, pcov = fit(X,c_moy, p0 = [0, 1.1, 5e-6] )

#---------------------------------------------

def get_kf_freq(file, t_tar=300,imin = 0, imax = 100, n=0):
    times , charge = pickle.load(open(file, 'rb'))
    #print(times)
  #  Ef = float(file.split('_')[-3].split('Ef')[-1])
    Ef = float(file.split('_')[-4].split('Ef')[-1])
   # print('Ef =', Ef)
    Q = get_freq2(times, charge,t_tar,imin = imin, imax = imax ) 
    kf = get_kf(W=1,Ef= Ef, n =n)
    return Q, kf, Ef
#--------------------------------------------
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
    
#--------------------------------------------------------------
def get_freq_amp(T,D,t_tar, imin = 0, imax = 150):
    
    W= 1
    N = np.argmin([abs(t-t_tar) for t in T])
   # f, ax = plt.subplots(1,1, figsize = (15,14))
   # ax.set_xlabel('position $x$', fontsize = 25)
   # ax.tick_params(axis= 'both', labelsize = 20)
   # ax.set_ylabel('density $n$',  fontsize = 25)
    L = int(len(D[0])/W)
    X = np.array([n for n in range(L)])[imin:imax]
    c = D[N] - D[0]
    c_moy = []
    mp = []
    mm = []
    mfreq = []
    for n in range(L-1):
        m = [c[n*W+i] for i in range(0,W)] 
        c_i = np.sum(m, axis = 0)/W
        c_moy.append(np.sum(m, axis = 0)/W)
        if c_i >= 0:
            mfreq.append(1)
        if c_i < 0:
            mfreq.append(-1)
        
    c_moy = np.array(c_moy[imin:imax])
    
  #  ax.plot(X,c_moy)
    amp = (abs(max(c_moy)) + abs(min(c_moy)))*0.5
    #spectrum = np.fft.fft(c_moy)
    #freq = np.fft.fftfreq(len(spectrum), d =1)
    #ax1.plot(freq,abs(spectrum))
    
    mfreq = np.array(mfreq[imin:imax])
    x = np.linspace(imin, imax, 400)
 #   ax.plot(mfreq*max(c_moy)*0.8)
    
    freq = find_freq_sq(mfreq)
    
        
 #   ax.plot(np.cos(2*np.pi*freq*X)*max(c_moy)*0.4)
        
        

    return 2*np.pi*freq, amp
    #ax.set_ylim(min(c_moy)*0.6, max(c_moy)*0.4)
  #  popt, pcov = fit(X,c_moy, p0 = [0, 1.1, 5e-6] )
#--------------------------------------------------------
def get_kf(Ef):
    return np.arccos((4-Ef)/2)
#------------------------------------------------------------
def get_kf_freqn(file, t_tar=300, imax = 100, QCP = 0):
    if QCP ==0:
        times , charge = pickle.load(open(file, 'rb'))
    if QCP == 1 :
        times , charge = pickle.load(open(file, 'rb'))
    #print(times)
    Ef = float(file.split('_')[-4].split('Ef')[-1])
    nbar = float(file.split('_')[-2].split('.n')[0].split('n')[-1])
   # print('Ef =', Ef)
    Q , amp= get_freq_amp(times, charge,t_tar, imax = imax ) 
    kf = get_kf(Ef= Ef)
    return Q, kf, Ef, nbar, amp
#---------------------------------------------------------
def abs_sinus(n,A):
    return A*abs(np.sin(np.pi*n))

def fit_sinus(NBAR,AMP):
    popt, pcov = curve_fit(abs_sinus, NBAR,AMP)
    return popt,pcov
#------------------------------------------------------
def plot_q_vs_nbar_T(files, QCP = 0, t_tar=300, imax = 100):
    Q0 = []
    KF0 = []
    EF0 = []
    NBAR0 = []
    AMP0 = []
    Q024 = []
    KF024 = []
    EF024 = []
    NBAR024 = []
    AMP024 = []
    Q04 = []
    KF04 = []
    EF04 = []
    NBAR04 = []
    AMP04 = []

    #ethr = [2.2,2.3,2.4,2.5,2.6]
 #   ethr = [0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9]
    for file in files:
        v0 =  float(file.split('_')[3].split('V0')[-1])
        print('v=', v0)
        if v0 == 0.2:
            q,k,e,n, amp = get_kf_freqn(file, QCP = QCP, t_tar = t_tar, imax = imax)
            Q024.append(q)
            KF024.append(k)
            EF024.append(e)
            NBAR024.append(n)
            AMP024.append(amp)
        if v0 == 0.4:
            q,k,e,n, amp = get_kf_freqn(file, QCP = QCP, t_tar = t_tar, imax = imax)
            Q04.append(q)
            KF04.append(k)
            EF04.append(e)
            NBAR04.append(n)
            AMP04.append(amp)
        if v0 == 3:
            q,k,e,n, amp = get_kf_freqn(file, QCP = QCP, t_tar = t_tar, imax = imax)
            Q0.append(q)
            KF0.append(k)
            EF0.append(e)
            NBAR0.append(n)
            AMP0.append(amp)

    AMP0 = [e for _,e in sorted(zip(NBAR0,AMP0))]
    NBAR0.sort()
    AMP0= [e for _,e in sorted(zip(NBAR0,AMP0))]
    AMP0 = np.array(AMP0)
    AMP024 = [e for _,e in sorted(zip(NBAR024,AMP024))]
    NBAR024.sort()
    AMP024= [e for _,e in sorted(zip(NBAR024,AMP024))]
    AMP024 = np.array(AMP024)
    AMP04 = [e for _,e in sorted(zip(NBAR04,AMP04))]
    NBAR04.sort()
    AMP04= [e for _,e in sorted(zip(NBAR04,AMP04))]
    AMP04 = np.array(AMP04)
   
    popt0,pcov0 = fit_sinus(NBAR0,AMP0)
    popt024,pcov024 = fit_sinus(NBAR024,AMP024)
    popt04,pcov04 = fit_sinus(NBAR04,AMP04)
    
    
    x = np.linspace(0,3,100)
    f1, ax1 = plt.subplots(1,1, figsize = (8,6), sharex = True)
    colormap = plt.cm.Blues

    a = ax1.scatter(NBAR0,AMP0*1e4, s=100, marker = 's', color = colormap(0.99), 
                label = '$T = 0$')
    ax1.plot(x,abs_sinus(x,popt0[0]*1e4),color = colormap(0.99))
    
    c = ax1.scatter(NBAR04,AMP04*1e4, s=100, marker = 's', color = colormap(0.75),
                label = '$T = 0.67$')
    
    ax1.plot(x,abs_sinus(x,popt04[0]*1e4),color = colormap(0.75))

    b = ax1.scatter(NBAR024,AMP024*1e4, s=100, marker = 's',color =colormap(0.4),
                label = '$T = 0.37$')
    ax1.plot(x,abs_sinus(x,popt024[0]*1e4),color = colormap(0.4))

 

    ax1.set_xlabel('$\overline{n}$', fontsize = 25)
    ax1.tick_params(axis= 'both', labelsize = 20)
    ax1.legend(fontsize = 15)
    ax1.set_ylabel('$A$ [$10^{-4}. a^{-1}$]',  fontsize = 22)
  #  ax.grid(1)
    ax1.grid(1)
    #f.suptitle('$E_F$ = {:1.1f}, $k_F$ = {:1.2f}'.format(EF[0], KF[0]), fontsize =25 )
    ax1.ticklabel_format(axis='y', style='sci', scilimits=(0,1))
    plt.savefig('plot_A_vs_nbar_vs_T.pdf', format = 'pdf', dpi = 70)
    
    f2, ax2 = plt.subplots(1,1, figsize = (8,6), sharex = True)
    
    trans = [0,0.37,0.67]
    amps = [popt0[0],popt04[0],popt024[0]]
    ax2.scatter(trans,amps, s = 100)
    ax2.set_ylabel('$A_{fit}$',fontsize = 22)
    ax2.set_xlabel('$T$',fontsize = 22)
    ax2.grid(1)
    ax2.ticklabel_format(axis='y', style='sci', scilimits=(0,1))
    ax2.tick_params(axis= 'both', labelsize = 20)
    ax2.legend(fontsize = 15)
    plt.savefig('plot_A_from_fit_vs_Transm.pdf', format = 'pdf', dpi = 70)
#-------------------------------------------------------------------------------
files = glob.glob('*npy')
plot_q_vs_nbar_T(files, QCP = 1,t_tar=3000, imax = 700)
#plt.savefig('A_vs_nbar_vs_transm.pdf', format = 'pdf', dpi = 80)