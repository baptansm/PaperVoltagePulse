#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 15 21:16:56 2020

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
# ----------------------------------------------------------------------------------
#KWANT SIMULATION    
def make_system_DC(a = 1, t = 1.0,L=50,W = 2, x0 = 20, y0 = -0.4,iP = 5,lx = 4, ly = 0.2,
                U = 5, V0 = 0,
               Vp = 0.005, tau = 100):
    lat = kwant.lattice.square(a=1, norbs = 1)
    syst = kwant.Builder()
 #   def onsite(site):
  #      (x,y) = site.pos
   #     Vx = (np.tanh((x - x0)/lx) + np.tanh(-(x - x0)/lx))*0.5*V0
    #    Vy = (np.tanh((y - y0)/ly) +  np.tanh(-(y + y0)/ly) +2)*0.5
      #  Vy =-(1-np.tanh((y-y0)/ly))
     #   return Vx*Vy + 4*t
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
    return syst
#--------------------------------------------------------
def get_kf(Ef):
    return np.arccos((4-Ef)/2)
#----------------------------------------
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
    
#-----------------------------------------------
def get_freq_amp(T,D,t_tar, imin = 0, imax = 150):
    
    W= 1
    N = np.argmin([abs(t-t_tar) for t in T])

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
    
   # ax.plot(X,c_moy)
    amp = (abs(max(c_moy)) + abs(min(c_moy)))*0.5
    #spectrum = np.fft.fft(c_moy)
    #freq = np.fft.fftfreq(len(spectrum), d =1)
    #ax1.plot(freq,abs(spectrum))
    
    mfreq = np.array(mfreq[imin:imax])
    x = np.linspace(imin, imax, 400)
    #f, ax = plt.subplots(1,1, figsize = (15,14))
    #ax.set_xlabel('position $x$', fontsize = 25)
    #ax.tick_params(axis= 'both', labelsize = 20)
    #ax.set_ylabel('density $n$',  fontsize = 25)
    #ax.plot(mfreq*max(c_moy)*0.8)
    
    freq = find_freq_sq(mfreq)
    
        
   # ax.plot(np.cos(2*np.pi*freq*X)*max(c_moy)*0.4)
        
        

    return 2*np.pi*freq, amp
#-------------------------------------------------------
  
def conductance(syst, energies):
    # Compute conductance
    data = []
    for energy in energies:
        smatrix = kwant.smatrix(syst, energy)
        data.append(smatrix.transmission(1, 0))
    return data

#-----------------------------------------------
def get_kf_freqn(file, t_tar=300, imax = 100, imin =0, QCP = 0):
    if QCP ==0:
        times , charge = pickle.load(open(file, 'rb'))
    if QCP == 1 :
        times , charge = pickle.load(open(file, 'rb'))
    #print(times)
    Ef = float(file.split('_')[-4].split('Ef')[-1])
   # nbar = float(file.split('_')[-1].split('.n')[0].split('n')[-1])
   # print('Ef =', Ef)
    nbar = 0.01
    Q , amp= get_freq_amp(times, charge,t_tar, imax = imax, imin=imin ) 
    kf = get_kf(Ef= Ef)
    return Q, kf, Ef, nbar, amp
#---------------------------------------------
def plot_amp_vs_G(files, QCP = 0, t_tar=300, imax = 100, imin = 0):
    Q = []
    KF = []
    EF = []
    NBAR = []
    AMP = []
    #fac = 4*np.pi*np.sqrt(np.log(2)/2)
    fac = 1
    #Ef = float(files[0].split('_')[-3].split('Ef')[-1])
    
    #ethr = [2.2,2.3,2.4,2.5,2.6]
 #   ethr = [0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9]
    for file in files:
        q,k,e,n, amp = get_kf_freqn(file, QCP = QCP, t_tar = t_tar, imax = imax, imin = imin)
      #  if e in ethr:
        Q.append(q)
        KF.append(k)
        EF.append(e)
        NBAR.append(n*fac)
        AMP.append(amp)
        print(e)
        print('Ef = ',e, "Q =", q, "kf = ", k, 'nbar =', n)
        
    f, (ax, ax1,ax2) = plt.subplots(3,1, figsize = (16,15), sharex =True)
   # ax.set_xlabel('$\overline{n}$', fontsize = 25)
    ax.set_xlim(2,6)
    ax.tick_params(axis= 'both', labelsize = 20)
    ax.set_ylabel('$q$',  fontsize = 25)
    ax.scatter(EF,Q, s = 200, marker = 's')
    ax1.scatter(EF,AMP, s = 100, marker = 's')
    ax1.set_ylim(0,0.4e-5)
  #  ax1.set_xlabel('$\overline{n}$', fontsize = 25)
    ax1.tick_params(axis= 'both', labelsize = 20)
    ax1.set_ylabel('$A$ [$a^{-1}$]',  fontsize = 25)
    ax.grid(1)
    ax1.grid(1)
    ax1.ticklabel_format(axis='y', style='sci', scilimits=(0,1))
   # f.suptitle('$E_F$ = {:1.1f}, $k_F$ = {:1.2f}'.format(EF[0], KF[0]), fontsize =25 )    #par = np.polyfit(KF,Q, deg = 1)
    syst = make_system_DC(W= 1, L = 3000 , V0 = 1.7, lx = 2.1).finalized()
    bands = kwant.physics.Bands(syst.leads[0])
    energies = np.linspace(0,6.1,300)    
    transm = conductance(syst, energies)
    transm_thr = conductance(syst, EF)
#    transm_thr2 = conductance(syst, ethr2)

    ax2.plot(energies, transm ,color = 'black')
    ax2.scatter(EF, transm_thr)
#    ax2.scatter(transm_thr2, ethr2, color = 'r')
    ax2.set_xlabel('$E_F$ [$\gamma$]', fontsize = 25)
    ax2.set_ylabel('$G$ [$e^2/h$]', fontsize = 25)
    ax2.tick_params(axis= 'both', labelsize = 20)
    plt.savefig('plot_G_AMP_Q.pdf', format = 'pdf', dpi = 80)
    f3, ax3 = plt.subplots(1,1, figsize = (8,6), sharex =True)
    ax3.scatter(transm_thr, AMP, s = 100)
    ax3.set_ylabel('$A$ [$a^{-1}$]',fontsize = 25)
    ax3.set_xlabel('$G$ [$e^2/h$]', fontsize = 25)
    ax3.tick_params(axis= 'both', labelsize = 20)
    plt.savefig('plot_A_vs_G.pdf', format = 'pdf', dpi = 80)
    return Q,KF,EF


files500= glob.glob('d*tau500*.npy')
plot_amp_vs_G(files500, QCP = 0, t_tar=3000, imax = 500, imin = 0)

