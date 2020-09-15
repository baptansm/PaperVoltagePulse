#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 19 13:50:23 2020

@author: baptiste
"""

import kwant
import tkwant
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erf
import kwant_spectrum
import functools

import matplotlib.path as mpath
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection


def make_system_DC(a = 1, t = 1.0,L=6,W = 3, V0 = 0.5, lx = 5):
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
    time_dependent_hopping = (lat(2, 0), lat(3, 0))
    #tkwant.leads.add_voltage(syst, 0, gaussian)

    colormap = plt.cm.inferno
    
    def region_colors(site):
        if V0 ==0:
            return 'k'
        else:
            (x,y) = site.pos
            maxpot = V0
            rap = (onsite(site)- 4*t)/maxpot
            return colormap(rap*0.8)
    
    kwant.plot(syst, site_color= region_colors,lead_color='grey', 
                hop_color=lambda a, b: 'red' if (a, b) in [time_dependent_hopping] else 'k',
                hop_lw=lambda a, b: 0.3 if (a, b) in [time_dependent_hopping] else 0.1,
                fig_size = (20,15), colorbar = 0)
 #   lead_fin = lead.finalized()
    return syst
lat = kwant.lattice.square(a=1, norbs=1)
syst = kwant.Builder()
V0 = 1
L = 16
x0 = 12
t = 1
lx = 1
color_w = 'darkblue'
color_qpc = 'darkred'
syst[(lat(x, 0) for x in range(L))] = 4*t
syst[lat.neighbors()] = -t

def coupling(site1, site2, fac, omega, time):
    return -1 + fac * sin(omega * time)

# add the time-dependent coupling element
time_dependent_hopping = (lat(0, 0), lat(1, 0))
syst[time_dependent_hopping] = coupling
    #tkwant.leads.add_voltage(syst, 0, gaussian)

colormap = plt.cm.inferno
def onsite(site):
        (x,y) = site.pos
        return np.exp(-((x-x0)/lx)**2)*V0 + 4*t
#def region_colors(site):
#    if V0 ==0:
#        return 'k'
#    else:
#        (x,y) = site.pos
#        maxpot = V0
#        rap = (onsite(site)- 4*t)/maxpot
#        return colormap(rap*0.8)
def region_colors(site):
    (x,y) = site.pos
    if -2.2*lx + x0<x<2.2*lx + x0:
        return color_qpc
    else:
        return 'k'
 
lead = kwant.Builder(kwant.TranslationalSymmetry((-1, 0)))
lead[lat(0, 0)] = 1
lead[lat.neighbors()] = -1
syst.attach_lead(lead)
syst.attach_lead(lead.reversed())

# plot the system

kwant.plot(syst, site_color=region_colors, lead_color='grey',
           hop_lw=lambda a, b: 0.2 if (a, b) in [time_dependent_hopping] else 0.1,
           hop_color=lambda a, b: color_w if (a, b) in [time_dependent_hopping] else 'k',
          fig_size = (20,15));
plt.box(False)
plt.xticks([])
plt.yticks([])
#plt.text(x0*0.92, 2, '$QPC$', fontsize = 20)
#plt.text(-0.5, 2, '$w(t)$', fontsize = 20)
#
##plt.text(0.5, 0.5, 'Voltage drop', fontsize = 20)
#x = np.linspace(-0.7,0.7)
#y = np.exp(-((x+0.)/0.3)**2)*1.1 + 0.4
#plt.arrow(x0, 2.3, 0, -1.4, head_width=0.2, color = 'k')
#plt.arrow(1, 0.9, 0.8, 0, head_width=0.2, color = 'k')
#plt.plot(x,y, 'k')



plt.text(x0*0.9, 2, '$V_{QPC}(x)$', fontsize = 20, color = color_qpc)
plt.text(-0.5, 2, '$w(t)$', fontsize = 20, color = color_w)

#plt.text(0.5, 0.5, 'Voltage drop', fontsize = 20)
x = np.linspace(-0.7,0.7)
y = np.exp(-((x+0.)/0.3)**2)*1.1 + 0.4
#plt.arrow(x0, 2.3, 0, -1.4, head_width=0.2, color = 'k')
plt.arrow(1, 1, 0.8, 0, head_width=0.2, color = 'k')
plt.plot(x,y, color = color_w)
xqpc = np.linspace(-2*lx,2*lx,100) +x0
yqpc = np.exp(-((xqpc-x0)/lx)**2)*1.1 +0.4
plt.plot(xqpc,yqpc, color = color_qpc)