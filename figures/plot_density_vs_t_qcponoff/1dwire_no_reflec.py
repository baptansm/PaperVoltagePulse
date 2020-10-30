
#QCP FOR MPI
import warnings
import numpy as np
#import matplotlib.pyplot as plt
from scipy.special import erf
import time as temps
import pickle
import kwant
import tkwant
import kwant_spectrum
import functools

t0 = temps.time()

W = 1
L = 2500
V0 = 0
Ef = 2.02
#Energies = [2.01,2.02,2.03,2.04,2.05,2.06,2.07,2.08]
u=0
nbar = 0.01
tau = 500
#prm_str = 'W{:d}_L{:d}_V0{:1.1f}_Ef{:1.2f}_u{:1.1f}_n{:1.1f}_tau{:3.1f}'.format(W,L,V0,Ef,u, nbar,tau)
path = './Plots_MPI/'
#warnings.simplefilter('ignore') 
comm = tkwant.mpi.get_communicator()

def am_master():
    return comm.rank == 0 # cache possibly performance warning

# -------------------------------------------------------------------------------------
#Make system for TKWANT SIMULATION

def make_system(a = 1, t = 1.0,L=8,W = 1, lx = 50, V0 = 4):
    lat = kwant.lattice.square(a=1, norbs = 1)
    syst = kwant.Builder()


   # def onsite(site):
    #    (x,y) = site.pos
    #    Vx = (np.tanh((x - x0)/lx) + np.tanh(-(x + x0)/lx))*0.5*V0
   #     Vy = (np.tanh((y - y0)/ly) +  np.tanh(-(y + y0)/ly) +2)*0.5
     #   Vy =-(1-np.tanh((y-y0)/ly))
      #  return Vx*Vy + 4*t
    x0 = L*0.8
    def onsite(site):
        (x,y) = site.pos
        return np.exp(-((x-x0)/lx)**2)*V0 + 4*t
      # area under voltage pulse
    
    time0 = 1200
    Vp = (nbar / tau) * np.sqrt(2/np.pi)
    
    #Vp = 2*(nbar / tau) * np.sqrt(2*np.pi)
#    Vp = 4 * (nbar / tau) * np.sqrt(np.log(2) * np.pi)

    def gaussian(time):
        A = Vp*tau*np.pi*np.sqrt(np.pi/2)
        sigma = tau/np.sqrt(2)
        phi = A*(1 + erf((time - time0)/sigma))
        return phi
        
    

    syst[(lat(x,y) for x in range(0,L) for y in range(0,W))] = onsite
    syst[lat.neighbors()] = -t  

    
    
    lead = kwant.Builder(kwant.TranslationalSymmetry((-a,0)))
    lead[(lat(0,y) for y in range(0,W))] = 4*t
    lead[lat.neighbors()] = -t
    syst.attach_lead(lead)
    syst.attach_lead(lead.reversed())
    added_sites = tkwant.leads.add_voltage(syst, 0, gaussian)
    

    return syst, added_sites


def main(W = 1,L = 200,V0 = 0, Ef = 4):

    # create system
    syst, added_sites = make_system(W=W,L=L,V0 = V0)
    syst = syst.finalized()
    
    interaction_wavefunc = tkwant.manybody.Interaction(added_sites, u=u, update_error_tol=1E-6)
    syst = tkwant.manybody.CustomSystem(syst, interaction=interaction_wavefunc)
    # plot the system and dispersion
    chemical_potential = Ef
    tmax = 15000
    #times = np.concatenate((np.array([0]), np.linspace(800,6000,10)) ,axis = 0)
    times = np.linspace(0,tmax,28)
    # define an observable
   
    
    density_operator = kwant.operator.Density(syst)

    # do the actual tkwant simulation
   # print('occupation')
    spectra = kwant_spectrum.spectra(syst.leads)
    occupation = tkwant.manybody.lead_occupation(chemical_potential)
    #intervals = tkwant.manybody.calc_intervals(spectra, occupation)
    Interval = functools.partial(tkwant.manybody.Interval, order =25,quadrature = 'gausslegendre' )
    intervals = tkwant.manybody.calc_intervals(spectra, occupation, interval_type=Interval)
    intervals = tkwant.manybody.split_intervals(intervals, number_subintervals = 120)
    tasks = tkwant.manybody.calc_tasks(intervals, spectra, occupation)
    emin, emax = tkwant.manybody.calc_energy_cutoffs(occupation)
    boundaries = tkwant.leads.automatic_boundary(spectra, tmax=tmax, emin=emin, emax=emax)
    psi_init = tkwant.manybody.calc_initial_state(syst, tasks, boundaries)
    solver = tkwant.manybody.WaveFunction(psi_init, tasks)
    interaction_wavefunc.set_solver(solver, syst)

    #print('solver fini')
 

    charge_density = []
    #current = []
    for time in times:
        

        interaction_wavefunc.evolve(time=time)
        d= interaction_wavefunc.evaluate(density_operator)
     #   i= interaction_wavefunc.evaluate(current_operator)

        charge_density.append(d)
   #     current.append(i)
        if am_master():
            print("t =", time)
            print('max n - n0 =', max(abs(d - charge_density[0])))



    #charge_density = np.array(charge_density)

    return times, charge_density


	
times, charge_density = main(W = W,L = L ,V0 = V0, Ef = Ef) 
prm_str = 'W{:d}_L{:d}_V0{:1.1f}_Ef{:1.2f}_u{:1.1f}_n{:1.1f}_tau{:3.1f}'.format(W,L,V0,Ef,u, nbar, tau)
if am_master():
    
    results = (times, charge_density)
    pickle.dump(results, open('dens_current'+prm_str+'.npy', "wb"))
