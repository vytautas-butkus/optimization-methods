# -*- coding: utf-8 -*-
"""
Created on Fri Nov 18 15:18:30 2016

@author: butkus
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from scipy.optimize import minimize, differential_evolution

def exp_dec (par):
    """calculates exponential decay vector, par[0] is decay constant tau"""
    tau = par[0]
    return A0 + A1*np.exp(-t/tau)
    
def deviation (par):
    """fit deviation from data"""
    sumsq = np.sum(np.power(df['rho11'].values-exp_dec(par), 2))
    return sumsq/len(df)
    
    
if __name__ == '__main__':
    filename = 'heom_allSB_2_10_2_dimer_de100_J40_l300_g100_300K.txt'
    
    df = pd.read_csv(filename, sep='\t', index_col=0)
    
    #Getting rid of columns, caused by \t\n at the end of line
    df = df[df.columns[~df.columns.str.contains('Unnamed:')]]
    
    t = np.array(df.index)
    t_min = t[0]
    t_max = t[-1]

    #stupid move: at infinity!=t_max exponential is zero
    A0 = df['rho11'][t_max]
    
    #smart move: initially, exponential is equal to 1
    A1 = df['rho11'][t_min] - A0
    
    fig1 = plt.figure('Density matrix')
    fig1.clf()

    #plotting a few guesses    
    ax1 = plt.subplot(2,1,1)
    ax = ax1.plot(df['rho11'], 'k.', label=r'$\rho_{11}$')
    ax1.set_title ('A few guesses')
    ax1.set_ylim ([0, 1])
    ax1.set_xlabel ('time, fs')
    ax1.set_ylabel ('population')
    
    for tau_init in [1000, 1200, 1500]:
        ax1.plot(t, exp_dec([tau_init]), label=r'$\tau=$'+str(tau_init)+' fs')
    
    ax1.legend()
    
       
    #fitting using BFGS method
    par_init = 100
    par_min = minimize (deviation, [par_init], method="BFGS")
    
    #ploting fitted curve
    ax2 = plt.subplot(2,1,2)
    
    ax2.set_title ('Fit')
    ax2.set_ylim ([0, 1])
    ax2.set_xlabel ('time, fs')
    ax2.set_ylabel ('population')
    
    
    ax2.plot(df['rho11'], 'k.', label=r'$\rho_{11}$')
    ax2.plot(t, exp_dec(par_min.x), label=r'$\tau=$'+str(round(par_min.x[0],0))+' fs (BFGS)')
    
    #fitting using conjugate gradient method
    par_min = minimize (deviation, [par_init], method="CG")
    ax2.plot(t, exp_dec(par_min.x), label=r'$\tau=$'+str(round(par_min.x[0],0))+' fs (CG)')
    
    #uneducated guess of initial point - using differential evolution method
    par_init_bad = 10000
    par_min = differential_evolution (deviation, [(-par_init_bad, par_init_bad)])
    
    print (r'With tau_init={:.2f}, tau_fitted={:.2f}'.format(par_init_bad, par_min.x[0]))
    
    ax2.plot(t, exp_dec(par_min.x), label=r'$\tau=$'+str(round(par_min.x[0],0))+' fs (DE)')
    
    ax2.legend()
    
    
