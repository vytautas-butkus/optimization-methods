# -*- coding: utf-8 -*-
"""
Created on Wed Nov  9 16:03:50 2016

@author: butkus
"""

import numpy as np
import matplotlib.pyplot as plt

#number of iterations for averaging
N_random = 100

#standard deviation of random energy off-set [cm-1]
E_sigma = 73

H = []

with open('hamiltonian.txt') as data_file:
    data = data_file.readlines()
    
H_lt = [[float(value) for value in line.strip().split()] for line in data]

N = len(H_lt[-1])

H = np.array([line+[0]*(N-len(line)) for line in H_lt])
    
E_disorderless, _ = np.linalg.eigh (H)

print ("Eigenvalues:", E_disorderless)


def gaussian (x, mu=0, sigma=1.0):
    """
    Normalized Gaussian function
    """
    return np.exp(-np.power(x-mu,2.0) / (2 * np.power(sigma, 2))) / (sigma * np.sqrt(2.0 * np.pi))
    
def density_of_states (E, ensemble, delta):
    """
    Density of ensemble E values, obtained by using given delta function
    """
    return np.sum([delta(E-e) for e in ensemble])/len(ensemble)
    
E_ensemble = []

for _ in range(N_random):
    E_ensemble = np.concatenate((E_ensemble, np.linalg.eigh(H + np.diag(np.random.normal(0, E_sigma, N)))[0]))
    
E = np.linspace(14500.0, 15500, num=100)
# Probability density
E_pdensity = [probability_density(x, E_ensemble, lambda x: gaussian(x, 0, 10.0)) for x in E]
# Energy histogram
E_histogram = np.histogram(E_ensemble, 20, normed=True)
width = E_histogram[1][1]-E_histogram[1][0]


plt.figure('Energy distribution density')
plt.bar(E_histogram[1][:-1]-width/2.0, E_histogram[0], width=width*0.9, edgecolor='none')
plt.plot(E, E_pdensity, 'r-')
plt.xlabel (r'Energy $E$, cm$^{-1}$')
plt.ylabel (r'Density of states $P(E)$')
plt.legend(['density of states', 'histogram'])


