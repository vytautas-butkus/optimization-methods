# -*- coding: utf-8 -*-
"""
Created on Wed Nov  9 16:03:50 2016

@author: butkus
"""

import numpy as np
import matplotlib.pyplot as plt

def read_file_lt (filename):
    """reads lower-triangular matrix from file"""
    with open(filename) as data_file:
        data = data_file.readlines()
        
    H_lt = [[float(value) for value in line.strip().split()] for line in data]
    
    N = len(H_lt[-1])
    
    H = np.array([line+[0]*(N-len(line)) for line in H_lt])
    return H

def gaussian (x, mu=0, sigma=1.0):
    """Normalized Gaussian function"""
    return np.exp(-np.power(x-mu,2.0) / (2 * np.power(sigma, 2))) / (sigma * np.sqrt(2.0 * np.pi))
    
def density_of_states (E, ensemble, delta):
    """Density of ensemble E values, obtained by using given delta function"""
    return np.sum([delta(E-e) for e in ensemble])/len(ensemble)
    
def energies_hamiltonian_disorder (H, iterations = 500, E_sigma = 73):
    """Density of states of given Hamiltonian and Gaussian disorder"""
    N = len(H[-1])    
    
    E_disorderless, _ = np.linalg.eigh (H)
    
    E = []
    
    for _ in range(iterations):
        E = np.concatenate((E, np.linalg.eigh(H + np.diag(np.random.normal(0, E_sigma, N)))[0]))
    
    return E

if __name__ == '__main__':
    
    H = read_file_lt ("hamiltonian.txt")
    
    E_ensemble = energies_hamiltonian_disorder (H);
        
    E = np.linspace(14500.0, 15500, num=100)
    E_pdensity = [density_of_states(x, E_ensemble, lambda x: gaussian(x, 0, 40.0)) for x in E]
    
    
    plt.figure('Energy distribution density')
    plt.plot(E, E_pdensity, 'r-')
    plt.xlabel (r'Energy $E$, cm$^{-1}$')
    plt.ylabel (r'Density of states $P(E)$')


