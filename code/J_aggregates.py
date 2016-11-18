# -*- coding: utf-8 -*-
"""
Created on Thu Nov 10 18:53:33 2016

@author: butkus
"""
import numpy as np
import matplotlib.pyplot as plt

k_BT = 0.69503476   #Boltzman constant in cm-1/K

def coupling_nn (J, m, n):
    """coupling only between nearest neighbors"""
    return J if np.fabs(m-n)==1 else 0
    
def coupling_full (J, m, n):
    """coupling between all molecules"""
    return J / np.power(np.fabs(m-n), 3.0) if m-n != 0 else 0

def J_hamiltonian (N, E, J_mn):
    H = np.zeros((N,N))
    
    for i in range(N):
        for j in range(i+1,N):
            H[i][j] = J_mn(i,j)
    
    return H + H.T + np.diag([E]*N)
    
def make_rho_eq_SB (E, U, temperature):
    """Equilibrium density matrix in eigenstate basis"""
    N = len(E)
    
    rho = np.exp(-E/k_BT/temperature)
    rho = np.diag(rho / sum(rho))
    
    Um1 = U.T.conj()
    
    return np.dot(np.dot(Um1, rho), U)


if __name__ == '__main__':
    from hamiltonian_eigenvalues import energies_hamiltonian_disorder, density_of_states, gaussian
    
    N = 50
    E = 15000
    J = 100
    iterations = 10
    temperature = 300
    E_sigma = 20
    
    H_linear_coupling = J_hamiltonian (N, E, lambda x, y: coupling_full(J, x, y))
    H_nn_coupling = J_hamiltonian (N, E, lambda x, y: coupling_nn(J, x, y))
    

    fig = plt.figure('J-aggregate')
    fig.clf();
    ax1 = fig.add_subplot(2,1,1)
    ax2 = [fig.add_subplot(2,2,3), fig.add_subplot(2,2,4)]
    ax1.set_xlabel (r'Energy $E$, cm$^{-1}$')
    ax1.set_ylabel (r'Density of states $P(E)$')
    
    colors = ['red','blue']
    
    for index, H in enumerate([H_linear_coupling,  H_nn_coupling]):
        E_disorderless, _ = np.linalg.eigh (H)
        
        E_ensemble = []
        rho_eq_EB = np.zeros((N,N));
    
        for _ in range(iterations):
            E, U = np.linalg.eig(H + np.diag(np.random.normal(0, E_sigma, N)))
            
            E_ensemble = np.concatenate((E_ensemble, E))
            
            rho_eq_EB = rho_eq_EB + make_rho_eq_SB(E, U, temperature)
            
    
        #density of states
        E = np.linspace(14700.0, 15300.0, num=200)
        E_pdensity = [density_of_states(x, E_ensemble, lambda x: gaussian(x, 0, 5.0)) for x in E]
        
        #histogram
        hist, bin_edges = np.histogram (E_ensemble, bins=50, range=(14700.0,15300.0), normed=True)
        
        #plot density of states, histogram
        ax1.plot(E, E_pdensity, color=colors[index])
        width=0.5*(bin_edges[1]-bin_edges[0])
        ax1.bar(index*width + bin_edges[:-1], hist, width=width, color=colors[index])
        
        
        #calculate rho_disorder
        pl = ax2[index].pcolor(np.fabs(rho_eq_EB), cmap='RdBu', vmin=-0.1, vmax=0.1)
        ax2[index].axis('square')
        ax2[index].set_xlabel(r'Site no.')
        ax2[index].set_ylabel(r'Site no.')
    
    fig.colorbar(pl)
    ax1.legend(['full coupling, J=100', 'nn coupling, J=100'])
    for i in range(2):
        ax2[i].set_title(r'$\rho_\mathrm{EB}^\mathrm{eq}$ '+ ['linear coupling, J=100', 'nn coupling, J=100'][i])
    
    