# -*- coding: utf-8 -*-
"""
Created on Fri Mar 24 16:38:28 2017

@author: butkus
"""
import pandas as pd
import matplotlib.pyplot as plt

def exp_dec_ls (tau, t, y):
    
    A = [np.exp(-np.array(t)/tau_val) for tau_val in tau]
    
    
    

def exp_dec_multi_fit (:
    



# read whole matrix
fln_data = pd.read_csv('FLN.txt', decimal='.', sep='\t', header=1)

# set index as first column and remove it from the dataframe
fln_data.set_index('0', inplace=True)

#==============================================================================
# PLOT
#==============================================================================
plt.figure()
plt.clf()
plt.plot(fln_data)
plt.legend
plt.xlabel("Time, ps")
plt.ylabel("FLN value, a.u.")
#==============================================================================

n = 2
