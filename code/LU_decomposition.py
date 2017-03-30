# -*- coding: utf-8 -*-
"""
Created on Mon Nov 21 19:43:14 2016

@author: butkus
"""

def lu_factorization (A, method="Crout"):
    """LU factorization of matrix A using given Crout/Doolittle method"""
    if not (method in ["Crout", "Doolittle"]):
        print ("Wrong method provided: Use \"Crout\" or \"Doolittle\"")
        return (None, None)
        
    n = len(A)
    L = np.zeros ((n,n))
    U = np.zeros ((n,n))
    
    if method == "Crout":
        for i in range(n):
            L[i,i] = 1
            U[0,i] = A[0,i]

        for i in range(1,n):
            L[i,0] = A[i,0]/A[0,0]

        for i in range(1,n):
            for j in range(i,n):
                U[i,j] = A[i,j] - np.dot(L[i,:i], U[:i,j])

            for j in range(i+1,n):
                L[j,i] = (A[j,i] - np.dot(L[j,:i], U[:i,i]))/U[i][i]

    return (L, U)

def pivot (A):
    """Pivoting of matrix A"""
    
    n = len(A)
    P = np.diag([1]*n)
    
    for i in range(n):
        j =  np.argmax(np.abs(A[i:,i])) + i
#        print(j, A[i3:,i])
        P_new = np.diag([1]*n)
        P_new[i][i] = 0
        P_new[j][j] = 0
        P_new[i][j] = 1
        P_new[j][i] = 1
        A = np.dot(P_new,A)
        P = np.dot(P_new,P)
        
        
    return (P,A)
    

if __name__ == "__main__":
    from scipy.linalg import lu, inv
    import numpy as np
    
    A = np.array([[ 0, -4, -6,  2],
     [ 1, -4, -7,  7],
     [ 2, -6, -9, 11],
     [-1,  4,  7, -9]])
    
    
    P, B = pivot(A)
     
    L, U = lu_factorization (B, method="Crout")
