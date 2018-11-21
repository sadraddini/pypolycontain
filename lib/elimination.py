#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 10 12:14:51 2018

@author: sadra
"""

import numpy as np
from scipy.linalg import block_diag

from pypolycontain.utils.redundancy_reduction import canonical_polytope
from polytope import polytope

      
def fourier_motzkin_eliminate_single(var_index,A,b,C=None,d=None,atol=10**-8):
    """
    Performs Fourier-Motzkin elimination method
    Inputs:
        var_index: an integer between [0,N], the index of variable to be removed, N is the number of variables
        A and b: A x <= b
        C and d (optional): Cx=d
        A,B,c,d should be numpy arrays
    Output:
        A_new and b_new such that A_new x_reduced <=b_new, where x_reduced does not include var_index
    """
    if var_index>A.shape[1]:
        raise("Error: %d is greather the number of variables. Choose a variable index between 0 and %d"%(var_index,A.shape[1]))
    if A.shape[0]!=b.shape[0]:
        raise("Error: number of rows in A: ",A.shape[0]," and b:",b.shape[0]," mismatch")
    if type(C)==type(np.array([1])):
        A=np.vstack((A,C,-C))
        b=np.vstack((b,d,-d))
        return fourier_motzkin_eliminate_single(var_index,A,b,None,None,atol)
    else:
        phi_positive=[i for i in range(A.shape[0]) if A[i,var_index]>=atol] # list of positive var entries
        phi_negative=[i for i in range(A.shape[0]) if A[i,var_index]<=-atol]  # list of negative var entries
        phi_core=[i for i in range(A.shape[0]) if abs(A[i,var_index])<atol]  # list of zero var entries
        s_smaller=np.diag(1/A[phi_positive,var_index]) # positive
        s_larger=np.diag(1/A[phi_negative,var_index]) # negative
        A_positive=np.dot(s_smaller,A[phi_positive,:]) # A of postives scaled by var entries
        b_positive=np.dot(s_smaller,b[phi_positive,:])
        A_negative=np.dot(s_larger,A[phi_negative,:])
        b_negative=np.dot(s_larger,b[phi_negative,:]) 
        """ We have A_positive x_other + x_r <= b_positive
        --> We have A_negative x_other + x_r >= b_negative
        --> We have b_postive - b_negative >= (A_neg - A _pos) * x_other (all combinations)
        """
        A_new=np.empty((0,A.shape[1]-1))
        b_new=np.empty((0,1))
        other=list(range(0,var_index))+list(range(var_index+1,A.shape[1]))
        for i in range(len(phi_positive)):
            for j in range(len(phi_negative)):
                alpha=(-A_negative[j,other]+A_positive[i,other]).reshape(1,len(other))
                beta=b_positive[i,:]-b_negative[j,:]
                A_new=np.vstack((A_new,alpha))
                b_new=np.vstack((b_new,beta))
        if phi_core!=[]:
            A_new=np.vstack((A_new,A[phi_core,:][:,other]))
            b_new=np.vstack((b_new,b[phi_core,:]))
        return canonical_polytope(A_new,b_new)

def project(T,A,b,C=None,d=None,atol=10**-8):
    """
    Finds the H-representation of T{Ax<=b, Cx=d}
    Inputs: T
    """
    (m,n)=T.shape # m: y, n: x, y=Tx
    if type(C)!=type(np.array([1])):
        A=np.hstack((np.zeros((A.shape[0],m)),A))
        b=b
        C=np.hstack((-np.eye(m),T))
        d=np.zeros((m,1))
        (A,b)=fourier_motzkin_eliminate_single(n+m-1,A,b,C,d,atol)
        for j in range(n-1):
            (A,b)=fourier_motzkin_eliminate_single(A.shape[1]-1,A,b,None,None,atol)
        return polytope(A,b)
    else:
        print("Projecting with A,b,C,d")
        A=np.hstack((np.zeros((A.shape[0],m)),A))
        b=b
        C=np.vstack((np.hstack((-np.eye(m),T)),np.hstack((np.zeros((C.shape[0],n)),C))))
        d=np.vstack((np.zeros((m,1)),d))
        (A,b)=fourier_motzkin_eliminate_single(A.shape[1]-1,A,b,C,d,atol)
        for j in range(n-1):
            (A,b)=fourier_motzkin_eliminate_single(A.shape[1]-1,A,b,None,None,atol)
        return polytope(A,b)
        
def convexhull(list_of_polytopes,atol=10**-8):
    """
    Computes a H-representation of the convex hull of a list of polytopes
    Inputs:
        list_of_polytopes: pairs of (H,h)
    """
    f=open("convexhull.log","w")
    n=list_of_polytopes[0].H.shape[1] # The space dimension
    N=len(list_of_polytopes)
    if N>2:
        p_2=convexhull(list_of_polytopes[0:2],atol)
        return convexhull([p_2]+list_of_polytopes[2:],atol)
    I_N=np.hstack(tuple([np.eye(n)]+([-np.eye(n)]*N)))
    X1=block_diag(np.zeros((2*n,n)),*[p.H for p in list_of_polytopes])
    X1[0:n,:]=I_N
    X1[n:2*n:]=-I_N
    X2=np.vstack( (np.zeros((2*n,N)),block_diag(*[-p.h for p in list_of_polytopes]), ))
    X3=np.vstack((-np.eye(N),np.ones((1,N))))
    # A and b
    A=np.block([[X1, X2], [np.zeros((X3.shape[0],X1.shape[1])), X3]])
    b=np.zeros((A.shape[0],1))
    b[A.shape[0]-1,0]=1
    # Fourier Motzkin Elimination Method
    for i in range(N*(n+1)):
        f.write("\nIteration %d from %d: A shape is (%d,%d) and b length is %d"%(i,N*(n+1),A.shape[0],A.shape[1],b.shape[0]))
        print "Iteration %d from %d: A shape is (%d,%d) and b length is %d"%(i,N*(n+1),A.shape[0],A.shape[1],b.shape[0])
        (A,b)=fourier_motzkin_eliminate_single(A.shape[1]-1,A,b,None,None,atol)
    f.close()
    return polytope(A,b)