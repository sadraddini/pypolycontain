#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 11 17:49:54 2019

@author: sadra
"""

import cvxpy as cp
import numpy as np

class cvxpy_constraints:
    def __init__(self):
        pass
    
    def add_variable(self,x):
        pass
    
def constraint_xy(x,y):
    """
    x^T * y 
    output:
        
    """
    s=cp.atoms.affine.vstack.vstack((x,y))
    n=s.shape[0]
    S = cp.Variable((n,n), symmetric=True)
    W_up=cp.atoms.affine.hstack.hstack((S,s))
    W_down=cp.atoms.affine.hstack.hstack((s.T,np.array([[1]])))
    W=cp.atoms.affine.vstack.vstack((W_up,W_down))
    C = [W >> 0]
    I=np.eye(n/2)
    A=np.vstack((np.hstack((I*0,I/2)),np.hstack((I/2,I*0))))
    return A,S,C,W

def constraint_XY(X,Y,G,Z):
    """
    Pose the constraint and variables for X*Y=G*Z
    X*Y + G*Z sign 0
    """
    raise NotImplementedError
    
    

n=3

# Define and solve the CVXPY problem.
# Create a symmetric matrix variable.
x = cp.Variable((n,1))
# The operator >> denotes matrix inequality.

G=np.eye(n)
g=np.ones((n,1))

constraints = [ x <= g, x >= g*0 ]

A,S,constraint_xy,W=constraint_xy(x,x)
constraints+=[ cp.trace(A*S) == 1] + constraint_xy 

c=np.array([1,-1,1]).reshape(n,1)

prob = cp.Problem(cp.Minimize(c.T *x),
                  constraints)
prob.solve(verbose=True)

# Print result.
print("The optimal value is", prob.value)
print("A solution X is")
print(x.value)