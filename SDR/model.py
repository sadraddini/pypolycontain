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

def constraint_XY(X,Y,G,Z,sign="="):
    """
    Pose the constraint and variables for the following quadratic Constraint:
    X*Y + G*Z sign 0
    """
    constraints=[]
    A,S,C,W={},{},{},{}
    for row in range(X.shape[0]):
        for column in range(Y.shape[1]):
            x_row=X[row,:]
            x_row=cp.atoms.affine.reshape.reshape(x_row,(x_row.shape[0],1))
            y_column=Y[:,column]
            y_column=cp.atoms.affine.reshape.reshape(y_column,(y_column.shape[0],1))
            A[row,column],S[row,column],C[row,column],W[row,column]=constraint_xy(x_row,y_column)
            if sign=="=":
                constraints+=[cp.trace(A[row,column]*S[row,column]) + sum([G[row,k]*Z[k,column] for k in range(G.shape[1])]) == 0 ] + C[row,column] 
            elif sign==">":
                constraints+=[cp.trace(A[row,column]*S[row,column]) + sum([G[row,k]*Z[k,column] for k in range(G.shape[1])]) >= 0 ] + C[row,column]  
            elif sign=="<":
                constraints+=[cp.trace(A[row,column]*S[row,column]) + sum([G[row,k]*Z[k,column] for k in range(G.shape[1])]) >= 0 ] + C[row,column]  
            else:
                raise ValueError(sign," not identifed")
    return A,S,C,W,constraints

n=4   
X=cp.Variable((n,n))
Y=cp.Variable((n,n))
G=np.random.random((n,n))
Z=cp.Variable((n,n))

constraint_all=[]
A,S,C,W,constraints_XY=constraint_XY(X,Y,G,Z,sign="=")
constraint_all+=constraints_XY
prob = cp.Problem(cp.Minimize( cp.trace(G*X)),constraint_all)
prob.solve(verbose=True)
print(("The optimal value is", prob.value))
print("A solution X is")
print((X.value))

for row in range(X.shape[0]):
    for column in range(Y.shape[1]):
        print(row,column,np.trace(np.dot(A[row,column],S[row,column].value))+np.dot(G,Z.value)[row,column])

#n=3
#
## Define and solve the CVXPY problem.
## Create a symmetric matrix variable.
#x = cp.Variable((n,1))
## The operator >> denotes matrix inequality.
#
#G=np.eye(n)
#g=np.ones((n,1))
#
#constraints = [ x <= g, x >= g*0 ]
#
#A,S,constraint_xy,W=constraint_xy(x,x)
#constraints+=[ cp.trace(A*S) == 1] + constraint_xy 
#
#c=np.array([1,-1,1]).reshape(n,1)
#
#prob = cp.Problem(cp.Minimize(c.T *x),
#                  constraints)
#prob.solve(verbose=True)
#
## Print result.
#print("The optimal value is", prob.value)
#print("A solution X is")
#print(x.value)