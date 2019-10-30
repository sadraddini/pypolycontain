#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Jun  8 14:39:26 2019

@author: sadra
"""

import numpy as np
from pypolycontain.lib.operations import to_AH_polytope

def subset(program,Q1,Q2):
    """
    Adds containment property Q1 subset Q2
    Inputs:
        Q1,Q2: either polytope, zonotope, or AH_polytope
    Output:
        No direct output, adds Q1 \subset Q2 to the model
    """
    Q1=to_AH_polytope(Q1)
    Q2=to_AH_polytope(Q2)
    Hx,Hy,hx,hy,X,Y,xbar,ybar=Q1.P.H,Q2.P.H,Q1.P.h,Q2.P.h,Q1.T,Q2.T,Q1.t,Q2.t
    qx,qy,nx,ny=Hx.shape[0],Hy.shape[0],X.shape[1],Y.shape[1]
    Lambda=program.NewContinuousVariables(qy,qx,'Lambda')
    Gamma=program.NewContinuousVariables(ny,nx,'Gamma')
    beta=program.NewContinuousVariables(ny,1,'beta')
    # Constraints
    program.AddBoundingBoxConstraint(0,np.inf,Lambda) # Lambda Non-Negative
    program.AddLinearConstraint(np.equal(X,np.dot(Y,Gamma),dtype='object').flatten()) #X=YGamma
    program.AddLinearConstraint(np.equal(ybar-xbar,np.dot(Y,beta),dtype='object').flatten()) 
    program.AddLinearConstraint(np.equal(np.dot(Lambda,Hx),np.dot(Hy,Gamma),dtype='object').flatten()) 
    program.AddLinearConstraint(np.less_equal(np.dot(Lambda,hx),hy+np.dot(Hy,beta),dtype='object').flatten())
    
def subset_soft(program,Q1,Q2,M=1000):
    """
    Adds containment property Q1 subset Q2
    Inputs:
        Q1,Q2: either polytope, zonotope, or AH_polytope
    Output:
        No direct output, adds Q1 \subset Q2 to the model
    """
    Q1=to_AH_polytope(Q1)
    Q2=to_AH_polytope(Q2)
    Hx,Hy,hx,hy,X,Y,xbar,ybar=Q1.P.H,Q2.P.H,Q1.P.h,Q2.P.h,Q1.T,Q2.T,Q1.t,Q2.t
    qx,qy,nx,ny=Hx.shape[0],Hy.shape[0],X.shape[1],Y.shape[1]
    Lambda=program.NewContinuousVariables(qy,qx,'Lambda')
    Gamma=program.NewContinuousVariables(ny,nx,'Gamma')
    beta=program.NewContinuousVariables(ny,1,'beta')
    D=program.NewContinuousVariables(ny,1,'D')
    # Constraints
    program.AddBoundingBoxConstraint(0,np.inf,Lambda) # Lambda Non-Negative
    program.AddLinearConstraint(np.equal(X,np.dot(Y,Gamma),dtype='object').flatten()) #X=YGamma
    program.AddLinearConstraint(np.equal(ybar-xbar,np.dot(Y,beta),dtype='object').flatten()) 
    program.AddLinearConstraint(np.equal(np.dot(Lambda,Hx),np.dot(Hy,Gamma),dtype='object').flatten()) 
    program.AddLinearConstraint(np.less_equal(np.dot(Lambda,hx),np.dot(hy,np.diag(D.reshape(ny)))+np.dot(Hy,beta),dtype='object').flatten())
    # Constraint and cost for D
    program.AddBoundingBoxConstraint(0,np.inf,D) # D Non-Negative
    program.AddLinearCost(sum([D[i,0]*M for i in range(ny)]))
    return D