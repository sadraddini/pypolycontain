#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May 30 15:49:45 2019

@author: sadra
"""
import numpy as np
# Pydrake
import pydrake.solvers.mathematicalprogram as MP
# Pypolycontain
from pypolycontain.lib.objects import H_polytope,zonotope,AH_polytope
from pypolycontain.lib.operations import Box,point_membership,directed_Hausdorff_distance
from pypolycontain.visualization.visualize_2D import visualize_2D_zonotopes as visZ

from pypolycontain.lib.hausdorff.hausdorff import Hausdorff_directed 
from time import time
#np.random.seed(0)
def test_memebership():
    # Test 1: # Random
    N,n,m=10,4,3
    H=np.random.random((N,n))-0.5
    h=np.random.random((N,1))+5
    T=np.random.random((m,n))
    t=np.random.random((m,1))*0
    H_P=H_polytope(H,h)
    P=AH_polytope(T,t,H_P)
    x=np.random.random((m,1))*0
    print point_membership(P,x,solver="gurobi")
    
    # Test 2: # Zonotope
    n=4
    P=zonotope(np.zeros((n,1)),np.eye(n),Box(n))
    x=np.random.random((n,1))
    print point_membership(P,x,solver="gurobi")
    
def test_hausdorff():
    n=2
    q1=12
    q2=12
    z1=zonotope(np.random.random((n,1)),np.random.random((n,q1))-0.5,color='red')
    z2=zonotope(np.random.random((n,1)),np.random.random((n,q2))-0.5,color='blue')
    start=time()
    D=directed_Hausdorff_distance(z1,z2,solver="gurobi")
    print "Mathematical Program:",D,"\t",time()-start
    start=time()
    print "Gurobipy:",Hausdorff_directed(z1,z2),"\t",time()-start
    z3=zonotope(z2.x,np.hstack((z2.G,D*np.eye(n))),color='green')
    visZ([z3,z1,z2],a=0.5,alpha=0.2)
    return D
    
def __main__():
    test_memebership()
    test_hausdorff()

__main__()