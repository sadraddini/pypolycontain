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
from pypolycontain.lib.operations import Box,point_membership

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
    print point_membership(P,x)
    
    # Test 2: # Zonotope
    n=4
    P=zonotope(np.zeros((n,1)),np.eye(n),Box(n))
    x=np.random.random((n,1))
    print point_membership(P,x)
    
def __main__():
    test_memebership()

__main__()