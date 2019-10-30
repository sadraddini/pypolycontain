#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 30 13:35:37 2019

@author: sadra
"""
import numpy as np
import warnings

#pycdd    
try:
    from pypolycontain.lib.containment.complete import extreme_rays_for_containment as erc
    from pypolycontain.lib.objects import H_polytope,AH_polytope
except:
    warnings.warn("You don't have pyplycontain not properly installed.")
    
n=8
q=5  
H_y=np.vstack((np.eye(n),-np.eye(n)))
Y=np.random.random((q,n))-0.5
#Y=np.eye(2)
circombody=AH_polytope(t=np.zeros((Y.shape[0],1)),T=Y,P=H_polytope(H=H_y,h=np.zeros((Y.shape[0],1))))
Theta=erc(circombody,N=1)
print(Theta.shape)
print(np.all(Theta>-10**(-8)))