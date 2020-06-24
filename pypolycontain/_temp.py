#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 23 23:35:22 2020

@author: sadra
"""

import numpy as np
from itertools import combinations
import pypolycontain as pp

G=np.random.random((2,8))-0.5
G=np.array([[1,0,0],[0,1,0]])
Z=pp.zonotope(G,color='red')


S=combinations(range(G.shape[1]),G.shape[0])
V=0
for s in S:
#    print(s)
    Gs=np.hstack([G[:,i:i+1] for i in s])
    print(Gs)
    V+=abs(np.linalg.det(Gs))
V*=2**G.shape[0]
print(V)

pp.visualize([Z,pp.zonotope(np.eye(2)*0.5)],title=r'volume=%0.002f'%V)

S=combinations(range(G.shape[1]),G.shape[0])
V_dot=np.zeros(G.shape)
for s in S:
    Gs=np.hstack([G[:,i:i+1] for i in s])
    D=np.linalg.det(Gs)
    if D!=0:
        adj_Gs=D*np.linalg.inv(Gs)
    else:
        print("Warning: we have a singular matrix")
        e=np.eye(G.shape[0])*10**(-5)
        adj_Gs=np.linalg.det(Gs+e)*np.linalg.inv(Gs+e)
    X=adj_Gs.T*np.sign(np.linalg.det(Gs+e))
    for i in range(len(s)):
        V_dot[:,s[i]]+=X[:,i]
print(V_dot)