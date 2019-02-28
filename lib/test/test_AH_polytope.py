#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 20 15:24:54 2019

@author: sadra
"""
import numpy as np

from pypolycontain.lib.polytope import polytope,Box
from pypolycontain.lib.AH_polytope import AH_polytope,minimum_distance,check_collision,to_AH_polytope
from pypolycontain.lib.zonotope import zonotope
from pypolycontain.visualization.visualize_2D import visualize_2D_zonotopes as visZ

from time import time


n_p=40
n=30
P=Box(n_p)
T=(np.random.random((n,n_p))-0.5)*50
t=(np.random.random((n,1))-0.5)
X=AH_polytope(T,t,P)
x=np.zeros((n,1))

method="Gurobi"
start=time()
print X.is_inside(x)
print "Gurobi:",time()-start

#import time

G_l=np.array([[1,0,0,1,1],[0,1,15,-1,-3]])*1
G_r=np.array([[1,0,1,1,5,3,-1,-2],[0,1,1,-1,3,2,3,2]])*1
x_l=np.array([0,15]).reshape(2,1)
x_r=np.array([10,-10]).reshape(2,1)

z1=zonotope(x_l,G_l)
z2=zonotope(x_r,G_r)
visZ([z2,z1],title="Zonotopes")
Q1=to_AH_polytope(z1)
Q2=to_AH_polytope(z2)

D_min=minimum_distance(Q1,Q2)
print check_collision(Q1,Q2)