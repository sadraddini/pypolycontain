#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 26 12:44:44 2019

@author: sadra
"""

import numpy as np

from pypolycontain.lib.AH_polytope import AH_polytope,to_AH_polytope
from pypolycontain.lib.zonotope import zonotope,zonotope_directed_distance
from pypolycontain.lib.polytope import Box
from pypolycontain.lib.hausdorff.hausdorff import Hausdorff_directed
from pypolycontain.visualization.visualize_2D import visualize_2D_zonotopes as visZ

#import time

G_l=np.array([[1,0,0,1,1],[0,1,0,-1,-3]])*1
G_r=np.array([[1,0,1,1,5,6,-9,-7],[0,1,1,-1,3,2,3,2]])*1
x_l=np.array([0,8]).reshape(2,1)
x_r=np.array([1,0]).reshape(2,1)

z1=zonotope(x_l,G_l)
z2=zonotope(x_r,G_r)
visZ([z2,z1],title="Zonotopes")
D12=zonotope_directed_distance(z1,z2)
Q1=to_AH_polytope(z1)
Q2=to_AH_polytope(z2)

#start=time.time()
#D12=zonotope_directed_distance(z1,z2)
#t1=time.time()-start
#
#start=time.time()
#D21=zonotope_directed_distance(z2,z1)
#t2=time.time()-start

#start=time.time()
D12_H=Hausdorff_directed(Q1,Q2)
#t3=time.time()-start

#start=time.time()
D21_H=Hausdorff_directed(Q2,Q1)
#t4=time.time()-start

#print D12,D12_H
#print D21,D21_H
#print t1,t2,t3,t4
print D12_H,D21_H

visZ([zonotope(z2.x,np.hstack((z2.G,D12_H*np.eye(z2.x.shape[0])))),z2,z1],title="")
visZ([zonotope(z1.x,np.hstack((z1.G,D21_H*np.eye(z1.x.shape[0])))),z1,z2],title="")