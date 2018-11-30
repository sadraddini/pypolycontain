# -*- coding: utf-8 -*-
"""
Created on Thu Nov 29 21:59:12 2018

@author: sadra
"""

# External imports:
import numpy as np

from pypolycontain.lib.zonotope import zonotope,zonotope_directed_distance
from pypolycontain.visualization.visualize_2D import visualize_2D_zonotopes as visZ

G_l=np.array([[1,0,0,1,1],[0,1,0,-1,-3]])*1
G_r=np.array([[1,0,1,1,1],[0,1,1,-1,3,]])*1
x_l=np.array([0,1]).reshape(2,1)
x_r=np.array([1,0]).reshape(2,1)

z1=zonotope(x_l,G_l)
z2=zonotope(x_r,G_r)
visZ([z2,z1],title="Zonotopes")
D12=zonotope_directed_distance(z1,z2)
D21=zonotope_directed_distance(z2,z1)

visZ([zonotope(z2.x,np.hstack((z2.G,D12*np.eye(z2.x.shape[0])))),z2,z1],title="")
visZ([zonotope(z1.x,np.hstack((z1.G,D21*np.eye(z1.x.shape[0])))),z1,z2],title="")