# -*- coding: utf-8 -*-
"""
Created on Fri Nov 30 16:44:03 2018

@author: sadra
"""

import numpy as np
from pypolycontain.lib.zonotope import zonotope,zonotope_directed_distance
from pypolycontain.visualization.visualize_2D import visualize_2D_zonotopes as visZ

G_1=np.array([[1,0,0,1,1],[0,1,0,-1,-3]])
G_2=np.array([[1,0,1,1,2,4,-4],[0,1,1,-1,3,3,1]])
x_1=np.array([0,1]).reshape(2,1)
x_2=np.array([1,0]).reshape(2,1)

z1=zonotope(x_1,G_1,color="red")
z2=zonotope(x_2,G_2,color="green")
visZ([z2,z1],title="Zonotopes")
D12=zonotope_directed_distance(z1,z2)
D21=zonotope_directed_distance(z2,z1) 