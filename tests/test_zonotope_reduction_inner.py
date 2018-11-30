# -*- coding: utf-8 -*-
"""
Created on Wed Nov 28 19:22:33 2018

@author: sadra
"""


import numpy as np

from pypolycontain.lib.zonotope import zonotope,zonotope_order_reduction_inner,zonotope_distance_cocenter
from pypolycontain.visualization.visualize_2D import visualize_2D_zonotopes as visZ

G=np.array([[1,0,1,1,1,-2,-3,3,2,-1,2,7],[0,1,1,-1,3,2,6,-1,-1,3,6,-1]])
x=np.array([0,1]).reshape(2,1)
Z=zonotope(x,G)

G_r=np.array([[1,0,1,2],[0,2,1,-1]])
(Z_r,Z_r_list)=zonotope_order_reduction_inner(Z,G_r)
for i in range(15):
    e=zonotope_distance_cocenter(Z,Z_r_list[i])
    fig=visZ([Z,Z,Z_r_list[i]],title=r"Iteration %d - Inner Approximation $d_H(\mathbb{Z}_r,\mathbb{Z})=%0.02f$"%(i,e),axis_limit=50)
    fig.savefig('figures/reduction_inner %d'%i,dpi=100)
