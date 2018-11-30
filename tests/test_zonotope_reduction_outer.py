# -*- coding: utf-8 -*-
"""
Created on Tue Nov 20 19:43:51 2018

@author: sadra
"""

import numpy as np

from pypolycontain.lib.zonotope import zonotope,zonotope_order_reduction_outer,zonotope_distance_cocenter
from pypolycontain.visualization.visualize_2D import visualize_2D_zonotopes as visZ

G=np.array([[1,0,1,1,1,-2,-3,3,2,-1,2,7],[0,1,1,-1,3,2,6,-1,-1,3,6,-1]])
x=np.array([0,1]).reshape(2,1)
Z=zonotope(x,G)

G_r=np.array([[1,0,1,2],[0,2,1,-1]])
Z_r,Z_r_list,eps=zonotope_order_reduction_outer(Z,G_r,delta=0.01,scale=1.06)
i_max=len(eps)
for i in range(i_max):
    print "iteration",i,
    e=zonotope_distance_cocenter(Z_r_list[i],Z)
    fig=visZ([Z_r_list[i],Z],title=r"Iteration %d - Outer Approximation $d_H(\mathbb{Z}_r,\mathbb{Z})=%0.02f$"%(i,e),axis_limit=50)
    fig.savefig('figures/reduction_outer %d'%i,dpi=100)

fig.close()
import matplotlib.pyplot as plt
plt.plot(eps[:i])
