# -*- coding: utf-8 -*-
"""
Created on Sat Jan  5 17:38:14 2019

@author: sadra
"""

import numpy as np

from pypolycontain.lib.orthogonal_projection import orthogonal_projection_fixed_Hx,orthogonal_projection_gradient_decent,gradient_decent
from pypolycontain.visualization.visualize_2D import visualize_2D as vis2
from pypolycontain.lib.polytope import polytope

H=np.array([[1,0],[0,1],[-1,0],[0,-1],[0,0],[0,0]])
F=np.array([[1],[0],[-1],[0],[1],[-1]])
g=np.array([[1,1,1,1,1,1]]).T

Hx_correct=np.array([[0.5,0],[0,1],[-0.5,0],[0,-1]])
#Hx=np.array([[0.5,0],[0,1],[-0.5,0],[0,-1]])*3
Hx=np.array([[0,-1],[1,1],[-1,1],[1,-1],[1,0]])*2
xbar=np.zeros((2,1))
#xbar=-np.array([1.8,0.8]).reshape(2,1)

delta=0.1
output=gradient_decent(Hx,H,F,g,xbar,delta,N=50)
import matplotlib.pyplot as plt
plt.plot([x[1] for x in output])
Hx=output[-1][0]

p_out=polytope(Hx_correct,np.ones((Hx_correct.shape[0],1)))
for i in range(len(output)):
    Hx=output[i][0]
    p_in=polytope(Hx,np.ones((Hx.shape[0],1)))
    fig=vis2([p_out,p_in],title=r"Orthogonal Projection %d $d_H=%f$"%(i,output[i][1]))
    fig.savefig('figures/orthogonal_projection %d'%i,dpi=100)
    
