# -*- coding: utf-8 -*-
"""
Created on Thu Nov  1 09:40:19 2018

@author: sadra
"""

# External imports:
import numpy as np
import scipy as sp
from gurobipy import Model

from pypolycontain.lib.inclusion_encodings import subset_zonotope_both
from pypolycontain.lib.zonotope import zonotope,zonotope_inside
from pypolycontain.visualization.visualize_2D import visualize_2D_zonotopes as visZ

from pypolycontain.utils.utils import PI,valuation,null_space
from pypolycontain.utils.utils import vertices_cube as vcube

model=Model("Zonotope")
G_l=np.array([[1,0,0,1,1],[0,1,0,-1,-3]])*1
G_r=np.array([[1,0,1,1,1,-2],[0,1,1,-1,3,2]])*1
x_l=np.array([0,1]).reshape(2,1)
x_r=np.array([1,0]).reshape(2,1)
(alpha,beta)=subset_zonotope_both(model,x_l,G_l,x_r,G_r)


visZ([zonotope(x_r,G_r),zonotope(x_l,G_l)],list_of_dimensions=[0,1],title="")


# Let's look at vertices
zono_l=zonotope(x_l,G_l)
zono_r=zonotope(x_r,G_r)
y=zono_l.x.T+np.dot(zono_l.G,vcube(zono_l.G.shape[1]).T).T
Table=np.empty((y.shape[0],1))
for i in range(y.shape[0]):
    Table[i]=zonotope_inside(zono_r,y[i,:].reshape(2,1))
    print y[i,:],Table[i,0]

model.optimize()
alpha=valuation(alpha)
beta=valuation(beta)
print Table.T

N=G_r.shape[1]
C1=np.dot(np.linalg.pinv(PI(N).T),G_r.T)
D1=sp.linalg.orth(C1)
C2=(PI(N).T)
D2=null_space(C2)
y=np.hstack((D1,D2))
print np.linalg.matrix_rank(y)