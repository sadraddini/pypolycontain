# -*- coding: utf-8 -*-
"""
Created on Thu Nov  1 09:40:19 2018

@author: sadra
"""

# External imports:
import numpy as np
import scipy as sp
from gurobipy import Model

from pypolycontain.lib.inclusion_encodings import subset_minkowski_right,subset_zonotope_both
from pypolycontain.lib.zonotope import zonotope,zonotope_inside
from pypolycontain.lib.elimination import project
from pypolycontain.visualization.visualize_2D import visualize_2D as vis
from pypolycontain.visualization.visualize_2D import visualize_2D_zonotopes as visZ

from pypolycontain.utils.utils import PI,valuation,null_space
from pypolycontain.utils.utils import vertices_cube as vcube

model=Model("Zonotope")
#G_l=np.array([[1,1,0,-1,0.5,-1],[0,1,0.5,-0.5,-1,2]])*1
#G_l=np.array([[1,0],[0,1]])*3.2
G_l=np.array([[1,0,0,1,1],[0,1,0,-1,-3]])*1
#G_r=np.array([[2,1,0,1,1,-1,-2,-1],[-1,3,0,1,-1,1,0,5],[0,0,1,1,-1,-1,2,1]])*1
G_r=np.array([[1,0,1,1,1,-2],[0,1,1,-1,3,2]])*1
x_l=np.array([0,1]).reshape(2,1)
x_r=np.array([1,0]).reshape(2,1)
(alpha,beta)=subset_zonotope_both(model,x_l,G_l,x_r,G_r)

#p_l=translate(project(G_l,PI(G_l.shape[1]),np.ones((2*G_l.shape[1],1))),x_l)
#p_r=translate(project(G_r,PI(G_r.shape[1]),np.ones((2*G_r.shape[1],1))),x_r)
#vis([p_r,p_l])
visZ([zonotope(x_r,G_r),zonotope(x_l,G_l)],list_of_dimensions=[0,1])
#visZ([zonotope(x_l,G_l),zonotope(x_r,G_r)],list_of_dimensions=[1,3])
#visZ([zonotope(x_l,G_l),zonotope(x_r,G_r)],list_of_dimensions=[2,3])
#visZ([zonotope(x_l,G_l),zonotope(x_r,G_r)],list_of_dimensions=[0,3])


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