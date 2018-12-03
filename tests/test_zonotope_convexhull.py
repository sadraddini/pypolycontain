# -*- coding: utf-8 -*-
"""
Created on Mon Dec  3 16:48:41 2018

@author: sadra
"""

# External imports:
import numpy as np
import scipy as sp
from gurobipy import Model,GRB,LinExpr

from pypolycontain.lib.containment_encodings import subset_zonotopes_convexhull,subset_zonotopes_disjunctive,add_Var_matrix
from pypolycontain.lib.zonotope import zonotope,zonotope_inside
from pypolycontain.visualization.visualize_2D import visualize_2D_zonotopes as visZ
from pypolycontain.visualization.visualize_2D import visualize_2D_zonotopes_convexhull as visZhull

from pypolycontain.utils.utils import PI,valuation,null_space
from pypolycontain.utils.utils import vertices_cube as vcube

model=Model("Zonotope")
#G_l=np.array([[1,0,0,1,1],[0,1,0,-1,-3]])*1
#G_r=np.array([[1,0,1,1,1,-2],[0,1,1,-1,3,2]])*1
G_l=np.array([[1,0,0,-3],[0,1,2,-1]])*1.5
G_r=np.array([[1,0,1,1,2,-2],[0,1,1,-1,5,2]])*1
x_l=np.array([0,10]).reshape(2,1)
x_r=np.array([-20,0]).reshape(2,1)
zono_l=zonotope(x_l,G_l)
zono_r=zonotope(x_r,G_r)

x=add_Var_matrix(model,(2,1))
G=add_Var_matrix(model,(2,2),delta=20)
(Lambda,beta,Gamma)=subset_zonotopes_convexhull(model,x,G,[zono_l,zono_r])
J=LinExpr()
J.add(G[0,0]+G[1,1])
model.addConstr(Lambda[zono_l]==0.5)
model.setObjective(J,GRB.MAXIMIZE)
model.optimize()
beta=valuation(beta)
Gamma=valuation(Gamma)
zono=zonotope(valuation(x),valuation(G)+np.eye(2)*0.01)


fig,ax=visZ([zono_r,zono_l,zono])
visZhull(fig,ax,[zono_r,zono_l],title="Convex Hull of Zonotopes")