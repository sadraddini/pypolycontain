# -*- coding: utf-8 -*-
"""
Created on Sun Nov  4 13:40:49 2018

@author: sadra
"""

# External imports:
import numpy as np
from gurobipy import Model,GRB,LinExpr,QuadExpr

from pypolycontain.lib.inclusion_encodings import subset_minkowski_right
from pypolycontain.lib.polytope import translate,polytope
from pypolycontain.lib.elimination import project
from pypolycontain.visualization.visualize_2D import visualize_2D as vis
from pypolycontain.utils.utils import PI,valuation

model=Model("Zonotope")
H=np.array([[1,1],[-1,1],[0,-1]])
h=np.array([[1,1,0]]).reshape(3,1)
p_triangle=polytope(H,h)

e=0.1
scale=0.68
H=np.array([[1,0],[-1,0],[0,-1],[1,1],[-1,1],[0,1]])
h=np.array([[1+e,1+e,1,1+e,1+e,1]]).reshape(6,1)
p_sum=polytope(H,h)
vis([p_sum],0.5)

H=np.array([[1,0],[-1,0],[0,1],[0,-1]])
h=np.array([[e,e,0,1]]).reshape(4,1)
p_line=polytope(H,h)
vis([p_triangle,p_line],0.5)

#H=np.array([[1,1],[-2,1],[0,1],[0,-1]])
#h=np.array([[0.3,0.5,0.5,0.7]]).reshape(4,1)
p_test=polytope(p_sum.H,p_sum.h*scale)

model=Model("Minkowski")
(t,T)=subset_minkowski_right(model,np.zeros((2,1)),np.eye(2),p_test,[p_triangle,p_line])
model.optimize()
t_triangle=valuation(t[p_triangle])
t_line=valuation(t[p_line])
T_triangle=valuation(T[p_triangle])
T_line=valuation(T[p_line])

p_sub_triangle=translate(project(T_triangle,p_triangle.H,p_triangle.h*scale),t_triangle)

p_sub_line=translate(project(T_line,p_line.H,p_line.h*scale),t_line)

vis([p_sum,p_sum,p_sum],0.5)

vis([p_triangle,p_line,p_sub_triangle,p_sub_line],0.5)