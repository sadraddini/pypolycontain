# -*- coding: utf-8 -*-
"""
Created on Thu Nov  1 09:40:19 2018

@author: sadra
"""

# External imports:
import numpy as np
import scipy as sp
from gurobipy import Model

from pypolycontain.lib.containment_encodings import subset_zonotopes
from pypolycontain.lib.zonotope import zonotope,zonotope_inside
from pypolycontain.visualization.visualize_2D import visualize_2D_zonotopes as visZ

from pypolycontain.utils.utils import PI,valuation,null_space
from pypolycontain.utils.utils import vertices_cube as vcube

model=Model("Zonotope")
#G_l=np.array([[1,0,0,1,1],[0,1,0,-1,-3]])*1
#G_r=np.array([[1,0,1,1,1,-2],[0,1,1,-1,3,2]])*1
G_l=np.array([[1,0,0,3],[0,1,2,-1]])*0.8
G_r=np.array([[1,0,1,1,2,-2],[0,1,1,-1,5,2]])*1
x_l=np.array([0,1]).reshape(2,1)
x_r=np.array([1,0]).reshape(2,1)
zono_l=zonotope(x_l,G_l)
zono_r=zonotope(x_r,G_r)
subset_zonotopes(model,zono_l,zono_r)
visZ([zono_r,zono_l],title="")
model.optimize()