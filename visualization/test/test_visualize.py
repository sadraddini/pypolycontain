#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 11:46:58 2019

@author: sadra
"""

import numpy as np

from pypolycontain.lib.objects import zonotope, AH_polytope
from pypolycontain.lib.operations import to_AH_polytope, convex_hull_of_point_and_polytope
from pypolycontain.visualization.visualize_2D import visualize_2D_zonotopes as visZ
from pypolycontain.visualization.visualize_2D import visualize_ND_AH_polytope as visND
from pypolycontain.visualization.visualize_2D import visualize_2D_AH_polytope as vis2D
import matplotlib.pyplot as plt

G=np.random.random((2,3))
xbar=np.random.random((2,1))

N=10
Z={i:zonotope(np.random.random((2,1))*5,np.random.random((2,3)),color=np.random.random(3)) for i in range(N)}
# visZ(list(Z.values()))

# Test visualize 2D AH-polytope

N=10
zono_list = []
ah_list = []
for i in range(N):
    random_x = np.random.random((7,1))*5
    random_G = np.random.random((7,5))
    random_color = np.random.random(3)
    zono_list.append(zonotope(random_x[0:2, :], random_G[0:2, :], color=random_color))
    ah_list.append(to_AH_polytope(zonotope(random_x,random_G, color=random_color)))
fig, ax = visND(ah_list, 0, 1, N=300)
fig2, ax2 = visZ(zono_list)
plt.show()