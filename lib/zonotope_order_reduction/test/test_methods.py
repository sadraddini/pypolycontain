#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  7 17:05:28 2019

@author: sadra
"""

import numpy as np
# Pypolycontain
from pypolycontain.lib.objects import zonotope
from pypolycontain.visualization.visualize_2D import visualize_2D_zonotopes as visZ

from pypolycontain.lib.zonotope_order_reduction.methods import G_cut
q=10
np.random.seed(0)
Z=zonotope(np.zeros((2,1)),np.random.random((2,q))-0.5,color='red')
D,G=G_cut(Z,9,solver='gurobi')
Z_red=zonotope(np.zeros((2,1)),G,color='blue')
visZ([Z_red,Z],alpha=0.5)
