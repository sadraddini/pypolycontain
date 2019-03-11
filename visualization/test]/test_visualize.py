#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 11:46:58 2019

@author: sadra
"""

import numpy as np

from pypolycontain.lib.zonotope import zonotope
from pypolycontain.visualization.visualize_2D import visualize_2D_zonotopes as visZ

G=np.random.random((2,3))
xbar=np.random.random((2,1))

N=10
Z={i:zonotope(np.random.random((2,1))*5,np.random.random((2,3)),color=np.random.random(3)) for i in range(N)}
visZ(Z.values())