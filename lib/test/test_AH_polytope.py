#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 20 15:24:54 2019

@author: sadra
"""
import numpy as np

from pypolycontain.lib.polytope import polytope,Box
from pypolycontain.lib.AH_polytope import AH_polytope

from time import time


n_p=40
n=30
P=Box(n_p)
T=(np.random.random((n,n_p))-0.5)*50
t=(np.random.random((n,1))-0.5)
X=AH_polytope(T,t,P)
x=np.zeros((n,1))

X.method="scipy"
start=time()
print X.is_inside(x)
print "scipy:",time()-start

X.method="Gurobi"
start=time()
print X.is_inside(x)
print "Gurobi:",time()-start