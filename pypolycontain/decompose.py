"""
@author: kasra
Testing Boxing and PCA order reduction methods by illustration in 2D dimension
"""
import numpy as np
import matplotlib.pyplot as plt
import pypolycontain as pp
from cvxopt import solvers, matrix, spdiag, log

def acent(A, b):
    m, n = A.size
    def F(x=None, z=None):
        if x is None: return 0, matrix(1.0, (n,1))
        if min(x) <= 0.0: return None
        f = -sum(log(x))
        Df = -(x**-1).T
        if z is None: return f, Df
        H = spdiag(z[0] * x**-2)
        return f, Df, H
    return solvers.cp(F, A=A, b=b)['x']


A = matrix(np.random.rand(3,6))
b = matrix(np.random.rand(3))
acent(A,b)
