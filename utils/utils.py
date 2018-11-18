# -*- coding: utf-8 -*-
"""
Created on Thu Nov  8 19:15:11 2018

@author: sadra
"""
from numpy import eye,vstack,ones,array,empty
import numpy as np

def valuation(x):
    """
    Description: given a set of Gurobi variables, output a similar object with values
    Input:
        x: dictionary or a vector, each val an numpy array, each entry a Gurobi variable
        output: x_n: dictionary with the same key as, each val an numpy array, each entry a float 
    """
    if type(x)==type(dict()):
        x_n={}
        for key,val in x.items():
            x_n[key]=ones(val.shape)
            (n_r,n_c)=val.shape
            for row in range(n_r):
                for column in range(n_c):
                    x_n[key][row,column]=x[key][row,column].X   
        return x_n
    elif type(x)==type(array([0])):
        x_n=empty(x.shape)
        for row in range(x.shape[0]):
            for column in range(x.shape[1]):
                x_n[row,column]=x[row,column].X
        return x_n
    else:
        raise("x is neither a dictionary or a numpy array")

def PI(n):
    return vstack((eye(n),-eye(n)))

def vertices_cube(T):
    """
    Description: 2**n * n array of vectors of vertices in unit cube in R^n
    """
    from itertools import product 
    v=list(product(*zip([-1]*T,[1]*T)))
    return array(v)
    

def null_space(A, rcond=None):
    """
    Construct an orthonormal basis for the null space of A using SVD
    Parameters
    ----------
    A : (M, N) array_like
        Input array
    rcond : float, optional
        Relative condition number. Singular values ``s`` smaller than
        ``rcond * max(s)`` are considered zero.
        Default: floating point eps * max(M,N).
    Returns
    -------
    Z : (N, K) ndarray
        Orthonormal basis for the null space of A.
        K = dimension of effective null space, as determined by rcond
    See also
    --------
    svd : Singular value decomposition of a matrix
    orth : Matrix range
    Examples
    --------
    One-dimensional null space:
    >>> from scipy.linalg import null_space
    >>> A = np.array([[1, 1], [1, 1]])
    >>> ns = null_space(A)
    >>> ns * np.sign(ns[0,0])  # Remove the sign ambiguity of the vector
    array([[ 0.70710678],
           [-0.70710678]])
    Two-dimensional null space:
    >>> B = np.random.rand(3, 5)
    >>> Z = null_space(B)
    >>> Z.shape
    (5, 2)
    >>> np.allclose(B.dot(Z), 0)
    True
    The basis vectors are orthonormal (up to rounding error):
    >>> Z.T.dot(Z)
    array([[  1.00000000e+00,   6.92087741e-17],
           [  6.92087741e-17,   1.00000000e+00]])
    """
    u, s, vh = np.linalg.svd(A, full_matrices=True)
    M, N = u.shape[0], vh.shape[1]
    if rcond is None:
        rcond = np.finfo(s.dtype).eps * max(M, N)
    tol = np.amax(s) * rcond
    num = np.sum(s > tol, dtype=int)
    Q = vh[num:,:].T.conj()
    return Q