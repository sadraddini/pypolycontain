#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May 30 10:28:04 2019

@author: sadra

We have three type of objects:
    * H-polytopes
    * Zonotopes
    * AH-polytopes
"""
# Numpy
import numpy as np
# Pypolycontain
from pypolycontain.utils.utils import unique_rows

class H_polytope():
    def __init__(self,H,h,symbolic=False,color="red"):
        """
        Class Polytope. A Polytope is defined as (x \in R^n | Hx <= h)
        """
        if h.shape[1]!=1:
            ValueError("Error: not appropriate h size, it is",h.shape)
        if H.shape[0]!=h.shape[0]:
            ValueError("Error: not consistent dimension of H: %d and h: %d"%(H.shape[0],h.shape[0]))
        self.H,self.h=H,h
        self.type="H_polytope"
        if not symbolic:
            self.H,self.h=unique_rows(H,h)
        self.n=H.shape[1]
        self.hash_value = None
        self.distance_program=None
#        self.color="red"

    def __repr__(self):
        return ("H_polytope in R^%d"%self.n)

    def __hash__(self):
        if self.hash_value is None:
            self.hash_value = hash(str(np.hstack([self.H, self.h])))
        return self.hash_value

    def if_inside(self,x,tol=10**-5):
        if x.shape[0]!=self.H.shape[1]:
            return ValueError("H and x dimensions mismatch")
        return all(np.dot(self.H,x)<=self.h+tol)
    
class zonotope():
    """
    Definition of a Zonotope
    """
    def __init__(self,x,G,name=None,color="green"):
        self.x=x
        self.G=G
        if name==None:
            self.name="zonotope"
        else:
            self.name=name
        self.color=color
        self.type="zonotope"
        self.hash_value = None
        self.distance_program=None
#        self.color="red"

    def __repr__(self):
        return self.name

    def __hash__(self):
        if self.hash_value is None:
            self.hash_value = hash(str(np.hstack([self.G, self.x])))  # FIXME: better hashing implementation
        return self.hash_value
    
class AH_polytope():
    """
    Affine Transformation of an H-polytope
    Attributes:
        P: The underlying H-polytope P:{x in R^q | Hx \le h}
        T: R^(n*q) matrix: linear transformation
        t: R^{n*1) vector: translation
    """
    def __init__(self,T,t,P,color='blue'):
        """
        Initilization: T,t,P. X=TP+t
        """
        self.T=T # Matrix n*n_p
        self.t=t # vector n*1
        self.P=P # Polytope in n_p dimensions
        self.n=T.shape[0]
        if T.shape[1]!=P.H.shape[1]:
            ValueError("Error: not appropriate T size, it is",T.shape[1],P.n)
        self.type="AH_polytope"
        self.method="Gurobi"
        self.hash_value = None
        self.distance_program=None
        self.vertices_2D=None
        self.color=color

    def __repr__(self):
        return "AH_polytope from R^%d to R^%d"%(self.P.n,self.n)

    def __hash__(self):
        if self.hash_value is None:
            self.hash_value = hash(self.P) + hash(str(np.hstack([self.T, self.t])))  # FIXME: better hashing implementation
        return self.hash_value
    
class hyperbox():
    def __init__(self,N=None,corners=[],d=1):
        """
        returns N-dimensional Box 
        corners=typle of 2 numpy arrays: lower_corner and upper_corner
        """
        if len(corners)==0:
            H=np.vstack((np.eye(N),-np.eye(N)))
            h=d*np.ones((2*N,1))
            self.l=-np.ones((N,1))*d
            self.u=-self.l
        else:
            l,u=corners[0:2]
            assert all(u>=l)
            N=l.shape[0]
            assert N==u.shape[0]
            H=np.vstack((np.eye(N),-np.eye(N)))
            if not all(u>=l):
                raise ValueError("Upper-right corner not totally ordering lower-left corner")
            h=np.vstack((u,-l))
            self.l=l
            self.u=u
        self.H_polytope=H_polytope(H,h)
        self.zonotope=zonotope((self.l+self.u)/2,np.diagflat((self.u-self.l)/2))
        self.n=N
    
def Box(N=None,d=1,corners=[]):
    """
    returns N-dimensional Box 
    corners=typle of 2 numpy arrays: lower_corner and upper_corner
    """
    if len(corners)==0:
        H=np.vstack((np.eye(N),-np.eye(N)))
        h=d*np.ones((2*N,1))
    else:
        l,u=corners[0:2]
        N=l.shape[0]
        assert N==u.shape[0]
        H=np.vstack((np.eye(N),-np.eye(N)))
        if not all(u>=l):
            raise ValueError("Upper-right corner not totally ordering lower-left corner")
        h=np.vstack((u,-l))
    return H_polytope(H,h)