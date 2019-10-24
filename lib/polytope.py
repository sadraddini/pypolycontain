# -*- coding: utf-8 -*-
"""
Created on Wed Oct 24 13:56:09 2018

@author: sadra
"""

import numpy as np
from pypolycontain.utils.utils import unique_rows


class polytope():
    def __init__(self,H,h):
        """
        Class Polytope. A Polytope is defined as (x \in R^n | Hx <= h)
        """
        if h.shape[1]!=1:
            ValueError("Error: not appropriate h size, it is",h.shape)
        if H.shape[0]!=h.shape[0]:
            ValueError("Error: not consistent dimension of H: %d and h: %d"%(H.shape[0],h.shape[0]))
        self.type="H-polytope"
        self.H,self.h=unique_rows(H,h)
        self.n=H.shape[1]
        self.hash_value = None

    def __repr__(self):
        return ("polytope in R^%d"%self.n)

    def __hash__(self):
        if self.hash_value is None:
            self.hash_value = hash(str(np.hstack([self.H, self.h])))
        return self.hash_value

    def show(self):
        print(self)
        print(("H=",self.H))
        print(("h=",self.h))
        
    def if_inside(self,x,tol=10**-5):
        if x.shape[0]!=self.H.shape[1]:
            return ValueError("H and x dimensions mismatch")
        return all(np.dot(self.H,x)<=self.h+tol)

    
def translate(p,d):
    """
    Given a polytope p, translate it by d 
    """
    return polytope(p.H,p.h+np.dot(p.H,d))

def Box(N,d=1,corners=None):
    """
    returns N-dimensional Box 
    corners=typle of 2 numpy arrays: lower_corner and upper_corner
    """
    H=np.vstack((np.eye(N),-np.eye(N)))
    if corners==None:
        h=d*np.ones((2*N,1))
    else:
        l,u=corners[0:2]
        if not all(u>=l):
            raise ValueError("Upper-right corner not totally ordering lower-left corner")
        h=np.vstack((u,-l))
    return polytope(H,h)