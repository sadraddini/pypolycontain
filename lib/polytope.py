# -*- coding: utf-8 -*-
"""
Created on Wed Oct 24 13:56:09 2018

@author: sadra
"""

import numpy as np

class polytope():
    def __init__(self,H,h):
        self.H=H
        self.h=h
        self.n=H.shape[1]
        if h.shape[1]!=1:
            ValueError("Error: not appropriate h size, it is",h.shape)
        if H.shape[0]!=h.shape[0]:
            ValueError("Error: not consistent dimension of H: %d and h: %d"%(H.shape[0],h.shape[0]))
    
    def __repr__(self):
        return ("polytope in R^%d"%self.n)
    
    def show(self):
        print(self)
        print("H=",self.H)
        print("h=",self.h)

def translate(p,d):
    """
    Given a polytope p, translate it by d 
    """
    return polytope(p.H,p.h+np.dot(p.H,d))