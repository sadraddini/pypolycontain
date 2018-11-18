# -*- coding: utf-8 -*-
"""
Created on Tue Nov  6 10:24:55 2018

@author: sadra
"""
import numpy as np
from gurobipy import Model

from pyinpolytope.utilities.inclusion_encodings import subset_zonotope_both,constraints_AB_eq_CD


class zonotope():
    """
    Definition of a Zonotope
    """
    def __init__(self,x,G):
        self.x=x
        self.G=G
        
def zonotope_distance(z1,z2,eps_min=0,eps_max=10,eps=0.05):
    """
    Using binary search for finding Hausdorff distance between polytopes
    """
    d=(eps_max+eps_min)/2.0
    if (eps_max-eps_min)<eps:
        return d
    else:
        print("searching between",eps_min,"and",eps_max)
        model=Model("Zonotope Distance")
        subset_zonotope_both(model,z1.x,z1.G,z2.x,np.hstack((z2.G,d*np.eye(z2.G.shape[0]))))
        model.optimize()
        if model.Status==3:
            return zonotope_distance(z1,z2,d,eps_max,eps)
        else:
            return zonotope_distance(z1,z2,eps_min,d,eps)

def zonotope_order_reduction(z,N):
    if N<z.G.shape[1]:
        print("Zonotope already has less columns that what you want!")
        return        
    else:
        n=z.G.shape[0]
        model=Model("Zonotope_order_reduction")
        X=np.array((N,1),dtype='object')
        c=model.addVar(lb=1)
        subset_zonotope_both(model,np.zeros((n,1)),z.G,np.zeros((n,1)),X,c)
        pass

def zonotope_inside(z,x):
    model=Model("inside")
    n=z.x.shape[0]
    p=np.empty((z.G.shape[1],1),dtype='object')
    for row in range(p.shape[0]):
        p[row,0]=model.addVar(lb=-1,ub=1)
    model.update()
    constraints_AB_eq_CD(model,np.eye(n),x-z.x,z.G,p)
    model.setParam('OutputFlag', 0)
    model.optimize()
    if model.Status==3:
        return False
    else:
        return True   