# -*- coding: utf-8 -*-
"""
Created on Fri Nov 16 08:04:10 2018

@author: sadra
"""


# External imports:
import numpy as np
import scipy as sp
from gurobipy import Model,GRB
import matplotlib.pyplot as plt

from pypolycontain.utils.utils import PI,valuation
from pypolycontain.lib.inclusion_encodings import subset_zonotope_both,constraints_AB_eq_CD
from pypolycontain.lib.zonotope import zonotope,zonotope_inside
from pypolycontain.utils.utils import vertices_cube as vcube

def zonotope_inside_scale(z,Y):
    """
    Maximum scaling of Zonotope such that contains all points in zonotope
    """
    model=Model("inside_scale")
    n,N=Y.shape
    p=np.empty((z.G.shape[1],N),dtype='object')
    scale=model.addVar(obj=1)
    for row in range(p.shape[0]):
        for column in range(N):
            p[row,column]=model.addVar(lb=-GRB.INFINITY,ub=GRB.INFINITY)
    model.update()
    for row in range(p.shape[0]):
        for column in range(N):
            model.addConstr(p[row,column]<=scale)
            model.addConstr(-p[row,column]<=scale)
    constraints_AB_eq_CD(model,np.eye(n),Y-z.x,z.G,p)
    model.setParam('OutputFlag', 0)
    model.optimize()
    return scale.X
    
def find_a_counterxample(n,tol=0.01):
    for N_r in range(n,13):
        for k in range(2):
            N_l=n+k
            max_iterations=100
            for i in range(max_iterations):
                model=Model("Zonotope")
                G_l=np.round(np.random.random((n,N_l))-0.5,1)*10.0
                G_r=np.round(np.random.random((n,N_r))-0.5,1)*10.0
                x_l=np.zeros((n,1))
                x_r=np.zeros((n,1))
                print G_l,G_r,np.linalg.matrix_rank(G_l),np.linalg.matrix_rank(G_r)
                G_l_var=np.empty((n,N_l),dtype='object')
                scale=model.addVar(obj=-1)
                for row in range(n):
                    for column in range(N_l):
                        G_l_var[row,column]=model.addVar(lb=-GRB.INFINITY,ub=GRB.INFINITY)
                model.update()
                for row in range(n):
                    for column in range(N_l):
                        model.addConstr(G_l_var[row,column]==G_l[row,column]*scale)
                (alpha,beta)=subset_zonotope_both(model,x_l,G_l_var,x_r,G_r)
                model.optimize()
                alpha=valuation(alpha)
                beta=valuation(beta)
                scale_theorem=scale.X
                print "scale_theorem is",scale_theorem
                if scale_theorem<=10**-5:
                    continue
                # Let's look at vertices
                zono_l=zonotope(x_l,G_l)
                zono_r=zonotope(x_r,G_r)
                y=zono_l.x.T+np.dot(zono_l.G,vcube(zono_l.G.shape[1]).T).T
                scale_vertices=1/zonotope_inside_scale(zono_r,y.T)        
                if (scale_vertices-scale_theorem)/(scale_vertices+10**-12)>tol:
                    print "find a counterexample!"
                    return (G_l,G_r,scale_vertices,scale_theorem)

(G_l,G_r,scale_vertices,scale_theorem)=find_a_counterxample(3,tol=10**-8)