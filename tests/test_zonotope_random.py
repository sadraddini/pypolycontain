# -*- coding: utf-8 -*-
"""
Created on Thu Nov 15 11:26:30 2018

@author: sadra
"""

# External imports:
import numpy as np
import scipy as sp
from gurobipy import Model,GRB
import matplotlib.pyplot as plt

from pypolycontain.lib.polytope import translate
from pypolycontain.lib.elimination import project
from pypolycontain.lib.inclusion_encodings import subset_zonotope_both,constraints_AB_eq_CD
from pypolycontain.lib.zonotope import zonotope,zonotope_inside
from pypolycontain.visualization.visualize_2D import visualize_2D as vis
from pypolycontain.visualization.visualize_2D import visualize_2D_zonotopes as visZ
from pypolycontain.utils.utils import PI,valuation,null_space
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

if True:    
    scale_table={}
    difference_table={}
    for n in range(3,11):
        for N_r in range(n,13):
            for k in range(2):
                N_l=n+k
                max_iterations=100
                scale_table[n,N_r,k]=np.empty((max_iterations,2))
                difference_table[n,N_r,k]=np.empty(max_iterations)
                for i in range(max_iterations):
                    model=Model("Zonotope")
                    G_l=np.random.random((n,N_l))-0.5
                    G_r=np.random.random((n,N_r))-0.5
                    x_l=np.zeros((n,1))
                    x_r=np.zeros((n,1))
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
                    # Let's look at vertices
                    zono_l=zonotope(x_l,G_l)
                    zono_r=zonotope(x_r,G_r)
                    y=zono_l.x.T+np.dot(zono_l.G,vcube(zono_l.G.shape[1]).T).T
                    scale_vertices=1/zonotope_inside_scale(zono_r,y.T)
            
                    print n,N_r,i,scale_vertices,scale_theorem
            
                    scale_table[n,N_r,k][i,0]=scale_vertices
                    scale_table[n,N_r,k][i,1]=scale_theorem
                    difference_table[n,N_r,k][i]=(scale_vertices-scale_theorem)/(scale_vertices+10**-12)
    #        plt.hist(difference_table[n,N_r])
    x_all=np.hstack([difference_table[key] for key in difference_table.keys()])
#plt.hist(x_all)
#plt.hist(x_all,30)
#plt.yscale('log')
#plt.xlabel(r'Loss')
#plt.ylabel(r'Number of Zonotopes')

from matplotlib import pyplot


for n in [3,5,8,10]:
    z_all=np.hstack([difference_table[n,N_r,k] for N_r in range(n,13)])
    pyplot.hist(z_all*1.0,10,alpha=0.6,label=r"$n=%d$"%n,normed=False)
    pyplot.yscale('log')
    #plt.show()
    pyplot.xlabel(r'Loss')
    pyplot.ylabel(r'Number of Zonotopes')    
#    pyplot.title(r'Zonotopes in $\mathbb{R}^{%d}$' %n)
pyplot.legend(loc='upper right')
pyplot.show()
