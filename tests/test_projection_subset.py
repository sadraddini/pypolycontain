# -*- coding: utf-8 -*-
"""
Created on Wed Nov 14 09:05:36 2018

@author: sadra
"""
import numpy as np
from gurobipy import Model

from pypolycontain.lib.polytope import polytope,translate
from pypolycontain.lib.elimination import project
from pypolycontain.visualization.visualize_2D import visualize_2D as vis

from pypolycontain.utils.utils import PI,valuation
from pypolycontain.lib.inclusion_encodings import subset_both_projection

model=Model("Projection")
H=np.array([[1,1],[-1,1],[0,-1]])
h=np.array([[1,1,0]]).reshape(3,1)
p_triangle=polytope(H,h)

H=np.array([[1,0],[-1,0],[0,-1],[1,1],[-1,1]])
h=np.array([[1,1,1,1,1]]).reshape(5,1)
p_sum=polytope(H,h)
vis([p_sum],0.5)

e=0.1
H=np.array([[1,0],[-1,0],[0,1],[0,-1]])
h=np.array([[e,e,0,1]]).reshape(4,1)
p_line=polytope(H,h)
vis([p_triangle,p_line],0.5)

n=2
H_S1=np.hstack((p_triangle.H,np.zeros((p_triangle.H.shape[0],n)),np.zeros((p_triangle.H.shape[0],n))))
H_S2=np.hstack((np.zeros((p_line.H.shape[0],n)),p_line.H,np.zeros((p_line.H.shape[0],n))))
H_I=np.hstack((np.eye(n),np.eye(n),-np.eye(n)))
H_S=np.vstack((H_S1,H_S2,H_I,-H_I))

h_S=np.vstack((p_triangle.h,p_line.h,np.zeros((n,1)),np.zeros((n,1))))

S=polytope(H_S,h_S)
G_r=np.hstack((np.zeros((n,n)),np.zeros((n,n)),np.eye(n)))
x_r=np.zeros((n,1))
x_l=np.zeros((n,1))
scale=0.65
G_l=np.eye(n)*scale
(alpha,beta,Lambda)=subset_both_projection(model,x_l,G_l,p_sum,x_r,G_r,S)
model.optimize()
alpha=valuation(alpha)
beta=valuation(beta)
Lambda=valuation(Lambda)

H=np.array([[1,1],[-2,1],[0,1],[0,-1]])
h=np.array([[0.3,0.5,0.5,0.7]]).reshape(4,1)
p_test=polytope(p_sum.H,p_sum.h*scale)
vis([p_sum,p_triangle,p_test,p_line],0.5)