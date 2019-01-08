#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  7 15:55:37 2019

@author: sadra
"""

import numpy as np

from pypolycontain.lib.orthogonal_projection import orthogonal_projection_fixed_Hx,orthogonal_projection_gradient_decent,gradient_decent
from pypolycontain.visualization.visualize_2D import visualize_2D as vis2
from pypolycontain.lib.polytope import polytope

from pyfomo.src.main import fourier_motzkin_eliminate_single


A=np.array([[1,0.2],[0,1]])
B=np.array([0,1]).reshape(2,1)
N=10
C_x=np.array([[1,0],[0,1],[-1,0],[0,-1]])
c_x=np.array([[1,1,1,1]]).T
C_goal=np.array([[1,0],[0,1],[-1,0],[0,-1]])
c_goal=np.array([[1,1,1,1]]).T*0.1
C_u=np.array([1,-1]).reshape(2,1)
c_u=np.array([1,1]).reshape(2,1)

#C_goal=np.eye(2)
#C_x=C_goal

H=np.dot(C_goal,np.linalg.matrix_power(A, N))
for i in range(N+1):
    H=np.vstack( (H,np.dot(C_x,np.linalg.matrix_power(A, N-i)) ))
    
F=np.hstack( [np.dot(C_goal,np.dot(np.linalg.matrix_power(A, N-i-1),B)) for i in range(N)])
for i in range(N):
    F_temp=np.zeros((C_x.shape[0],N*B.shape[1]))
    F_temp[:,0:(N-i)*B.shape[1]]=np.hstack( [np.dot(C_x,np.dot(np.linalg.matrix_power(A, N-i-j-1),B)) for j in range(N-i)])
    F=np.vstack((F,F_temp))

F_temp=np.zeros((C_x.shape[0],N*B.shape[1]))
F=np.vstack((F,F_temp))

H=np.vstack ( ( H,np.zeros((N*C_u.shape[0],B.shape[0])) ) )
F=np.vstack ( ( F, np.kron(np.eye(N),C_u)) )

g=np.vstack((c_goal,np.vstack([c_x for i in range(N+1)]),np.vstack([c_u for i in range(N)])))

if False:
    A=np.hstack((H,F))
    b=g
    
    for i in range(10):
        var_index=A.shape[1]-1
        print "i"
        (A,b)=fourier_motzkin_eliminate_single(var_index,A,b,atol=10**-8)
    #    print("A_new=",A_new)
    #    print("b_new=",b_new)
    p_real=polytope(A,b)

#H=np.array([[1,0],[0,1],[-1,0],[0,-1],[0,0],[0,0]])
#F=np.array([[1],[0],[-1],[0],[1],[-1]])
#g=np.array([[1,1,1,1,1,1]]).T
#
Hx_correct=np.array([[1,0],[0,1],[-1,0],[0,-1]])
#Hx=np.array([[0.5,0],[0,1],[-0.5,0],[0,-1]])*3
Hx=np.array([[1,0],[-1,0],[0,-1],[0,1]])*2
xbar=np.zeros((2,1))
#xbar=-np.array([1.8,0.8]).reshape(2,1)

delta=0.05
output=gradient_decent(Hx,H,F,g,xbar,delta,N=50)
import matplotlib.pyplot as plt
plt.plot([x[1] for x in output])
Hx=output[-1][0]

p_out=polytope(Hx_correct,np.ones((Hx_correct.shape[0],1)))
for i in range(len(output)):
    Hx=output[i][0]
    p_in=polytope(Hx,np.ones((Hx.shape[0],1)))
    fig=vis2([p_out,p_real,p_in],title=r"Orthogonal Projection Iteration %d - $d_H=%f$"%(i,output[i][1]))
    fig.savefig('figures/orthogonal_projection_MPC %d'%i,dpi=100)