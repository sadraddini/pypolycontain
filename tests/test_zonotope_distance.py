# -*- coding: utf-8 -*-
"""
Created on Thu Nov  8 19:51:06 2018

@author: sadra
"""

import numpy as np

from pypolycontain.lib.polytope import translate
from pypolycontain.lib.elimination import project
from pypolycontain.visualization.visualize_2D import visualize_2D as vis

from pypolycontain.utils.utils import PI,valuation
from pypolycontain.lib.zonotope import zonotope,zonotope_distance


G_l=np.array([[1,1,0,-1,0.5,-1,1],[0,1,0.5,-0.5,-1,-5,4]])*1
#G_l=np.array([[1,0],[0,1]])*5
G_r=np.array([[1,0,1,-1,-2,-1,-3,1],[0,2,-1,3,3,-1,2,2]])*1
x_l=np.array([0,0]).reshape(2,1)
x_r=np.array([0,0]).reshape(2,1)

z_l=zonotope(x_l,G_l)
z_r=zonotope(x_r,G_r)

d_lr=zonotope_distance(z_l,z_r,eps=0.01)
d_rl=zonotope_distance(z_r,z_l,eps=0.01)

p_l=translate(project(G_l,PI(G_l.shape[1]),np.ones((2*G_l.shape[1],1))),x_l)
p_r=translate(project(G_r,PI(G_r.shape[1]),np.ones((2*G_r.shape[1],1))),x_r)
p_rd=translate(project(np.hstack((G_r,d_lr*np.eye(G_r.shape[0]))),PI(G_r.shape[1]+G_r.shape[0]),np.ones((2*G_r.shape[1]+2*G_r.shape[0],1))),x_r)
p_ld=translate(project(np.hstack((G_l,d_rl*np.eye(G_l.shape[0]))),PI(G_l.shape[1]+G_l.shape[0]),np.ones((2*G_l.shape[1]+2*G_l.shape[0],1))),x_l)

vis([p_r,p_l])
vis([p_rd,p_r,p_l])
vis([p_ld,p_r,p_l])