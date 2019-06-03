#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May 30 15:49:45 2019

@author: sadra
"""
import numpy as np
# Pydrake
import pydrake.solvers.mathematicalprogram as MP
# Pypolycontain
from pypolycontain.lib.objects import H_polytope,zonotope,AH_polytope,hyperbox
from pypolycontain.lib.operations import Box,point_membership,directed_Hausdorff_distance,\
        check_non_empty,distance_polytopes,bounding_box,\
        distance_hyperbox,directed_Hausdorff_hyperbox
from pypolycontain.visualization.visualize_2D import visualize_2D_zonotopes as visZ
from pypolycontain.visualization.visualize_2D import visualize_2D_zonotopes_ax as visZ_ax


from pypolycontain.lib.hausdorff.hausdorff import Hausdorff_directed
from pypolycontain.lib.AH_polytope import minimum_distance

from time import time
#np.random.seed(0)
def test_memebership():
    # Test 1: # Random
    N,n,m=10,4,3
    H=np.random.random((N,n))-0.5
    h=np.random.random((N,1))+5
    T=np.random.random((m,n))
    t=np.random.random((m,1))*0
    H_P=H_polytope(H,h)
    P=AH_polytope(T,t,H_P)
    x=np.random.random((m,1))*0
    print point_membership(P,x,solver="gurobi")
    
    # Test 2: # Zonotope
    n=4
    P=zonotope(np.zeros((n,1)),np.eye(n),Box(n))
    x=np.random.random((n,1))
    print point_membership(P,x,solver="gurobi")
 
def test_emptyness():
    N,n=10,4
    H=np.random.random((N,n))-0.5
    h=np.random.random((N,1))+5  
    H_P=H_polytope(H,h)
    print check_non_empty(H_P)
    
def test_hausdorff():
    n=2
    q1=7
    q2=8
    z1=zonotope(np.random.random((n,1))*1,np.random.random((n,q1))-0.5,color='red')
    z2=zonotope(np.random.random((n,1))*1,np.random.random((n,q2))-0.5,color='blue')
    start=time()
    D=directed_Hausdorff_distance(z1,z2,solver="gurobi")
    print "Mathematical Program:",D,"\t",time()-start
    start=time()
    print "Gurobipy:",Hausdorff_directed(z1,z2),"\t",time()-start
    z3=zonotope(z2.x,np.hstack((z2.G,D*np.eye(n))),color='green')
    visZ([z3,z1,z2],a=0.5,alpha=0.2)
    return D

def test_distance():
    n=10
    q1=40
    q2=40
    z1=zonotope(np.random.random((n,1))*50,np.random.random((n,q1))-0.5,color='red')
    z2=zonotope(np.random.random((n,1))*70,np.random.random((n,q2))-0.5,color='blue')
    start=time()
    D=distance_polytopes(z1,z2,solver="gurobi")
    print "Mathematical Program:",D,"\t",time()-start
    start=time()
    print "Gurobipy:",minimum_distance(z1,z2),"\t",time()-start
#    import matplotlib.pyplot as plt
#    fig, ax = plt.subplots() # note we must use plt.subplots, not plt.subplot
#    visZ_ax(ax,[z1,z2],a=0.5,alpha=0.2)
    
def test_bounding_box():
    n,q=2,4
    z=zonotope(np.random.random((n,1))*1,np.random.random((n,q))-0.5,color='red')
    l,u=bounding_box(z)
    b=hyperbox(corners=(l,u))
    visZ([z,b.zonotope],a=0.5,alpha=0.2)
    
def test_box_distances():
    l1,l2=np.random.random((2,1)),np.random.random((2,1))*2
    u1,u2=l1+np.random.random((2,1)),l2+np.random.random((2,1))*3
    B1=hyperbox(corners=(l1,u1))
    B2=hyperbox(corners=(l2,u2))
    print "directed Hausroff B1,B2" ,directed_Hausdorff_hyperbox(B1,B2)
    print "directed Hausroff B2,B1 ",directed_Hausdorff_hyperbox(B2,B1)
    print "distance",distance_hyperbox(B1,B2),distance_hyperbox(B2,B1)
    B1.zonotope.color="red"
    B2.zonotope.color="blue"
    visZ([B1.zonotope,B2.zonotope],a=0.5,alpha=0.8)
    
    
def __main__():
#    test_memebership()
#    test_emptyness()
#    test_hausdorff()
#    test_distance()
#    test_bounding_box()
    test_box_distances()
    

if __name__=="__main__":
    __main__()