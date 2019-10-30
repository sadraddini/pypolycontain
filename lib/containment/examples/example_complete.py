#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 30 13:35:37 2019

@author: sadra
"""
import numpy as np
import warnings


# Pydrake
try:
    import pydrake.solvers.mathematicalprogram as MP
    import pydrake.solvers.gurobi as Gurobi_drake
    import pydrake.solvers.osqp as OSQP_drake
    # use Gurobi solver
    global gurobi_solver,OSQP_solver, license
    gurobi_solver=Gurobi_drake.GurobiSolver()
    license = gurobi_solver.AcquireLicense()
    OSQP_solver=OSQP_drake.OsqpSolver()
except:
    warnings.warn("You don't have pydrake installed properly. Methods that rely on optimization may fail.")
    
#pypolycontain    
try:
    from pypolycontain.lib.containment.complete import extreme_rays_for_containment as erc
    from pypolycontain.lib.containment.complete import subset
    from pypolycontain.lib.objects import H_polytope,AH_polytope
    from pypolycontain.visualization.visualize_2D import visualize_2D as vis
    from pypolycontain.visualization.visualize_2D import visualize_2D_AH_polytope as visAH
    from pypolycontain.lib.operations import minkowski_sum
except:
    warnings.warn("You don't have pyplycontain not properly installed.")
    

def test_extreme_rays_generation():
    n=5
    q=3  
    H_y=np.vstack((np.eye(n),-np.eye(n)))
    Y=np.random.random((q,n))-0.5
    #Y=np.eye(2)
    circombody=AH_polytope(t=np.zeros((Y.shape[0],1)),T=Y,P=H_polytope(H=H_y,h=np.zeros((Y.shape[0],1))))
    Theta=erc(circombody,N=0)
    print(Theta.shape)
    print(np.all(Theta>-10**(-8)))

def test_triangle():
    H=np.array([[1,1],[-1,1],[0,-1]])
    h=np.array([[1,1,0]]).reshape(3,1)
    p_triangle=H_polytope(H,h)
    
    e=0.1
    scale=0.71
    H=np.array([[1,0],[-1,0],[0,-1],[1,1],[-1,1],[0,1]])
    h=np.array([[1+e,1+e,1,1+e,1+e,1]]).reshape(6,1)
    p_sum=H_polytope(H,h,color='green')
    #vis([p_sum],0.5)
    
    H=np.array([[1,0],[-1,0],[0,1],[0,-1]])
    h=np.array([[e,e,0,1]]).reshape(4,1)
    p_line=H_polytope(H,h,color='blue')
    #vis([p_triangle,p_line],0.5)
    
    #H=np.array([[1,1],[-2,1],[0,1],[0,-1]])
    #h=np.array([[0.3,0.5,0.5,0.7]]).reshape(4,1)
    p_test=H_polytope(p_sum.H,p_sum.h*scale)
    
    Q=minkowski_sum(p_line,p_triangle)
    Q.color='cyan'
    #visAH([Q],N=100,a=0.2)
    
    N=2
    Theta=erc(Q,N)
    
    prog=MP.MathematicalProgram()
    subset(prog,p_test,Q,N=N)
    result=gurobi_solver.Solve(prog,None,None)
    if result.is_success():
        print("subset test successfull")
    else:
        print("not a subset")