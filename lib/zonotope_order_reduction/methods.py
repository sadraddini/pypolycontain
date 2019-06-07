#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  7 16:29:04 2019

@author: sadra
"""

# Numpy ans scipy
import numpy as np
import scipy.linalg as spa
# pydrake
import pydrake.solvers.mathematicalprogram as MP
import pydrake.solvers.gurobi as Gurobi_drake
import pydrake.solvers.osqp as OSQP_drake
# Pypolycontain
from pypolycontain.lib.objects import zonotope
# use Gurobi solver
global gurobi_solver,OSQP_solver
gurobi_solver=Gurobi_drake.GurobiSolver()
OSQP_solver=OSQP_drake.OsqpSolver()

def G_cut(myzonotope,number_of_columns_wanted,solver="gurobi"):
    q_i,q_f=myzonotope.G.shape[1],number_of_columns_wanted
    G_norm=np.linalg.norm(myzonotope.G,np.inf, axis=0)#-np.linalg.norm(myzonotope.G,np.inf, axis=0)
    G_sorted=myzonotope.G[:,np.argsort(G_norm)]
    G_cut=G_sorted[:,q_i-q_f:]
    # Now encode G=G_cut * Gamma
    prog=MP.MathematicalProgram()
    D=prog.NewContinuousVariables(q_f,1,"D")
    Gamma=prog.NewContinuousVariables(q_f,q_i,"Gamma")
    Gamma_abs=prog.NewContinuousVariables(q_f,q_i,"Gamma_abs")
    prog.AddLinearConstraint(np.greater_equal(Gamma_abs,Gamma,dtype='object').flatten())
    prog.AddLinearConstraint(np.greater_equal(Gamma_abs,-Gamma,dtype='object').flatten())
    prog.AddLinearConstraint(np.less_equal(np.dot(Gamma_abs,np.ones((q_i,1))),D,dtype='object'))
    prog.AddLinearConstraint(np.equal(myzonotope.G,np.dot(G_cut,Gamma),dtype='object').flatten())
    if solver=="gurobi":
        prog.AddLinearCost(np.sum(D))
        result=gurobi_solver.Solve(prog,None,None)
    elif solver=="osqp":
        prog.AddQuadraticCost(np.eye(q_f),-1*np.ones(D.shape),D)
        result=OSQP_solver.Solve(prog,None,None)
    else:
        result=MP.Solve(prog)
    if result.is_success():
        D_num=result.GetSolution(D)
#        print D_num
        G_num=np.dot(G_cut,np.diag(D_num))
        return D_num,G_num