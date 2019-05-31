#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May 30 10:45:14 2019

@author: sadra
"""

import numpy as np
# Pydrake
import pydrake.solvers.mathematicalprogram as MP
import pydrake.solvers.gurobi as Gurobi_drake
# Pypolycontain
from pypolycontain.lib.objects import AH_polytope,Box
# use Gurobi solver
global gurobi_solver
gurobi_solver=Gurobi_drake.GurobiSolver()

def to_AH_polytope(P):
    if P.type=="AH_polytope":
        return P
    elif P.type=="H-polytope":
        n=P.H.shape[1]
        return AH_polytope(np.eye(n),np.zeros((n,1)),P)
    elif P.type=="zonotope":
        q=P.G.shape[1]
        return AH_polytope(P.G,P.x,Box(q))
    else:
        raise ValueError("P type not understood:",P.type)
        
def point_membership(Q,x,tol=10**-5,solver=gurobi_solver):
    if Q.type=="H_polytope":
        return Q.if_inside(x,tol)
    else:
        Q=to_AH_polytope(Q)
        prog=MP.MathematicalProgram()
        zeta=prog.NewContinuousVariables(Q.P.H.shape[1],1,"zeta")
        prog.AddLinearConstraint(A=Q.P.H,ub=Q.P.h+tol,lb=-np.inf*np.ones((Q.P.h.shape[1],1)),vars=zeta)
        prog.AddLinearEqualityConstraint(Q.T,x-Q.t,zeta)
        if solver==gurobi_solver:
            result=solver.Solve(prog,None,None)
        else:
            result=MP.Solve(prog)
    return result.is_success()

