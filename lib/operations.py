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
import pydrake.solvers.osqp as OSQP_drake

# Pypolycontain
from pypolycontain.lib.objects import AH_polytope,Box
# use Gurobi solver
global gurobi_solver,OSQP_solver
gurobi_solver=Gurobi_drake.GurobiSolver()
OSQP_solver=OSQP_drake.OsqpSolver()

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
        
def point_membership(Q,x,tol=10**-5,solver="gurobi"):
    if Q.type=="H_polytope":
        return Q.if_inside(x,tol)
    else:
        Q=to_AH_polytope(Q)
        prog=MP.MathematicalProgram()
        zeta=prog.NewContinuousVariables(Q.P.H.shape[1],1,"zeta")
        prog.AddLinearConstraint(A=Q.P.H,ub=Q.P.h+tol,lb=-np.inf*np.ones((Q.P.h.shape[1],1)),vars=zeta)
        prog.AddLinearEqualityConstraint(Q.T,x-Q.t,zeta)
        if solver=="gurobi":
            result=gurobi_solver.Solve(prog,None,None)
        elif solver=="osqp":
            prog.AddQuadraticCost(np.eye(zeta.shape[0]),np.zeros(zeta.shape),zeta)
            result=OSQP_solver.Solve(prog,None,None)
        else:
            result=MP.Solve(prog)
    return result.is_success()

def directed_Hausdorff_distance(Q1,Q2,ball="infinty_norm",solver="gurobi"):
    """
    ***************************************************************************
    Computes the directed Hausdorff distance of Q_1 and Q_2
    Minimum epsilon such that 
                                    Q1 \subset Q2+epsilon(Ball)
    zero if and only if Q1 subset Q2. The method is based on Sadraddini&Tedrake
    , 2019, CDC (available on ArXiv)
    We solve the following problem:
        D*ball+Q1 subset Q2
    We solve the following linear program:
        min     D
        s.t.    Lambda_1 H_1=H_2 Gamma_1
                Lambda_2 H_1=H_ball Gamma_2
                Lambda_1 h_1<=h_2 + H_2 beta_1
                Lambda_2 h_2<=D h_ball + H_ball beta_2
    ***************************************************************************
    """
    Q1,Q2=to_AH_polytope(Q1),to_AH_polytope(Q2)
    n=Q1.t.shape[0]
    if ball=="infinty_norm":
        HB=np.vstack((np.eye(n),-np.eye(n)))
        hB=np.vstack((np.ones((n,1)),np.ones((n,1))))
    prog=MP.MathematicalProgram()
    # Variables
    D=prog.NewContinuousVariables(1,1,"D")
    Lambda_1=prog.NewContinuousVariables(Q2.P.H.shape[0],Q1.P.H.shape[0],"Lambda_1")
    Lambda_2=prog.NewContinuousVariables(HB.shape[0],Q1.P.H.shape[0],"Lambda2")
    Gamma_1=prog.NewContinuousVariables(Q2.P.H.shape[1],Q1.P.H.shape[1],"Gamma1")
    Gamma_2=prog.NewContinuousVariables(HB.shape[1],Q1.P.H.shape[1],"Gamma1")
    beta_1=prog.NewContinuousVariables(Q2.P.H.shape[1],1,"beta1")
    beta_2=prog.NewContinuousVariables(HB.shape[1],1,"beta1")
    # Constraints
    # Lambda_1 and Lambda_2 positive
    [prog.AddLinearConstraint(v=Lambda_1[:,i],ub=np.inf*np.ones((Lambda_1.shape[0],1)),
                              lb=np.zeros((Lambda_1.shape[0],1))) for i in range(Lambda_1.shape[1])]
    [prog.AddLinearConstraint(v=Lambda_2[:,i],ub=np.inf*np.ones((Lambda_2.shape[0],1)),
                              lb=np.zeros((Lambda_2.shape[0],1))) for i in range(Lambda_2.shape[1])] 
    # Lambda_1 H_1
#    f=np.equal(np.dot(Lambda_1,Q1.P.H),np.dot(Q2.P.H,Gamma_1),dtype=object)
#    [prog.AddLinearConstraint(f[:,i]) for i in range(f.shape[1])]
    [prog.AddLinearEqualityConstraint(np.hstack((Q1.P.H[:,j].T,-Q2.P.H[i,:])).reshape(1,Lambda_1.shape[1]+Gamma_1.shape[0]),
                                     np.zeros((1)),
                                     np.hstack((Lambda_1[i,:],Gamma_1[:,j])))\
            for i in range(Lambda_1.shape[0]) for j in range(Gamma_1.shape[1])]
    # Lambda_2 H_1
#    f=np.equal(np.dot(Lambda_2,Q1.P.H),np.dot(HB,Gamma_2),dtype=object)
#    [prog.AddLinearConstraint(f[:,i]) for i in range(f.shape[1])]
    [prog.AddLinearEqualityConstraint(np.hstack((Q1.P.H[:,j].T,-HB[i,:])).reshape(1,Lambda_2.shape[1]+Gamma_2.shape[0]),
                                     np.zeros((1)),
                                     np.hstack((Lambda_2[i,:],Gamma_2[:,j])))\
            for i in range(Lambda_2.shape[0]) for j in range(Gamma_2.shape[1])]
    # Lambda_1 h_1
    f=np.less_equal(np.dot(Lambda_1,Q1.P.h),Q2.P.h+np.dot(Q2.P.H,beta_1),dtype=object)
    prog.AddLinearConstraint(f)
    # Lambda_2 h_1
    f=np.less_equal(np.dot(Lambda_2,Q1.P.h),np.dot(hB,D)+np.dot(HB,beta_2),dtype=object)
    prog.AddLinearConstraint(f)
    # X2 beta_1   
#    f=np.equal(Q2.t-np.dot(Q2.T,beta_1)-beta_2,Q1.t,dtype=object)
#    prog.AddLinearConstraint(f)
    prog.AddLinearEqualityConstraint(-np.hstack((Q2.T,np.eye(n))),Q1.t-Q2.t,np.vstack((beta_1,beta_2)))
    # X2 Gamma_1
    f=np.equal(np.dot(Q2.T,Gamma_1)+Gamma_2,Q1.T,dtype=object)
    [prog.AddLinearConstraint(f[:,i]) for i in range(f.shape[1])]
    # Cost
    # Optimize
    if solver=="gurobi":
            prog.AddLinearCost(D[0,0])
            result=gurobi_solver.Solve(prog,None,None)
    elif solver=="osqp":
        prog.AddQuadraticCost(D[0,0]*D[0,0])
        result=OSQP_solver.Solve(prog,None,None)
    else:
        result=MP.Solve(prog)
    if result.is_success():
        return np.asscalar(result.GetSolution(D))
    
    
"""
Pydrake Mathematical Program Helper: Matrix based Constraints
"""
def AddMatrixInequalityConstraint_classical(mathematical_program,A,X,B):
    raise NotImplementedError    