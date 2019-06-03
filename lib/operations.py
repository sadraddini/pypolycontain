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

import time

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
    start=time.time()
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
    print "Variables and initial",time.time()-start
    # Constraints
    # Lambda_1 and Lambda_2 positive
    start=time.time()
    positive_matrix(prog,Lambda_1)
    positive_matrix(prog,Lambda_2)
    print "Lambda Positive",time.time()-start
    # Lambda_1 H_1
    start=time.time()
#    f=np.equal(np.dot(Lambda_1,Q1.P.H),np.dot(Q2.P.H,Gamma_1),dtype=object)
#    [prog.AddLinearConstraint(f[:,i]) for i in range(f.shape[1])]
    Lambda_H_Gamma(prog,Lambda_1,Q1.P.H,Q2.P.H,Gamma_1)
#    [prog.AddLinearEqualityConstraint(np.hstack((Q1.P.H[:,j].T,-Q2.P.H[i,:])).reshape(1,Lambda_1.shape[1]+Gamma_1.shape[0]),
#                                     np.zeros((1)),
#                                     np.hstack((Lambda_1[i,:],Gamma_1[:,j])))\
#            for i in range(Lambda_1.shape[0]) for j in range(Gamma_1.shape[1])]
    # Lambda_2 H_1
#    f=np.equal(np.dot(Lambda_2,Q1.P.H),np.dot(HB,Gamma_2),dtype=object)
#    [prog.AddLinearConstraint(f[:,i]) for i in range(f.shape[1])]
    Lambda_H_Gamma(prog,Lambda_2,Q1.P.H,HB,Gamma_2)
#    [prog.AddLinearEqualityConstraint(np.hstack((Q1.P.H[:,j].T,-HB[i,:])).reshape(1,Lambda_2.shape[1]+Gamma_2.shape[0]),
#                                     np.zeros((1)),
#                                     np.hstack((Lambda_2[i,:],Gamma_2[:,j])))\
#            for i in range(Lambda_2.shape[0]) for j in range(Gamma_2.shape[1])]
    print "Lambda Gamma (auxilary)",time.time()-start
    start=time.time()
    # Lambda_1 h_1
#    f=np.less_equal(np.dot(Lambda_1,Q1.P.h),Q2.P.h+np.dot(Q2.P.H,beta_1),dtype=object)
#    prog.AddLinearConstraint(f)
    Lambda_h_Inequality(prog,Lambda_1,beta_1,Q2.P.H,Q1.P.h,Q2.P.h)
    # Lambda_2 h_1
#    f=np.less_equal(np.dot(Lambda_2,Q1.P.h),np.dot(hB,D)+np.dot(HB,beta_2),dtype=object)
#    prog.AddLinearConstraint(f)
    Lambda_h_Inequality_D(prog,Lambda_2,beta_2,HB,Q1.P.h,hB,D)
    print "Lambda h (optimized)",time.time()-start
    # X2 beta_1   
#    f=np.equal(Q2.t-np.dot(Q2.T,beta_1)-beta_2,Q1.t,dtype=object)
#    prog.AddLinearConstraint(f)
    start=time.time()
    prog.AddLinearEqualityConstraint(-np.hstack((Q2.T,np.eye(n))),Q1.t-Q2.t,np.vstack((beta_1,beta_2)))
    print "X beta (optimized)",time.time()-start
    # X2 Gamma_1
    start=time.time()
#    f=np.equal(np.dot(Q2.T,Gamma_1)+Gamma_2,Q1.T,dtype=object)
#    [prog.AddLinearConstraint(f[:,i]) for i in range(f.shape[1])]
#    print Q2.T.shape,Gamma_1.shape,Gamma_2.shape,Q1.T.shape
    Aeq=np.hstack((Q2.T,np.eye(Q2.T.shape[0])))
    for i in range(Gamma_1.shape[1]):
        beq=Q1.T[:,i]
        var=np.hstack((Gamma_1[:,i],Gamma_2[:,i]))
#        print Aeq.shape,beq.shape,var.shape
        prog.AddLinearEqualityConstraint(Aeq,beq,var)
    print "X Gamma (optimized)",time.time()-start
#    print "X Gamma",time.time()-start
    # Cost
    # Optimize
    start=time.time()
    if solver=="gurobi":
            prog.AddLinearCost(D[0,0])
            result=gurobi_solver.Solve(prog,None,None)
    elif solver=="osqp":
        prog.AddQuadraticCost(D[0,0]*D[0,0])
        result=OSQP_solver.Solve(prog,None,None)
    else:
        result=MP.Solve(prog)
    if result.is_success():
        print "optimize",time.time()-start
        return np.asscalar(result.GetSolution(D))
    
    
"""
Pydrake Mathematical Program Helper: Matrix based Constraints
"""
def AddMatrixInequalityConstraint_classical(mathematical_program,A,X,B):
    raise NotImplementedError    
    
def Lambda_h_Inequality(mathematical_program,Lambda,beta,H,h_1,h_2):
    """
    Adds Lambda H-1 \le h_2 + H beta to the Mathematical Program
    """
    M1=np.kron(np.eye(h_2.shape[0]),h_1.T)
    M=np.hstack((M1,-H))
    v1=Lambda.reshape((Lambda.shape[0]*Lambda.shape[1],1))
    v=np.vstack((v1,beta))
    mathematical_program.AddLinearConstraint(A=M,ub=h_2,lb=-np.inf*np.ones(h_2.shape),vars=v)
    
def Lambda_h_Inequality_D(mathematical_program,Lambda,beta,H,h_1,h_2,D):
    """
    Adds Lambda H-1 \le h_2 D + H beta to the Mathematical Program
    """
    M1=np.kron(np.eye(h_2.shape[0]),h_1.T)
    M=np.hstack((M1,-H,-h_2))
    v1=Lambda.reshape((Lambda.shape[0]*Lambda.shape[1],1))
    v=np.vstack((v1,beta,D))
    mathematical_program.AddLinearConstraint(A=M,ub=np.zeros(h_2.shape),lb=-np.inf*np.ones(h_2.shape),vars=v)
    
def positive_matrix(mathematical_program,Lambda):
    """
    All elements are non-negative 
    """
    q=Lambda.shape[0]
    [mathematical_program.AddLinearConstraint(A=np.eye(q),vars=Lambda[:,i],
                                              ub=np.inf*np.ones((q,1)),lb=np.zeros((q,1)))
                                                for i in range(Lambda.shape[1])]
#    q=Lambda.shape[0]*Lambda.shape[1]
#    mathematical_program.AddLinearConstraint(A=np.eye(q),vars=Lambda.reshape(q),ub=np.inf*np.ones((q)),lb=np.zeros((q)))
        
def Lambda_H_Gamma(mathematical_program,Lambda,H_1,H_2,Gamma):
    """
    Lambda H_1 = H_2 Gamma
    """
    v_1=Lambda.reshape((Lambda.shape[0]*Lambda.shape[1],1))
    for j in range(Gamma.shape[1]):
        M1=np.kron(np.eye(H_2.shape[0]),H_1[:,j].reshape(1,H_1.shape[0]))
        M2=-H_2
        M=np.hstack((M1,M2))
        v=np.vstack((v_1,Gamma[:,j].reshape(Gamma.shape[0],1)))
        mathematical_program.AddLinearEqualityConstraint(M,np.zeros((M.shape[0],1)),v)