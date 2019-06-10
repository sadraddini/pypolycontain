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
from pypolycontain.lib.objects import AH_polytope,Box,hyperbox
# use Gurobi solver
global gurobi_solver,OSQP_solver
gurobi_solver=Gurobi_drake.GurobiSolver()
OSQP_solver=OSQP_drake.OsqpSolver()

def to_AH_polytope(P):
    if P.type=="AH_polytope":
        return P
    elif P.type=="H_polytope" or P.type=="H-polytope":
        n=P.H.shape[1]
        return AH_polytope(np.eye(n),np.zeros((n,1)),P)
    elif P.type=="zonotope":
        q=P.G.shape[1]
        return AH_polytope(P.G,P.x,Box(N=q))
    else:
        raise ValueError("P type not understood:",P.type)
        
def point_membership(Q,x,tol=10**-5,solver="gurobi"):
    if Q.type=="H_polytope":
        return Q.if_inside(x,tol)
    else:
        Q=to_AH_polytope(Q)
        prog=MP.MathematicalProgram()
        zeta=prog.NewContinuousVariables(Q.P.H.shape[1],1,"zeta")
        prog.AddLinearConstraint(A=Q.P.H,ub=Q.P.h+tol,lb=-np.inf*np.ones((Q.P.h.shape[0],1)),vars=zeta)
        prog.AddLinearEqualityConstraint(Q.T,x-Q.t,zeta)
        if solver=="gurobi":
            result=gurobi_solver.Solve(prog,None,None)
        elif solver=="osqp":
            prog.AddQuadraticCost(np.eye(zeta.shape[0]),np.zeros(zeta.shape),zeta)
            result=OSQP_solver.Solve(prog,None,None)
        else:
            result=MP.Solve(prog)
    return result.is_success()

def check_non_empty(Q,tol=10**-5,solver="gurobi"):
    Q=to_AH_polytope(Q)
    prog=MP.MathematicalProgram()
    zeta=prog.NewContinuousVariables(Q.P.H.shape[1],1,"zeta")
    prog.AddLinearConstraint(A=Q.P.H,ub=Q.P.h+tol,lb=-np.inf*np.ones((Q.P.h.shape[0],1)),vars=zeta)
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
    Computes the directed Hausdorff distance of Q_1 and Q_2 (AH_polytopes)
    
                        Minimize    epsilon  
                        such that   Q1 \subset Q2+epsilon(Ball)
                        
    It is zero if and only if Q1 subset Q2. The method is based on 
                
                    Sadraddini&Tedrake, 2019, CDC (available on ArXiv)
                    
    We solve the following problem:
        D*ball+Q1 subset Q2
    We solve the following linear program:
        min     D
        s.t.    Lambda_1 H_1=H_2 Gamma_1
                Lambda_2 H_1=H_ball Gamma_2
                Lambda_1 h_1<=h_2 + H_2 beta_1
                Lambda_2 h_2<=D h_ball + H_ball beta_2
                x_2 - X_2 beta_1 - beta_2 = x_1
                X_2 Gamma_1 + Gamma_2 = X_1
    ***************************************************************************
    """
    Q1,Q2=to_AH_polytope(Q1),to_AH_polytope(Q2)
    n=Q1.t.shape[0]
    if ball=="infinty_norm":
        HB=np.vstack((np.eye(n),-np.eye(n)))
        hB=np.vstack((np.ones((n,1)),np.ones((n,1))))
    elif ball=="l1":
        HB,hb=make_ball(ball)
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
    prog.AddBoundingBoxConstraint(0,np.inf,Lambda_1)
    prog.AddBoundingBoxConstraint(0,np.inf,Lambda_2)
    # Lambda_1 H_1
    Lambda_H_Gamma(prog,Lambda_1,Q1.P.H,Q2.P.H,Gamma_1)
    # Lambda_2 H_1
    Lambda_H_Gamma(prog,Lambda_2,Q1.P.H,HB,Gamma_2)
    # Lambda_1 h_1
    Lambda_h_Inequality(prog,Lambda_1,beta_1,Q2.P.H,Q1.P.h,Q2.P.h)
    # Lambda_2 h_1
    Lambda_h_Inequality_D(prog,Lambda_2,beta_2,HB,Q1.P.h,hB,D)
    # X2 beta_1   
    prog.AddLinearEqualityConstraint(-np.hstack((Q2.T,np.eye(n))),Q1.t-Q2.t,np.vstack((beta_1,beta_2)))
    # X2 Gamma_1
    Aeq=np.hstack((Q2.T,np.eye(Q2.T.shape[0])))
    for i in range(Gamma_1.shape[1]):
        beq=Q1.T[:,i]
        var=np.hstack((Gamma_1[:,i],Gamma_2[:,i]))
        prog.AddLinearEqualityConstraint(Aeq,beq,var)
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

def Hausdorff_distance(Q1,Q2,ball="infinty_norm",solver="gurobi"):
    return max(directed_Hausdorff_distance(Q1,Q2,ball,solver),directed_Hausdorff_distance(Q2,Q1,ball,solver))
    
def distance_polytopes(Q1,Q2,ball="infinity",solver="Gurobi"):
    Q1,Q2=to_AH_polytope(Q1),to_AH_polytope(Q2)
    n=Q1.n
    prog=MP.MathematicalProgram()
    zeta1=prog.NewContinuousVariables(Q1.P.H.shape[1],1,"zeta1")
    zeta2=prog.NewContinuousVariables(Q2.P.H.shape[1],1,"zeta2")
    delta=prog.NewContinuousVariables(n,1,"delta")
    prog.AddLinearConstraint(A=Q1.P.H,ub=Q1.P.h,lb=-np.inf*np.ones((Q1.P.h.shape[0],1)),vars=zeta1)
    prog.AddLinearConstraint(A=Q2.P.H,ub=Q2.P.h,lb=-np.inf*np.ones((Q2.P.h.shape[0],1)),vars=zeta2)
    prog.AddLinearEqualityConstraint( np.hstack((Q1.T,-Q2.T,np.eye(n))),Q2.t-Q1.t,np.vstack((zeta1,zeta2,delta)) )
    if ball=="infinity":
        delta_abs=prog.NewContinuousVariables(1,1,"delta_abs")
        prog.AddBoundingBoxConstraint(0,np.inf,delta_abs)
        prog.AddLinearConstraint(np.greater_equal( np.dot(np.ones((n,1)),delta_abs),delta,dtype='object' ))
        prog.AddLinearConstraint(np.greater_equal( np.dot(np.ones((n,1)),delta_abs),-delta,dtype='object' ))
        cost=delta_abs
    elif ball=="l1":
        delta_abs=prog.NewContinuousVariables(n,1,"delta_abs")
        prog.AddBoundingBoxConstraint(0,np.inf,delta_abs)
        prog.AddLinearConstraint(np.greater_equal( delta_abs,delta,dtype='object' ))
        prog.AddLinearConstraint(np.greater_equal( delta_abs,-delta,dtype='object' ))
        cost=np.dot(np.ones((1,n)),delta_abs)
    else:
        raise NotImplementedError
    if solver=="gurobi":
        prog.AddLinearCost(cost[0,0])
        result=gurobi_solver.Solve(prog,None,None)
    elif solver=="osqp":
        prog.AddQuadraticCost(cost[0,0]*cost[0,0])
        result=OSQP_solver.Solve(prog,None,None)
    else:
        prog.AddLinearCost(cost[0,0])
        result=MP.Solve(prog)
    if result.is_success():
        return np.sum(result.GetSolution(delta_abs)),\
            np.dot(Q1.T,result.GetSolution(zeta1).reshape(zeta1.shape[0],1))+Q1.t,\
            np.dot(Q2.T,result.GetSolution(zeta2).reshape(zeta2.shape[0],1))+Q2.t
    
def bounding_box(Q,solver="Gurobi"):
    Q=to_AH_polytope(Q)
    prog=MP.MathematicalProgram()
    zeta=prog.NewContinuousVariables(Q.P.H.shape[1],1,"zeta")
    x=prog.NewContinuousVariables(Q.n,1,"x")
    prog.AddLinearConstraint(A=Q.P.H,ub=Q.P.h,lb=-np.inf*np.ones((Q.P.h.shape[0],1)),vars=zeta)
    prog.AddLinearEqualityConstraint(np.hstack((Q.T,np.eye(Q.n))),Q.t,np.vstack((zeta,x)))
    lower_corner=np.zeros((Q.n,1))
    upper_corner=np.zeros((Q.n,1))
    c=prog.AddLinearCost(np.dot(np.ones((1,Q.n)),x)[0,0])
    if solver=="Gurobi":
        solver=gurobi_solver
    else:
        raise NotImplementedError
    a=np.zeros((Q.n,1))
    # Lower Corners
    for i in range(Q.n):
        e=c.evaluator()
        a[i,0]=1
        e.UpdateCoefficients(a.reshape(Q.n))
        result=solver.Solve(prog,None,None)
        assert result.is_success()
        lower_corner[i,0]=result.GetSolution(x)[i]
        a[i,0]=0
    # Upper Corners
    for i in range(Q.n):
        e=c.evaluator()
        a[i,0]=-1
        e.UpdateCoefficients(a)
        result=solver.Solve(prog,None,None)
        assert result.is_success()
        upper_corner[i,0]=result.GetSolution(x)[i]
        a[i,0]=0
    return hyperbox(corners=(lower_corner,upper_corner))
        
        
def directed_Hausdorff_hyperbox(b1,b2):
    """
    The directed Hausdorff hyperbox 
    min epsilon such that b1 \in b2+epsilon
    """       
    return max(0,np.max(np.hstack((b1.u-b2.u,b2.l-b1.l))))           
    
def distance_hyperbox(b1,b2):
    """
    The distance between boxes
    """
    return max(0,np.max(np.hstack((b1.l-b2.u,b2.l-b1.u))))      
    

def make_ball(n,norm):
    if norm=="l1":
        pass
    elif norm=="infinity":
        pass
    return 


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
#    q=Lambda.shape[0]
    mathematical_program.AddBoundingBoxConstraint(0,np.inf,Lambda)
#    [mathematical_program.AddLinearConstraint(A=np.eye(q),vars=Lambda[:,i],
#                                              ub=np.inf*np.ones((q,1)),lb=np.zeros((q,1)))
#                                                for i in range(Lambda.shape[1])]
#    q=Lambda.shape[0]*Lambda.shape[1]
#    mathematical_program.AddLinearConstraint(A=np.eye(q),vars=Lambda.reshape(q),ub=np.inf*np.ones((q)),lb=np.zeros((q)))
        
def Lambda_H_Gamma(mathematical_program,Lambda,H_1,H_2,Gamma):
    """
    Lambda H_1 = H_2 Gamma
    """
#    v_1=Lambda.reshape((Lambda.shape[0]*Lambda.shape[1],1))
#    for j in range(Gamma.shape[1]):
#        M1=np.kron(np.eye(H_2.shape[0]),H_1[:,j].reshape(1,H_1.shape[0]))
#        M2=-H_2
#        M=np.hstack((M1,M2))
#        v=np.vstack((v_1,Gamma[:,j].reshape(Gamma.shape[0],1)))
#        mathematical_program.AddLinearEqualityConstraint(M,np.zeros((M.shape[0],1)),v)
    for i in range(Lambda.shape[0]):
        for j in range(Gamma.shape[1]):
            M=np.hstack((H_1[:,j],-H_2[i,:]))
            v=np.hstack((Lambda[i,:],Gamma[:,j]))
            mathematical_program.AddLinearEqualityConstraint(M.reshape(1,M.shape[0]),np.zeros(1),v)