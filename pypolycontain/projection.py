#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  3 19:57:12 2020

@author: sadra
"""

import warnings
import numpy as np
from itertools import combinations

# Scipy
try:
    import scipy.linalg as spa
    from scipy.spatial import ConvexHull
    from scipy.linalg import block_diag
except:
    warnings.warn("You don't have scipy package installed. You may get error while using some feautures.")

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
    import pydrake.solvers.scs as SCS
    scs_solver=SCS.ScsSolver()
except:
    warnings.warn("You don't have pydrake installed properly. Methods that rely on optimization may fail.")
    

# Pypolycontain
try:
    import pypolycontain as pp
except:
    warnings.warn("You don't have pypolycontain properly installed. Can not import objects")



def ray_shooting_hyperplanes_older(Q,N=0,H_rays=None):
    """
    Ray Shooting to find an outer-approximation of the AH-polytope
    """
    prog=MP.MathematicalProgram()
    Q=pp.to_AH_polytope(Q)

    if type(H_rays)==type(None):
        if N==0:
            N=2**(Q.n-1) # This many hyperplanes I want :) 
        H_rays=np.random.normal(size=(N,Q.n))
    else:
        N=H_rays.shape[0]
        assert H_rays.shape[1]==Q.n
    h_y=prog.NewContinuousVariables(2*N,1,"hy")
    H_y=np.vstack(( H_rays , -H_rays ))
    Y=pp.H_polytope(H_y,h_y)
    pp.subset(prog,Q,Y)
    prog.AddLinearCost(np.ones(2*N),np.array([0]),h_y)
    result=gurobi_solver.Solve(prog,None,None)
    if result.is_success():
        h_y_n=result.GetSolution(h_y).reshape(2*N,1)
        return pp.H_polytope(H_y,h_y_n)
    else:
        print("The polytope you gave me seems unbounded or \
              there is another error")

def ray_shooting_hyperplanes_old(Q,N=0,H_y=None):
    """
    Ray Shooting to find an outer-approximation of the AH-polytope
    """
    prog=MP.MathematicalProgram()
    Q=pp.to_AH_polytope(Q)
    if type(H_y)==type(None):
        if N==0:
            N=2**(Q.n-1) # This many hyperplanes I want :) 
        H_rays=np.random.normal(size=(N,Q.n))
        H_y=np.vstack(( H_rays , -H_rays ))
    else:
        N=H_rays.shape[0]
        assert H_rays.shape[1]==Q.n
    h_n=np.zeros((2*N,1))
    zeta=prog.NewContinuousVariables(Q.P.H.shape[1],1,"zeta")
    prog.AddLinearConstraint(A=Q.P.H,ub=Q.P.h,lb=-np.inf*np.ones((Q.P.h.shape[0],1)),vars=zeta)
    a=np.dot(H_rays[0,:],Q.T)
    b=np.dot(H_rays[0,:],Q.t)
    cost=prog.AddLinearCost(a,b,zeta)
    for i in range(2*N):
        new_a=np.dot(H_y[i,:],Q.T)
        new_b=np.dot(H_y[i,:],Q.t)
        cost.evaluator().UpdateCoefficients( -new_a, new_b)
        result=gurobi_solver.Solve(prog,None,None)
        if result.is_success():
            _s=result.GetSolution(zeta)
            h_n[i,0]=np.dot(new_a,_s)+new_b
        else:
            print("The polytope you gave me seems unbounded or \
                  there is another error")    
    return pp.H_polytope(H_y,h_n)

def ray_shooting_hyperplanes(Q,N=0,H_y=None,tol=1e-2):
    """
    Ray Shooting to find an outer-approximation of the AH-polytope
    """
    prog=MP.MathematicalProgram()
    Q=pp.to_AH_polytope(Q)
    if type(H_y)==type(None):
        if N==0:
            N=2**(Q.n) # This many hyperplanes I want :) 
        H_rays=np.random.normal(size=(N,Q.n))
        H_y=np.vstack(( H_rays , -H_rays ))
    else:
        assert H_y.shape[1]==Q.n
    # h_n=np.zeros((2*N,1))
    zeta=prog.NewContinuousVariables(Q.P.H.shape[1],1,"zeta")
    prog.AddLinearConstraint(A=Q.P.H,ub=Q.P.h,lb=-np.inf*np.ones((Q.P.h.shape[0],1)),vars=zeta)
    a=np.dot(H_y[0,:],Q.T)
    b=np.dot(H_y[0,:],Q.t)
    cost=prog.AddLinearCost(a,b,zeta)
    Y=pp.H_polytope(np.zeros((0,Q.n)),np.zeros((0,1)))
    for i in range(H_y.shape[0]):
        new_a=np.dot(H_y[i,:],Q.T)
        new_b=np.dot(H_y[i,:],Q.t)
        cost.evaluator().UpdateCoefficients( -new_a, new_b)
        result=gurobi_solver.Solve(prog,None,None)
        if result.is_success():
            _s=result.GetSolution(zeta)
            new_h=np.dot(new_a,_s)+new_b
            new_H=H_y[i,:].reshape(1,Q.n)
            if not _check_if_new_hyperplane_is_redundant(Y,new_H,new_h,tol=tol):
                Y.H=np.vstack((  Y.H,new_H ))
                Y.h=np.vstack((  Y.h,new_h ))
        else:
            print("The polytope you gave me seems unbounded or \
                  there is another error")    
    return Y

def _check_if_new_hyperplane_is_redundant(P,new_H,new_h,tol=1e-2):
    """
    P: H_polytope
    """
    assert P.type=='H_polytope'
    prog=MP.MathematicalProgram()
    q,n=P.H.shape
    if q==0:
        return False
    zeta=prog.NewContinuousVariables(n,1,"zeta")
    prog.AddLinearConstraint(A=P.H,ub=P.h,lb=-np.inf*np.ones((q,1)),vars=zeta)
    prog.AddLinearCost(-new_H.reshape(n),np.array([new_h]),zeta)
    result=gurobi_solver.Solve(prog,None,None)
    if result.is_success():
        _s=result.GetSolution(zeta)
#        print( np.asscalar( np.dot(new_H,_s) ) )
#        print( new_h )
#        print(np.asscalar( np.dot(new_H,_s) ) <=new_h + tol)
#        print("*")
        return np.asscalar( np.dot(new_H,_s) )  <= new_h + tol
    else:
        return False
    
def inner_optimization(Q,X=None,N=100,k=-1,method="alternate",iterations=5,tol=1e-2):
    """
    Q= AH_polytope
    X= H_polytope Candidate
    """
    # Sanity Checks
    assert Q.type=='AH_polytope' or Q.type=='V_polytope'
    Q=pp.to_AH_polytope(Q)
    if type(X)==type(None):
        X=ray_shooting_hyperplanes(Q,N=N,tol=tol)
    else:
        assert X.type=='H_polytope'
    # Program
    n=Q.n
    prog=MP.MathematicalProgram()
    # T=prog.NewSymmetricContinuousVariables(Q.n,'T')
    T=prog.NewContinuousVariables(n,n,"T")
    t=prog.NewContinuousVariables(n,1,"t")
#    prog.AddPositiveSemidefiniteConstraint(T)
    Y=pp.AH_polytope(T=T,t=t,P=X) 
    pp.subset(prog,Y,Q,k=k,verbose=True)
    if method=="SDP":
        prog.AddMaximizeLogDeterminantSymmetricMatrixCost(T)
        result=scs_solver.Solve(prog,None,None)
    elif method=="alternate":
        result=volume_maximization(prog, T,np.eye(n) ,iterations)
    if result.is_success():
        print("success")
        T_n= result.GetSolution(T)
        t_n= result.GetSolution(t).reshape(n,1)
        print("determinent=",np.linalg.det(T_n))
        return pp.affine_map( T=T_n, P=X, t=t_n),np.linalg.det(T_n)
    else:
        print("not succesfull") 
        
def outer_optimization(Q,X=None,N=100,k=-1):
    """
    Q= AH_polytope
    X= H_polytope Candidate
    """
    # Sanity Checks
    assert Q.type=='AH_polytope' or Q.type=='V_polytope'
    Q=pp.to_AH_polytope(Q)
    if type(X)==type(None):
        X=ray_shooting_hyperplanes(Q,N=N)
    else:
        assert X.type=='H_polytope'
    # Program
    n=Q.n
    prog=MP.MathematicalProgram()
#    T=prog.NewSymmetricContinuousVariables(Q.n,'T')
    T=prog.NewContinuousVariables(n,n,"T")
    t=prog.NewContinuousVariables(n,1,"t")
#    prog.AddPositiveSemidefiniteConstraint(T)
    prog.AddMaximizeLogDeterminantSymmetricMatrixCost(T)
    Q_new=pp.AH_polytope(T=np.dot(T,Q.T),t=np.dot(T,Q.t)+t,P=Q.P) 
    pp.subset(prog,Q_new,X,k=k,verbose=True)
    result=scs_solver.Solve(prog,None,None)
    if result.is_success():
        print("success")
        T_n= result.GetSolution(T)
        t_n= result.GetSolution(t).reshape(n,1)
        Tinv=np.linalg.inv(T_n)
        t_new=np.dot(-Tinv,t_n)
        print("determinent=",np.linalg.det(T_n))
#        Q_new_n=pp.AH_polytope(T=np.dot(T_n,Q.T),t=np.dot(T_n,Q.t)+t_n,P=Q.P) 
#        return Q_new_n
        return pp.affine_map( T=Tinv, P=X, t=t_new),np.linalg.det(Tinv)
#        return pp.AH_polytope(T=Tinv,t=t_new,P=X) 
    else:
        print("not succesfull") 
        
        
def inner_optimization_new(Q,X=None,N=100,k=-1):
    """
    Q= AH_polytope
    X= H_polytope Candidate
    """
    # Sanity Checks
    assert Q.type=='AH_polytope' or Q.type=='V_polytope'
    Q=pp.to_AH_polytope(Q)
    if type(X)==type(None):
        X=ray_shooting_hyperplanes(Q,N=N)
    else:
        assert X.type=='H_polytope'
    # Program
    n=Q.n
    prog=MP.MathematicalProgram()
    # T=prog.NewSymmetricContinuousVariables(Q.n,'T')
    T=prog.NewContinuousVariables(n,n,"T")
    t=prog.NewContinuousVariables(n,1,"t")
#    prog.AddPositiveSemidefiniteConstraint(T)
    Y=pp.AH_polytope(T=T,t=t,P=X) 
    pp.subset(prog,Y,Q,k=k,verbose=True)
    result=volume_maximization(prog, T,0.1*np.eye(n)+0.1 )
    if result.is_success():
        print("success")
        T_n= result.GetSolution(T)
        t_n= result.GetSolution(t).reshape(n,1)
        print("determinent=",np.linalg.det(T_n))
        return pp.affine_map( T=T_n, P=X, t=t_n),np.linalg.det(T_n)
    else:
        print("not succesfull")         
   
def volume_maximization(program,G,G_0,iterations=5,tol=0.01):
    x=G.reshape(-1)
    print("*"*20,"\n")
    print("\t\t\t Alternating Convex Program for Determinent Maximization")
    print("\t This is often faster than SDP solving. If not, use SDP")
    G_grad=_volume_gradient(G_0).reshape(-1)
    cost=program.AddLinearCost( -G_grad, np.array([0]), x )
    result=gurobi_solver.Solve(program,None,None)
    det=np.linalg.det(G_0)
    for i in range(iterations):
        if result.is_success():
            print(i,"det=",np.linalg.det(G_0))
            G_0=result.GetSolution(G)
            r=np.linalg.det(G_0)/det
            if r<1+tol and r>1-tol:
                print('converged')
                return result
            else:
                det=np.linalg.det(G_0)
                # print(G_0)
                new_a=_volume_gradient(G_0).reshape(-1)
                cost.evaluator().UpdateCoefficients( -new_a, np.array([0]))
                result=gurobi_solver.Solve(program,None,None)
    print('Error of convergence')
    return 
        

def _volume_gradient(G):
    S=combinations(range(G.shape[1]),G.shape[0])
    V_dot=np.zeros(G.shape)
    for s in S:
        Gs=np.hstack([G[:,i:i+1] for i in s])
        D=np.linalg.det(Gs)
        if D!=0:
            adj_Gs=D*np.linalg.inv(Gs)
            X=adj_Gs.T*np.sign(np.linalg.det(Gs))
        else:
#                print("Warning: we have a singular matrix")
            e=np.eye(G.shape[0])*10**(-5)
            adj_Gs=np.linalg.det(Gs+e)*np.linalg.inv(Gs+e)
            X=adj_Gs.T*np.sign(np.linalg.det(Gs+e))
#                X=adj_Gs.T*(np.linalg.det(Gs+e))
        for i in range(len(s)):
            V_dot[:,s[i]]+=X[:,i]
    return V_dot