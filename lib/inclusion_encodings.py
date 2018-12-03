#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 23 16:03:39 2018

@author: sadra
"""

# External imports:
import numpy as np
from gurobipy import Model,GRB,LinExpr,QuadExpr

def subset_LP(model,x,G,P,S):
    """
    Description: Add Farkas lemma constraints for subset inclusion of x+GP subset S
    Inputs: 
        model: Gurobi optimization model
        G: n * n_g generator matrix
        P:  primitive polytope
        S: polytope
        x: shift vector for GP
    Output:
        no direct output. Adds constraints to the model. 
    """
    (n,n_g)=G.shape
    (n_p,n_g)=P.H.shape
    assert(G.shape[1]==P.H.shape[1])
    (n_h,n)=S.H.shape
    assert(S.H.shape[1]==G.shape[0])
    Lambda=np.empty((n_h,n_p),dtype='object')
    for row in range(n_h):
        for column in range(n_p):
            Lambda[row,column]=model.addVar(lb=0)
    model.update()
    # Lambda * Pi = H * G
    for row in range(n_h):
        for column in range(n_g):
            s_left=LinExpr()
            s_right=LinExpr()
            for k in range(n_p):
                s_left.add(Lambda[row,k]*P.H[k,column])
            for k in range(n):
                s_right.add(S.H[row,k]*G[k,column])
            model.addConstr(s_left==s_right)
    # Lambda * P.h <= S.h - S.H*x
    for row in range(n_h):
        s_left=LinExpr()
        s_right=LinExpr()
        for k in range(n_p):
            s_left.add(Lambda[row,k]*P.h[k,0])
        for k in range(n):
            s_right.add(S.H[row,k]*x[k,0])
        model.addConstr(s_left<=S.h[row,0]-s_right) 

def subset_both_projection(model,x_l,G_l,P_l,x_r,G_r,P_r):
    alpha=np.empty((G_r.shape[1],G_l.shape[1]),dtype='object')
    alpha=add_Var_matrix(model,alpha)    
    beta=np.empty((G_r.shape[1],1),dtype='object')
    beta=add_Var_matrix(model,beta)
    Lambda=np.empty((P_r.H.shape[0],P_l.H.shape[0]),dtype='object')
    Lambda=add_Var_matrix(model,Lambda,pos=1)
    d=np.empty((G_r.shape[1],1),dtype='object')
    d=add_Var_matrix(model,d)
    model.update()
    constraints_AB_eq_CD(model,np.eye(G_l.shape[1]),G_l,G_r,alpha) # X=Y*alpha
    for row in range(G_r.shape[0]):
        model.addConstr(d[row,0]==x_r[row,0]-x_l[row,0])
    constraints_AB_eq_CD(model,np.eye(G_l.shape[1]),d,G_r,beta) # d=Y*beta
    # Now the main constraint!
    constraints_AB_eq_CD(model,Lambda,P_l.H,P_r.H,alpha)
    for row in range(Lambda.shape[0]):
        s_left=LinExpr()
        s_right=LinExpr()
        for k in range(Lambda.shape[1]):
            s_left.add(Lambda[row,k]*P_l.h[k,0])
        for k in range(beta.shape[0]):
            s_right.add(P_r.H[row,k]*beta[k,0])
        model.addConstr(s_left<=P_r.h[row,0]+s_right)
    return (alpha,beta,Lambda)
        
    
    
def subset_minkowski_right(model,x,G,P,list_of_polytopes):
    """
    Description: Add inclusion constraints for subset inclusion of x+GP subset S
    Inputs: 
        model: Gurobi optimization model
        G: n * n_g generator matrix
        P:  primitive polytope
        list_of_polytopes: polytope
        x: shift vector for GP
    Output:
        no direct output. Adds constraints to the model. 
    """
    (n,n_g)=G.shape
    (n_p,n)=P.H.shape
    T={}
    Lambda={}
    t={}
    for poly in list_of_polytopes:
        (n_h,n)=poly.H.shape
        Lambda[poly]=np.empty((n_h,n_p),dtype='object')
        T[poly]=np.empty((n,n_g),dtype='object')
        t[poly]=np.empty((n,1),dtype='object')
        for row in range(n_h):
            for column in range(n_p):
                Lambda[poly][row,column]=model.addVar(lb=0)
        for row in range(n):
            for column in range(n_g):
                T[poly][row,column]=model.addVar(lb=-GRB.INFINITY,ub=GRB.INFINITY)
        for row in range(n):
            t[poly][row,0]=model.addVar(lb=-GRB.INFINITY,ub=GRB.INFINITY)
    model.update()
    t_sum={}
    T_sum={}
    for row in range(n):
        for column in range(n_g):
            T_sum[row,column]=LinExpr()
        t_sum[row,0]=LinExpr()
    for poly in list_of_polytopes:
        # Lambda P.H= poly.H T[poly]
        (n_h,n)=poly.H.shape
        for row in range(n_h):
            for column in range(n_g):
                s_left=LinExpr()
                s_right=LinExpr()
                for k in range(n_p):
                    s_left.add(Lambda[poly][row,k]*P.H[k,column])
                for k in range(n):
                    s_right.add(poly.H[row,k]*T[poly][k,column])
                model.addConstr(s_left==s_right)
        # Lambda * P.h <= poly.h - poly.H * t[poly]
        for row in range(n_h):
            s_left=LinExpr()
            s_right=LinExpr()
            for k in range(n_p):
                s_left.add(Lambda[poly][row,k]*P.h[k,0])
            for k in range(n):
                s_right.add(poly.H[row,k]*t[poly][k,0])
            model.addConstr(s_left<=poly.h[row,0]-s_right)
        for row in range(n):
            for column in range(n_g):
                T_sum[row,column].add(T[poly][row,column])
            t_sum[row,0].add(t[poly][row,0])
    for row in range(n):
        for column in range(n_g):
            model.addConstr(T_sum[row,column]==G[row,column])
        model.addConstr(t_sum[row,0]==x[row,0])
    return (t,T)
    

def subset_zonotope_both(model,x,G,y,Z,c=1):
    """
    Description: Add inclusion constraints for subset inclusion of <x,G> subset <y,Z>
    Inputs: 
        model: Gurobi optimization model
        variable: x: n * 1 shift vector
        variable: G: n * n_g zonotope generator matrix
        variable: y: n * 1 shift vector 
        constant: Z: n * n_Z generator matrix of zonotope
    Output:
        no direct output. Adds constraints to the model.    
    """
    (n,n_G)=G.shape
    (n,n_Z)=Z.shape
    alpha=np.empty((n_Z,n_G),dtype='object')
    beta=np.empty((n_Z,1),dtype='object')
    d=np.empty((n,1),dtype='object')
    alpha_abs=np.empty(alpha.shape,dtype='object')
    beta_abs=np.empty(beta.shape,dtype='object')
    alpha=add_Var_matrix(model,alpha)
    beta=add_Var_matrix(model,beta)
    alpha_abs=add_Var_matrix(model,alpha_abs)
    beta_abs=add_Var_matrix(model,beta_abs)
    d=add_Var_matrix(model,d)
    model.update()
    constraints_AB_eq_CD(model,np.eye(n),G,Z,alpha)
    constraints_AB_eq_CD(model,np.eye(n),d,Z,beta)
    for row in range(n):
        model.addConstr(d[row,0]==x[row,0]-y[row,0])
    for row in range(n_Z):
        for column in range(n_G):
            model.addConstr(alpha_abs[row,column]>=alpha[row,column])
            model.addConstr(alpha_abs[row,column]>=-alpha[row,column])
    for row in range(n_Z):
        model.addConstr(beta_abs[row,0]>=beta[row,0])
        model.addConstr(beta_abs[row,0]>=-beta[row,0])
    for row in range(n_Z):
        sum_terms=LinExpr()
        sum_terms.add(beta_abs[row,0])
        for column in range(n_G):
            sum_terms.add(alpha_abs[row,column])
        model.addConstr(sum_terms<=c)
    return (alpha,beta)
    
def subset_zonotope_convexhullOFzonotopes(model,x,G,list_of_zonotopes):
    """
    a zonotope being inside the convexhull of a family of zonotopes
    Note: this is only a sufficient condition, not necessay
    WARNING: NOT WORKING YET
    """
    (n,n_G)=G.shape
    alpha={}
    alpha_abs={}
    beta={}
    beta_abs={}
    Lambda={}
    F={}
    d={}
    y={}
    for zono in list_of_zonotopes:
        (n,n_Z)=zono.shape
        F[zono]=np.empty((n,n_G),dtype='object')
        alpha[zono]=np.empty((n_Z,n_G),dtype='object')
        beta[zono]=np.empty((n_Z,1),dtype='object')
        y[zono]=np.empty((n,1),dtype='object')
        alpha_abs[zono]=np.empty(alpha.shape,dtype='object')
        beta_abs[zono]=np.empty(beta.shape,dtype='object')
        # Add Matrices
        F[zono]=add_Var_matrix(model,F[zono])
        alpha[zono]=add_Var_matrix(model,alpha[zono])
        beta[zono]=add_Var_matrix(model,beta[zono])
        alpha_abs[zono]=add_Var_matrix(model,alpha_abs[zono])
        beta_abs[zono]=add_Var_matrix(model,beta_abs[zono])
        y[zono]=add_Var_matrix(model,d[zono])
        Lambda[zono]=np.empty((1,1),dtype='object')
        Lambda[zono][0,0]=model.addVar(lb=0,ub=1)
        # Model Update
        model.update()
        # Constraints on 
        constraints_AB_eq_CD(model,np.eye(n),F[zono],zono.G,alpha)
        constraints_AB_eq_CD(model,np.eye(n),y[zono],zono.G,beta)
        for row in range(n_Z):
            for column in range(n_G):
                model.addConstr(alpha_abs[zono][row,column]>=alpha[zono][row,column])
                model.addConstr(alpha_abs[zono][row,column]>=-alpha[zono][row,column])
        for row in range(n_Z):
            model.addConstr(beta_abs[zono][row,0]>=beta[zono][row,0])
            model.addConstr(beta_abs[zono][row,0]>=-beta[zono][row,0])
        for row in range(n_Z):
            sum_terms=LinExpr()
            sum_terms.add(beta_abs[zono][row,0])
            for column in range(n_G):
                sum_terms.add(alpha_abs[zono][row,column])
            model.addConstr(sum_terms<=Lambda[zono])
    # Summation for F:
    constraints_sum(model,G,[F[zono] for zono in list_of_zonotopes])
    constraints_sum(model,x,[d[zono] for zono in list_of_zonotopes])
    constraints_sum(model,np.array([1]).reshape(1,1),[Lambda[zono] for zono in list_of_zonotopes])
    for row in range(n):
        model.addConstr(d[row,0]==x[row,0]-y[row,0])  
    pass
    

    
    
def Q_inside_conv_P_s(s,model,x,G,Q,list_of_polytopes):
    # WARNING: INCOMPLETE CODE
    """
    WARNING: INCOMPLETE CODE
    Given:
        s: PWA system that is considered
        Numbers:
            Q= {HQ x \le Hq} \subset R^q
            P_i={H_i x \le H_i} \subset R^n, i=1,...,N, 
        Symbols:
            model= Gurobi model
            T= matrix in R^{n*q}: model variable
            d= vector in R^n: model variable
    
    What the function does:
        Model with Constraints required to have TQ+d in Convexhull(Ps)
    
    Output:
        None
    """
    pass
    Lambda={}
    delta={}
    for polytope in list_of_polytopes:
        (n_H,n_HQ)=(polytope.H.shape[0],Q.H.shape[0])
        delta[polytope]=model.addVar(lb=0,ub=1)
        Lambda[polytope]=np.empty((n_H,n_Q))
        
    (n,n_g)=G.shape
    (n_p,n)=s.Pi.shape
    (n_h,n)=(2*s.n,s.n)
    z_pol={}
    x_pol={}
    G_pol={}
    G_bound=100
    for polytope in list_of_polytopes:
        Lambda[polytope]=np.empty((n_h,n_p),dtype='object')
        x_pol[polytope]=np.empty((n,1),dtype='object')
        G_pol[polytope]=np.empty((n,n_g),dtype='object')
        for row in range(n_h):
            for column in range(n_p):
                Lambda[polytope][row,column]=model.addVar(lb=0)
        z_pol[polytope]=model.addVar(vtype=GRB.BINARY)
        for row in range(n):
            x_pol[polytope][row,0]=model.addVar(lb=-G_bound,ub=G_bound)
        for row in range(n):
            for column in range(n_g):
                G_pol[polytope][row,column]=model.addVar(lb=-G_bound,ub=G_bound)                
    model.update()
    z_sum=LinExpr()
    G_sum=np.empty((n,n_g),dtype='object')
    x_sum=np.empty((n,1),dtype='object')
    for row in range(n):
        x_sum[row,0]=LinExpr()
        for column in range(n):
            G_sum[row,column]=LinExpr()
    for polytope in list_of_polytopes:
        z_sum.add(z_pol[polytope])
        for row in range(n):
            x_sum[row,0].add(x_pol[polytope][row,0])
            for column in range(n_g):
                G_sum[row,column].add(G_pol[polytope][row,column])        
        H=np.dot(s.Pi,polytope.G_eps_inv)
        h=np.ones((s.Pi.shape[0],1))+np.dot(H,polytope.x)
        factorize=np.amax(abs(H),1).reshape(s.Pi.shape[0],1)
        H=np.divide(H,factorize)
        h=np.divide(h,factorize)
        for row in range(n_h):
            for column in range(n_g):
                s_left=LinExpr()
                s_right=LinExpr()
                for k in range(n_p):
                    s_left.add(Lambda[polytope][row,k]*s.Pi[k,column])
                for k in range(n):
                    s_right.add(H[row,k]*G_pol[polytope][k,column])
                model.addConstr(s_left==s_right)
        # Lambda * 1 <= H*
        for row in range(n_h):
            s_left=LinExpr()
            s_right=LinExpr()
            for k in range(n_p):
                s_left.add(Lambda[polytope][row,k])
            for k in range(n):
                s_right.add(H[row,k]*x_pol[polytope][k,0])
            model.addConstr(s_left<=h[row,0]*z_pol[polytope]-s_right)
    model.addConstr(z_sum==1)
    for row in range(n):
        model.addConstr(x_sum[row,0]==x[row,0])
        for column in range(n_g):
            model.addConstr(G_sum[row,column]==G[row,column])
    return z_pol

def PI(n):
    return np.vstack((np.eye(n),-np.eye(n)))
    
def constraints_AB_eq_CD(model,A,B,C,D):
    """
    Add constraint A*B=C*D to gurobi model
    """
    if B.shape[1]!=C.shape[0] or A.shape[1]!=B.shape[0] or C.shape[1]!=D.shape[0] or A.shape[0]!=C.shape[0]:
        ValueError("Dimensions mistmatch")
    for row in range(A.shape[0]):
        for column in range(B.shape[1]):
            lhs=LinExpr()
            rhs=LinExpr()
            for k in range(A.shape[1]):
                lhs.add(A[row,k]*B[k,column])
            for k in range(C.shape[1]):
                rhs.add(C[row,k]*D[k,column])
            model.addConstr(rhs==lhs)
            
def add_Var_matrix(model,A,pos=0,delta=None):
    for row in range(A.shape[0]):
        for column in range(A.shape[1]):
            if pos==0:
                A[row,column]=model.addVar(lb=-GRB.INFINITY,ub=GRB.INFINITY)
            elif pos==1:
                A[row,column]=model.addVar(lb=0,ub=GRB.INFINITY)
            if delta!=None:
                A[row,column]=model.addVar(lb=-delta,ub=delta)
    model.update()
    return A

    
def constraints_sum(model,A,list_of_A):
    """
    Simple Summation
    """
    for row in range(A.shape[0]):
        for column in range(A.shape[1]):
            sum_term=LinExpr()
            for a in list_of_A:
                sum_term.add(a[row,column])
            model.addConstr(A[row,column]==sum_term)            