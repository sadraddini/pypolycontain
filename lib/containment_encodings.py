# -*- coding: utf-8 -*-
"""
Created on Mon Dec  3 09:47:34 2018

@author: sadra
"""

# External imports:
import numpy as np
from gurobipy import Model,GRB,LinExpr,QuadExpr

from pypolycontain.lib.AH_polytope import to_AH_polytope

def subset_generic(model,Q1,Q2):
    """
    Adds containment property Q1 subset Q2
    Inputs:
        Q1,Q2: either polytope, zonotope, or AH_polytope
    Output:
        No direct output, adds Q1 \subset Q2 to the model
    """
    Q1=to_AH_polytope(Q1)
    Q2=to_AH_polytope(Q2)
    Gamma=tupledict_to_array(model.addVars(range(Q2.T.shape[1]),range(Q1.T.shape[1]),lb=-GRB.INFINITY,ub=GRB.INFINITY,name="Gamma"))
    Lambda=tupledict_to_array(model.addVars(range(Q2.P.H.shape[0]),range(Q1.P.H.shape[0]),lb=0,ub=GRB.INFINITY,name="Lambda"))
    beta=tupledict_to_array(model.addVars(range(Q2.T.shape[1]),[0],lb=-GRB.INFINITY,ub=GRB.INFINITY,name="beta"))
    model.update()
    n=Q1.T.shape[0]
    assert n==Q2.T.shape[0]
    constraints_list_of_tuples(model,[(np.eye(Q1.T.shape[1]),Q1.T),(-Q2.T,Gamma)],sign="=")
    constraints_list_of_tuples(model,[(np.eye(n),Q2.t),(-np.eye(n),Q1.t),(-Q2.T,beta)],sign="=")
    constraints_list_of_tuples(model,[(Lambda,Q1.P.H),(-Q2.P.H,Gamma)],sign="=")
    constraints_list_of_tuples(model,[(Lambda,Q1.P.h),(-np.eye(Q2.P.h.shape[0]),Q2.P.h),(-Q2.P.H,beta)],sign="=")
    
    
    
    

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
    """
    Description: Add polytope containment constraints for subset inclusion of x_l+G_l P_l subset x_r+G_r P_r
    Inputs: 
        model: Gurobi optimization model
        x: vector, G: matrix, P: polytopes
    Output:
        no direct output. Adds constraints to the model. 
    Method:
        matrices alpha, beta, Lambda > 0
        add constraints
            x_r-x_l=G_r*beta
            G_l=G_r*alpha
            Lambda * H_l = H_r * alpha
            Lambda * h_l <= h_r + H_r * beta
    """
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

def subset_LP_disjunctive(model,x,G,P,list_of_polytopes):
    """
    Description: Add disjucntive x+GP subset list_of_polytopes
    Inputs: 
        model: Gurobi optimization model
        G: n * n_g generator matrix
        P:  primitive polytope
        x: shift vector for GP
    Output:
        no direct output. Adds constraints to the model. 
    """
    Lambda={}
    beta={}
    Delta={}
    X={}
    for poly in list_of_polytopes:
        Lambda[poly]=add_Var_matrix(model,(poly.H.shape[0],P.H.shape[1]),pos=1)
        beta[poly]=add_Var_matrix(model,(poly.H.shape[0],1))
        Delta[poly]=model.addVar(vtype=GRB.BINARY)
        X[poly]=add_Var_matrix(model,G.shape)
    model.update()
    for poly in list_of_polytopes:
        constraints_AB_eq_CD(model,Lambda[poly],P.H,poly.H,X[poly])
        for row in range(Lambda.shape[0]):
            s_left=LinExpr()
            s_right=LinExpr()
            for k in range(Lambda.shape[1]):
                s_left.add(Lambda[row,k]*P.h[k,0])
            for k in range(beta.shape[0]):
                s_right.add(-poly.H[row,k]*beta[poly][k,0])
            model.addConstr(s_left<=Delta[poly]*poly.h[row,0]+s_right)
    R=LinExpr()
    for poly in list_of_polytopes:
        R.add(Delta[poly])
    model.addConstr(R==1)
    constraints_sum(model,G,[X[poly] for poly in list_of_polytopes])
    constraints_sum(model,x,[beta[poly] for poly in list_of_polytopes])


def subset_zonotopes(model,z_l,z_r):
    """
    Description: Add inclusion constraints for subset inclusion of <x,G> subset <y,Z>
    Inputs: 
        model: Gurobi optimization model
        z_l zonotope on the left
        z_r zonotope on the left
    Output:
        no direct output. Adds constraints to the model.  
    Method:
        find matrix alpha and beta such that
            G_x = G_y * alpha
            x-y = g_y * beta  
            ||(alpha,beta)||_infty<=1
    """
    (n,n_l)=z_l.G.shape
    (n,n_r)=z_r.G.shape
    Gamma=add_Var_matrix(model,(n_r,n_l))
    beta=add_Var_matrix(model,(n_r,1))
    d=add_Var_matrix(model,(n,1))
    model.update()
    for row in range(n):
        model.addConstr(d[row,0]==z_l.x[row,0]-z_r.x[row,0])
    constraints_AB_eq_CD(model,np.eye(n),z_l.G,z_r.G,Gamma)
    constraints_AB_eq_CD(model,np.eye(n),d,z_r.G,beta)
    infinity_norm_list(model,[Gamma,beta],1)
    return (Gamma,beta,d)


def subset_zonotopes_convexhull(model,x,G,list_of_zonotopes):
    """
    Description: Add inclusion constraints for subset inclusion of <x,G> subset convexhull of a list of zonotopes
    Inputs: 
        model: Gurobi optimization model
        variable: x : n * 1 shift vector
        variable: G : n * n_g zonotope generator matrix
    Output:
        no direct output. Adds constraints to the model.  
    Method:
        find matrix Gamma_i and beta_i such that
            G = sum over G_i * Gamma_i
            x = sum over G_i * beta_i + lambda_i * x_i            
    """    
    Gamma={}
    Lambda={}
    beta={}
    for zono in list_of_zonotopes:
        Lambda[zono]=model.addVar(lb=0)
        Gamma[zono]=add_Var_matrix(model,(zono.G.shape[1],G.shape[1]))
        beta[zono]=add_Var_matrix(model,(zono.G.shape[1],1))
    model.update()
    for zono in list_of_zonotopes:
        infinity_norm_list(model,[Gamma[zono],beta[zono]],Lambda[zono])
    for row in range(G.shape[0]):
        for column in range(G.shape[1]):
            R=LinExpr()
            for zono in list_of_zonotopes:
                for k in range(zono.G.shape[1]):
                    R.add(zono.G[row,k]*Gamma[zono][k,column])
            model.addConstr(R==G[row,column])
    for row in range(G.shape[0]):
        R=LinExpr()
        for zono in list_of_zonotopes:
            for k in range(zono.G.shape[1]):
                R.add(zono.G[row,k]*beta[zono][k,0])
            R.add(Lambda[zono]*zono.x[row,0])
        model.addConstr(R==x[row,0])
    R=LinExpr()
    for zono in list_of_zonotopes:
        R.add(Lambda[zono])
    model.addConstr(R==1)
    return (Lambda,beta,Gamma)
        
def subset_zonotopes_disjunctive(model,x,G,list_of_zonotopes):
    """
    Description: Add inclusion constraints for subset inclusion of <x,G> in disjunction of zonotopes
    Inputs: 
        model: Gurobi optimization model
        variable: x : n * 1 shift vector
        variable: G : n * n_g zonotope generator matrix
    Output:
        no direct output. Adds constraints to the model.  
    Method:
        find matrix Gamma_i and beta_i such that
            G = sum over G_i * Gamma_i
            x = sum over G_i * beta_i + lambda_i * x_i            
    """    
    Gamma={}
    Delta={}
    beta={}
    for zono in list_of_zonotopes:
        Delta[zono]=model.addVar(vtype=GRB.BINARY)
        Gamma[zono]=add_Var_matrix(model,(zono.G.shape[1],G.shape[1]))
        beta[zono]=add_Var_matrix(model,(zono.G.shape[1],1))
    model.update()
    for zono in list_of_zonotopes:
        infinity_norm_list(model,[Gamma[zono],beta[zono]],Delta[zono])
    for row in range(G.shape[0]):
        for column in range(G.shape[1]):
            R=LinExpr()
            for zono in list_of_zonotopes:
                for k in range(zono.G.shape[1]):
                    R.add(zono.G[row,k]*Gamma[zono][k,column])
            model.addConstr(R==G[row,column])
    for row in range(G.shape[0]):
        R=LinExpr()
        for zono in list_of_zonotopes:
            for k in range(zono.G.shape[1]):
                R.add(zono.G[row,k]*beta[zono][k,0])
            R.add(Delta[zono]*zono.x[row,0])
        model.addConstr(R==x[row,0])
    R=LinExpr()
    for zono in list_of_zonotopes:
        R.add(Delta[zono])
    model.addConstr(R==1)
    return (Delta,beta,Gamma)










"""
Auxilary Gurobi Shortcut Functions
used for convenience
"""

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
            
def add_Var_matrix(model,size,pos=0,delta=None):
    A=np.empty(size,dtype='object')
    for row in range(size[0]):
        for column in range(size[1]):
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
            
            
def infinity_norm(model,A_abs,alpha):
    for row in range(A_abs.shape[0]):
        lin=LinExpr()
        for column in range(A_abs.shape[1]):
            lin.add(A_abs[row,column])
        model.addConstr(lin<=alpha)

def absolute_value(model,A,A_abs):
    for row in range(A_abs.shape[0]):
        for column in range(A_abs.shape[1]):
            model.addConstr(A[row,column]<=A_abs[row,column])
            model.addConstr(-A[row,column]<=A_abs[row,column])

def infinity_norm_list(model,list_of_matrices,alpha):
    """
    inputs:
        list_of_matrices=A_0,A_1,...A_N, stacked horizontally (same # of rows)
        alpha: scalar or variable
    Output:
        adds ||A_0,A_1,...A_N||_infty<=alpha to the model
    """
    A_abs={}
    for i in range(len(list_of_matrices)):
        A=list_of_matrices[i]
        A_abs[i]=add_Var_matrix(model,A.shape)
        model.update()
        absolute_value(model,A,A_abs[i])
    # Now the infinity norm
    A0=list_of_matrices[0]
    for row in range(A0.shape[0]):
        lin=LinExpr()
        for i in range(len(list_of_matrices)):
            for column in range(list_of_matrices[i].shape[1]):
                lin.add(A_abs[i][row,column])
        model.addConstr(lin<=alpha)

            
def sum_matrix_equality(model,A,B,C):
    """
    add A=B+C
    """
    for row in range(A.shape[0]):
        for column in range(A.shape[1]):
            model.addConstr(A[row,column]==B[row,column]+C[row,column])
            
def tupledict_to_array(mytupledict):
    # It should be 2D
    n,m=max(mytupledict.keys())
    n+=1
    m+=1
    array=np.empty((n,m),dtype="object")
    for i in range(n):
        for j in range(m):
            array[i,j]=mytupledict[i,j]
    return array

def constraints_list_of_tuples(model,mylist,sign="="):
    term_0=mylist[0]
    ROWS,COLUMNS=term_0[0].shape[0],term_0[1].shape[1]
    for row in range(ROWS):
        for column in range(COLUMNS):
            expr=LinExpr()
            for term in mylist:
                q,qp=term[0].shape[1],term[1].shape[0]
                if q!=qp:
                    raise ValueError(term,"q=%d qp=%d"%(q,qp))
                if type(term[1][0,column])==type(model.addVar()):
                    expr.add(LinExpr([(term[0][row,k],term[1][k,column]) for k in range(q)]))
                elif type(term[0][row,0])==type(model.addVar()):
                    expr.add(LinExpr([(term[1][k,column],term[0][row,k]) for k in range(q)]))
                else:
                    expr.addConstant(sum([term[1][k,column]*term[0][row,k] for k in range(q)]))
            if sign=="<":
                model.addConstr(expr<=0)
            elif sign=="=":
                model.addConstr(expr==0)
            elif sign==">=":
                model.addConstr(expr>=0)
            else:
                raise "sign indefinite"