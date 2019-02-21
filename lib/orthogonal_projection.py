# -*- coding: utf-8 -*-
"""
Created on Fri Jan  4 14:27:01 2019

@author: sadra


Orthogonal Projection
"""

import numpy as np
from gurobipy import Model,LinExpr,QuadExpr,GRB

from pypolycontain.utils.utils import valuation


class ortho_polytope:
    def __init__(self,H,F,g):
        self.H=H
        self.F=F
        self.g=g
    
def orthogonal_projection_fixed_Hx(Hx,H,F,g,xbar):
    """
    Problem: find the inner-approximation of 
    Inputs: q is the user-defined intger for the number of hyperplanes in X
    """
    model=Model("Orthogonal Projection with fixed H")
    Lambda_0=add_Var_matrix(model,(H.shape[0],Hx.shape[0]),pos=1)
    Gamma=add_Var_matrix(model,(F.shape[1],Hx.shape[1]))
    beta_u=add_Var_matrix(model,(F.shape[1],1))
    beta_x=add_Var_matrix(model,(Hx.shape[1],1))
    epsilon=add_Var_matrix(model,(1,1))
    n=H.shape[1]
    Hball=np.vstack((np.eye(n),-np.eye(n)))
    hball=np.ones((2*n,1))
    # Other side    
    X_1=add_Var_matrix(model,(Hx.shape[1],H.shape[1]))
    X_2=add_Var_matrix(model,(Hx.shape[1],F.shape[1]))
    Lambda_1=add_Var_matrix(model,(Hx.shape[0],H.shape[0]),pos=1)
    Lambda_2=add_Var_matrix(model,(Hball.shape[0],F.shape[0]),pos=1)
    model.update()
    constraints_list_of_tuples(model,[(Lambda_0,Hx),(-np.eye(H.shape[0]),H),(-F,Gamma)])
    constraints_list_of_tuples(model,[(Lambda_0,np.ones((Lambda_0.shape[1],1))),\
        (-np.eye(g.shape[0]),g),(H,xbar),(-F,beta_u)],sign="<")
    constraints_list_of_tuples(model,[(Lambda_1,H),(-Hx,X_1)])
    constraints_list_of_tuples(model,[(Lambda_1,F),(-Hx,X_2)])    
    constraints_list_of_tuples(model,[(Lambda_2,H),(np.eye(Hball.shape[0]),-Hball),(Hball,X_1)])
    constraints_list_of_tuples(model,[(Lambda_2,F),(Hball,X_2)])
    constraints_list_of_tuples(model,[(Lambda_1,g),(-np.eye(Lambda_1.shape[0]),\
        np.ones((Lambda_1.shape[0],1))),(Hx,beta_x)],sign="<")
    constraints_list_of_tuples(model,[(Lambda_2,g),(-hball,epsilon),(Hball,xbar),(-Hball,beta_x)],sign="<")    
    model.setObjective(epsilon[0,0])
    model.optimize()
    if model.Status!=2:
        return False,False,False,False
#    print "epsilon=",epsilon[0,0].X
#    print "Lambda_0",valuation(Lambda_0)
#    print "Lambda_1",valuation(Lambda_1)
#    print "Lambda_2",valuation(Lambda_2)
#    print "Gamma",valuation(Gamma)
#    print "beta_x",valuation(beta_x)
#    print "beta_u",valuation(beta_u)   
#    print "X_1",valuation(X_1)
#    print "X_2",valuation(X_2)  
    return valuation(Lambda_0),valuation(X_1),valuation(X_2),valuation(beta_x)


def orthogonal_projection_gradient_decent(Hx,Lambda_0,X_1,X_2,beta_x,H,F,g,xbar,delta):
    """
    Gradient Decent for H_x, X_1, X_2, Lambda_0
    """
    model=Model("Orthogonal Projection Projected Gradient Descent")
    delta_Lambda_0=add_Var_matrix(model,(H.shape[0],Hx.shape[0]),delta=delta)
    Gamma=add_Var_matrix(model,(F.shape[1],Hx.shape[1]))
    beta_u=add_Var_matrix(model,(F.shape[1],1))
    delta_beta_x=add_Var_matrix(model,(Hx.shape[1],1),delta=delta)
    epsilon=add_Var_matrix(model,(1,1))
    n=H.shape[1]
    Hball=np.vstack((np.eye(n),-np.eye(n)))
    hball=np.ones((2*n,1))  
    delta_Hx=add_Var_matrix(model,Hx.shape,delta=delta)
    # Other side
    delta_X_1=add_Var_matrix(model,(Hx.shape[1],H.shape[1]),delta=delta)
    delta_X_2=add_Var_matrix(model,(Hx.shape[1],F.shape[1]),delta=delta)
    Lambda_1=add_Var_matrix(model,(Hx.shape[0],H.shape[0]),pos=1)
    Lambda_2=add_Var_matrix(model,(Hball.shape[0],F.shape[0]),pos=1)
    model.update()    
    # ------ Constraints -------
    constraints_list_of_tuples(model,[(Lambda_0,Hx),(delta_Lambda_0,Hx),(Lambda_0,delta_Hx),(-np.eye(H.shape[0]),H),(-F,Gamma)])
    constraints_list_of_tuples(model,[(Lambda_0,np.ones((Lambda_0.shape[1],1))),(delta_Lambda_0,np.ones((Lambda_0.shape[1],1))),\
        (-np.eye(g.shape[0]),g),(H,xbar),(-F,beta_u)],sign="<")
    constraints_list_of_tuples(model,[(Lambda_1,H),(-Hx,X_1),(delta_Hx,-X_1),(-Hx,delta_X_1)])
    constraints_list_of_tuples(model,[(Lambda_1,F),(-Hx,X_2),(delta_Hx,-X_2),(-Hx,delta_X_2)])    
    constraints_list_of_tuples(model,[(Lambda_2,H),(np.eye(Hball.shape[0]),-Hball),(Hball,X_1),(Hball,delta_X_1)])
    constraints_list_of_tuples(model,[(Lambda_2,F),(Hball,X_2),(Hball,delta_X_2)])
    constraints_list_of_tuples(model,[(Lambda_1,g),(-np.eye(Lambda_1.shape[0]),\
        np.ones((Lambda_1.shape[0],1))),(Hx,beta_x),(delta_Hx,beta_x),(Hx,delta_beta_x)],sign="<")    
    constraints_list_of_tuples(model,[(Lambda_2,g),(-hball,epsilon),(Hball,xbar),(-Hball,beta_x),(-Hball,delta_beta_x)],sign="<") 
    # Lambda_0 should remain non-negative
#    print "*"*20,"sizes",Lambda_0.shape,delta,Lambda_0.shape
    constraints_list_of_tuples(model,[(-np.eye(Lambda_0.shape[0]),Lambda_0),(-np.eye(Lambda_0.shape[0]),delta_Lambda_0)],sign="<") 
    # Solve
    model.setObjective(epsilon[0,0])
    model.setParam("OutputFlag",False)
    model.optimize()
    if model.Status!=2:
        return False,False,False,False,False,False
#    print "epsilon=",epsilon[0,0].X
#    print "delta_Lambda_0",valuation(delta_Lambda_0)
#    print "delta_X_1",valuation(delta_X_1)
#    print "delta_X_2",valuation(delta_X_2)  
    return valuation(delta_Lambda_0),valuation(delta_X_1),valuation(delta_X_2),valuation(delta_beta_x),valuation(delta_Hx),epsilon[0,0].X

def gradient_decent(Hx,H,F,g,xbar,delta,N=10):
    history=[]
    for i in range(N):
        print "iteration","*"*30,i,"*"*30
        Lambda_0,X_1,X_2,beta_x=orthogonal_projection_fixed_Hx(Hx,H,F,g,xbar)
        if type(Lambda_0)==type(False):
            break
        delta_Lambda_0,delta_X_1,delta_X_2,delta_beta_X,delta_Hx,epsilon=orthogonal_projection_gradient_decent(Hx,Lambda_0,X_1,X_2,beta_x,H,F,g,xbar,delta)
        history.append((Hx,epsilon))
        print "epsilon","+"*30,"+"*30,epsilon
        if epsilon==False:
            break
        Lambda_0=Lambda_0+delta_Lambda_0
        X_1=X_1+delta_X_1
        X_2=X_2+delta_X_2
        beta_x=beta_x+delta_beta_X
        Hx=Hx+delta_Hx
    return history
    
    
    
    
    






"""
Auxilary Gurobi Shortcut Functions
used for convenience
"""

def constraints_list_of_tuples(model,mylist,sign="="):
    term_0=mylist[0]
    ROWS=term_0[0].shape[0]
    COLUMNS=term_0[1].shape[1]
    for row in range(ROWS):
        for column in range(COLUMNS):
            expr=LinExpr()
            i=0
            for term in mylist:
                i+=1
                q,qp=term[0].shape[1],term[1].shape[0]
                assert q==qp
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