#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 26 11:15:25 2019

@author: sadra

Computing Hausdoff Distance Between AH-Polytopes
"""

import numpy as np

from gurobipy import Model,LinExpr,QuadExpr,GRB

from pypolycontain.lib.AH_polytope import AH_polytope,to_AH_polytope



def Hausdorff_directed(Q1,Q2,ball="infinty_norm"):
    """
    Minimum epsilon such that Q1 \subset Q2+epsilon(Ball)
    zero if and only if Q1 subset Q2
    """
    Q1,Q2=to_AH_polytope(Q1),to_AH_polytope(Q2)
    model=Model("Hausdorff Distance")
    D=model.addVar(lb=0,obj=1)
    n=Q1.t.shape[0]
    assert n==Q2.t.shape[0]
    assert n==Q1.T.shape[0]
    assert n==Q2.T.shape[0]
    if ball=="infinty_norm":
        HB=np.vstack((np.eye(n),-np.eye(n)))
        hB=np.vstack((np.ones((n,1)),np.ones((n,1))))
    Lambda_1=tupledict_to_array(model.addVars(range(Q2.P.H.shape[0]),range(Q1.P.H.shape[0]),lb=0,ub=GRB.INFINITY,name="Lambda1"))
    Lambda_2=tupledict_to_array(model.addVars(range(HB.shape[0]),range(Q1.P.H.shape[0]),lb=0,ub=GRB.INFINITY,name="Lambda2"))
    Gamma_1=tupledict_to_array(model.addVars(range(Q2.P.H.shape[1]),range(Q1.P.H.shape[1]),lb=-GRB.INFINITY,ub=GRB.INFINITY,name="Gamma1"))
    Gamma_2=tupledict_to_array(model.addVars(range(HB.shape[1]),range(Q1.P.H.shape[1]),lb=-GRB.INFINITY,ub=GRB.INFINITY,name="Gamma2"))
    beta_1=tupledict_to_array(model.addVars(range(Q2.P.H.shape[1]),[0],lb=-GRB.INFINITY,ub=GRB.INFINITY,name="beta1"))
    beta_2=tupledict_to_array(model.addVars(range(HB.shape[1]),[0],lb=-GRB.INFINITY,ub=GRB.INFINITY,name="beta2"))
    model.update()
    constraints_list_of_tuples(model,[(Lambda_1,Q1.P.H),(-Q2.P.H,Gamma_1)],sign="=")
    constraints_list_of_tuples(model,[(Lambda_2,Q1.P.H),(-HB,Gamma_2)],sign="=")
    constraints_list_of_tuples(model,[(Lambda_1,Q1.P.h),(-Q2.P.H,beta_1),(-np.eye(Q2.P.h.shape[0]),Q2.P.h)],sign="<")
    constraints_list_of_tuples(model,[(Lambda_2,Q1.P.h),(-hB,np.array([D]).reshape(1,1)),(-HB,beta_2)],sign="<")
    # T Constraints come here!
    constraints_list_of_tuples(model,[(-Q2.T,beta_1),(np.eye(n),Q2.t-Q1.t),(-np.eye(n),beta_2)],sign="=")
    constraints_list_of_tuples(model,[(Q2.T,Gamma_1),(np.eye(n),Gamma_2),(-np.eye(n),Q1.T)],sign="=")
    # Optimize
    model.setParam("OutputFlag",False)
    model.optimize()
    return D.X
    
   
    















"""
Auxilary Gurobi Shortcut Functions
used for convenience
"""

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
                    raise ValueError(term,"q=%d qp=%d"%(q,qp),"*\n the constrainst where",mylist)
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