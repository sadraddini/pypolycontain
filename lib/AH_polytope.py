#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 20 15:16:48 2019

@author: sadra
"""

import numpy as np

from scipy.optimize import linprog as lp
from gurobipy import Model,LinExpr,QuadExpr,GRB



class AH_polytope():
    """
    Affine Transformation of an H-polytope
    """
    def __init__(self,T,t,P):
        """
        Initilization: T,t,P
        """
        self.T=T # Matrix n*n_p
        self.t=t # vector n*1
        self.P=P # Polytope in n_p dimensions
        self.n=T.shape[0]
        if T.shape[1]!=P.H.shape[1]:
            ValueError("Error: not appropriate T size, it is",T.shape[1],P.n)
        self.type="AH_polytope"
        self.method="Gurobi"
    
    def __repr__(self):
        return "AH_polytope from R^%d to R^%d"%(self.P.n,self.n)
    
    def is_inside(self,x):
        """
        Return if x is inside AH_polytope by checking the feasibility of a linear program
        """
        if self.method=="Gurobi":
            model=Model("check_if_inside")
            p=tupledict_to_array(model.addVars(range(self.P.n),[0],name="p"))
            model.update()
            constraints_list_of_tuples(model,[(self.T,p),(np.eye(self.n),self.t-x)],sign="=")
            constraints_list_of_tuples(model,[(self.P.H,p),(-np.eye(self.P.h.shape[0]),self.P.h)],sign="<")
            model.setParam('OutputFlag', False)
            model.optimize()
            if model.Status==3:
                return False
            else:
                return True  
        elif self.method=="scipy":
            sol=lp(np.ones((self.P.n,1)), A_ub=self.P.H, b_ub=self.P.h, A_eq=self.T, b_eq=x-self.t)
    #        print sol.message
            return sol.status==0
        else:
            raise ValueError("method %s not recognized"%method)
            
            


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