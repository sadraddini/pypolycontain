#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 10 13:52:06 2018

@author: sadra
"""
import numpy as np
from gurobipy import Model, GRB, LinExpr

def check_empty_polytope(p):
    model=Model("Row Redundancy Check")
    n=p.H.shape[1]
    x=np.empty((n,1),dtype='object')
    for row in range(n):
        x[row,0]=model.addVar(lb=-GRB.INFINITY,ub=GRB.INFINITY)
    model.update()
    for row in range(p.H.shape[0]):
        Hx=LinExpr()
        for column in range(n):
            Hx.add(p.H[row,column]*x[column,0])
        model.addConstr(Hx<=p.h[row,0])
#    model.setParam('OutputFlag',False)
    model.optimize()
    if model.Status==2:
            return True # It is not empty
    else:
        # It is empty
        return False

def canonical_polytope(H,h,flag=None,atol=10**-8):
    """
    Given a polytope in form {H x <= h}, provide canonical polytope by finding and removing redundant rows
    Also scale H such that its largest absolute element is 1
    """
    row=0
    while row<H.shape[0]:
        if check_redundancy_row(H,h,row,atol)==True: # The row should be removed
            (H,h)=remove_row(H,h,row)
            row=row
            print("solved a linear program and removed a row")
        else:
            row+=1
            print("solved a linear program but no row removal happened")
    return normalize(H,h)
    
def remove_row(H,h,row):
    """
    Given {x| Hx <= h}, remove the row'th row of H and h
    """
    remaining_rows=list(range(0,row))+list(range(row+1,H.shape[0]))
    H_new=H[remaining_rows,:]
    h_new=h[remaining_rows,:]
    assert H_new.shape[0]==H.shape[0]-1
    return (H_new,h_new)
                    

def check_redundancy_row_old(H,h,ROW,atol=10**-8):
    """
    Solve a linear program to check if the ROW'th row of {x|Hx<=h} is redundant
    """
    model=Model("Row Redundancy Check")
    n=H.shape[1]
    x=np.empty((n,1),dtype='object')
    for row in range(n):
        x[row,0]=model.addVar(lb=-GRB.INFINITY,ub=GRB.INFINITY)
    model.update()
    for row in [r for r in range(H.shape[0]) if r!=ROW]:
        Hx=LinExpr()
        for column in range(n):
            Hx.add(H[row,column]*x[column,0])
        model.addConstr(Hx<=h[row,0])
    J=LinExpr()
    for column in range(n):
        J.add(H[ROW,column]*x[column,0])
    model.setObjective(J, GRB.MAXIMIZE)
    model.setParam('OutputFlag',False)
    model.optimize()
    if model.Status==2:
        if J.getValue()>h[ROW,0]+atol:
            return False # It is NOT redundant
        else:
            return True # It is redudant
    else:
        # Model becomes infeasible (unbounded)
        return False
    
def check_redundancy_row(H,h,ROW,atol=10**-8):
    """
    Solve a linear program to check if the ROW'th row of {x|Hx<=h} is redundant
    """
    model=Model("Row Redundancy Check")
    n=H.shape[1]
    x=np.empty((n,1),dtype='object')
    for row in range(n):
        x[row,0]=model.addVar(lb=-GRB.INFINITY,ub=GRB.INFINITY)
    model.update()
    for row in [r for r in range(H.shape[0]) if r!=ROW]:
        Hx=LinExpr()
        for column in range(n):
            Hx.add(H[row,column]*x[column,0])
        model.addConstr(Hx<=h[row,0])
    J=LinExpr()
    for column in range(n):
        J.add(H[ROW,column]*x[column,0])
    model.setObjective(J, GRB.MAXIMIZE)
    model.setParam('OutputFlag',False)
    model.optimize()
    if model.Status==2:
        if J.getValue()>h[ROW,0]+atol:
            return False # It is NOT redundant
        else:
            return True # It is redudant
    else:
        # Model becomes infeasible (unbounded)
        return False
    
    
def normalize(H,h):
    H_max=np.amax(abs(H),axis=1)
    H=np.dot(np.diag(1/H_max),H)
    h=np.dot(np.diag(1/H_max),h)
    return (H,h)