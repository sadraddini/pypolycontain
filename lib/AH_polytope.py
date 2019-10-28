#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 20 15:16:48 2019

@author: sadra
"""

import numpy as np

from scipy.optimize import linprog as lp
from scipy.linalg import block_diag as blk
#from gurobipy import Model,LinExpr,QuadExpr,GRB

from pypolycontain.lib.polytope import polytope,Box
from pypolycontain.utils.utils import valuation
from pypolycontain.lib.operations import distance_point_polytope

class AH_polytope():
    """
    Affine Transformation of an H-polytope
    Attributes:
        P: The underlying H-polytope P:{x in R^q | Hx \le h}
        T: R^(n*q) matrix: linear transformation
        t: R^{n*1) vector: translation
    """
    def __init__(self,T,t,P):
        """
        Initilization: T,t,P. X=TP+t
        """
        self.T=T # Matrix n*n_p
        self.t=t # vector n*1
        self.P=P # Polytope in n_p dimensions
        self.n=T.shape[0]
        if T.shape[1]!=P.H.shape[1]:
            ValueError("Error: not appropriate T size, it is",T.shape[1],P.n)
        self.type="AH_polytope"
        self.method="Gurobi"
        self.hash_value = 0
        self.distance_program = None
        self.color = 'r'
    def __repr__(self):
        return "AH_polytope from R^%d to R^%d"%(self.P.n,self.n)

#    def __hash__(self):
#        if self.hash_value is None:
#            self.hash_value = hash(self.P) + hash(str(np.hstack([self.T, self.t])))  # FIXME: better hashing implementation
#        return self.hash_value

    def is_inside(self,x):
        """
        Return if x is inside AH_polytope by checking the feasibility of a linear program
        """
        if self.method=="Gurobi":
            model=Model("check_if_inside")
            p=tupledict_to_array(model.addVars(list(range(self.P.n)),[0],lb=-GRB.INFINITY,ub=GRB.INFINITY,name="p"))
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
            raise NotImplementedError # fix the -infty bounds first! 
            sol=lp(np.ones((self.P.n,1)), A_ub=self.P.H, b_ub=self.P.h, A_eq=self.T, b_eq=x-self.t)
    #        print sol.message
            return sol.status==0
        else:
            raise ValueError("Method %s not recognized"%self.method)
            
    def is_nonempty(self,tol=10**-6):
        model=Model("check_if_inside")
        p=tupledict_to_array(model.addVars(list(range(self.P.n)),[0],lb=-GRB.INFINITY,ub=GRB.INFINITY,name="p"))
        x=tupledict_to_array(model.addVars(list(range(self.T.shape[0])),[0],lb=-GRB.INFINITY,ub=GRB.INFINITY,name="x"))
        model.update()
        constraints_list_of_tuples(model,[(self.T,p),(np.eye(self.n),self.t),(-np.eye(self.n),x)],sign="=")
        constraints_list_of_tuples(model,[(self.P.H,p),(-np.eye(self.P.h.shape[0]),(self.P.h+tol))],sign="<")
        model.setParam('OutputFlag', False)
        model.optimize()
        if model.Status==3:
#            print("AH-polytope is empty")
            return False
        elif model.Status==2:
#            print("AH-polytope is not empty and a point is",np.array([x[i,0].X for i in range(x.shape[0])]))
            return True
        else:
            return "Model status is %d"%model.Status
        
            
def to_AH_polytope(P):
    if P.type=="AH_polytope":
        return P
    elif P.type=="H_polytope" or P.type=="H-polytope":
        n=P.H.shape[1]
        return AH_polytope(np.eye(n),np.zeros((n,1)),P)
    elif P.type=="zonotope":
        q=P.G.shape[1]
        return AH_polytope(P.G,P.x,Box(q))
    else:
        raise ValueError("P type not understood:",P.type)
        
        
def distance_point(P,x,norm="L2"):
    """
    Distance from a point x to the closet point in AH_polytope 
    Arguments:
        poly: polytope, zonotope, or AH_polytope
        x= numpy array
        norm: choice of L1, L2, or infinity norm
    Returns:
        Updated (June 10)
        float --> distance
        array --> closest point in poly
    """
    poly=to_AH_polytope(P)
    n=poly.T.shape[0]
    x=x.reshape(n,1)
    model=Model("distance to point")
    p=tupledict_to_array(model.addVars(list(range(poly.T.shape[1])),[0],lb=-GRB.INFINITY,ub=GRB.INFINITY,name="p")) # The point in H-polytope
    delta=tupledict_to_array(model.addVars(list(range(n)),[0],lb=-GRB.INFINITY,ub=GRB.INFINITY,name="delta")) # The distance
    model.update()
    constraints_list_of_tuples(model,[(np.eye(n),x),(-np.eye(n),delta),(-np.eye(n),poly.t),(-poly.T,p)],sign="=")
    constraints_list_of_tuples(model,[(poly.P.H,p),(-np.eye(poly.P.h.shape[0]),poly.P.h)],sign="<")
    model.setParam('OutputFlag', False)
    if norm in ["infinity","Linfinity","Inf","Linf"]:
        delta_max=model.addVar(lb=0,obj=1)
        model.update()
        for i in range(n):
            model.addConstr(delta_max>=delta[i,0])
            model.addConstr(delta_max>=-delta[i,0])
        model.optimize()
        raise NotImplementedError
        return delta_max.X #TODO: implement projection point
    elif norm in [2,"L2","Euclidean"]:
        J=QuadExpr()
        for i in range(n):
            J.add(delta[i,0]*delta[i,0])
        model.setObjective(J)
        model.optimize()
        #convert from Gurobi variable to float
        p_float = np.zeros(p.shape)
        for i, p_i in enumerate(np.ndarray.flatten(p)):
            p_float[i] = p_i.X
        # print(np.dot(poly.T,p)+poly.t)
        return model.getObjective().getValue()**0.5,np.dot(poly.T,p_float)+poly.t
    elif norm in [1,"L1","l1"]:
        delta_max=tupledict_to_array(model.addVars(list(range(n)),[0],lb=0,ub=GRB.INFINITY,name="delta_max",obj=1))
        for i in range(n):
            model.addConstr(delta_max[i,0]>=delta[i,0])
            model.addConstr(delta_max[i,0]>=-delta[i,0])
        model.optimize()
        p_num=np.array([p[i,0] for i in range(p.shape[0])]).reshape(p.shape)
        return sum([delta_max[i,0].X for i in range(n)]),np.dot(poly.T,p)+poly.t
    else:
        raise ValueError("The following norm: %s is not identifed"%str(norm))      
    
def is_inside(poly,x,tol=10**-6):
    """
    Boolean answer to whether x is in poly
    """
    if poly.type=="H-polytope":
        return all(np.dot(poly.H,x)<=poly.h)
    else:
        Q=to_AH_polytope(poly)
        return distance_point_polytope(Q,x)<=tol
        
def minimum_distance(poly_1,poly_2,norm="infinity"):
    """
    find the closets points in poly_1 and poly_2
    """
    poly_1,poly_2=to_AH_polytope(poly_1),to_AH_polytope(poly_2)
    n=poly_1.T.shape[0]
    assert n==poly_2.T.shape[0]
    model=Model("minimum_distance")
    delta=tupledict_to_array(model.addVars(list(range(n)),[0],lb=-GRB.INFINITY,ub=GRB.INFINITY,name="delta"))
    x=tupledict_to_array(model.addVars(list(range(n)),[0],lb=-GRB.INFINITY,ub=GRB.INFINITY,name="x"))
    p_1=tupledict_to_array(model.addVars(list(range(poly_1.P.H.shape[1])),[0],lb=-GRB.INFINITY,ub=GRB.INFINITY,name="p_1"))
    p_2=tupledict_to_array(model.addVars(list(range(poly_2.P.H.shape[1])),[0],lb=-GRB.INFINITY,ub=GRB.INFINITY,name="p_2"))
    model.update()
    constraints_list_of_tuples(model,[(np.eye(n),x),(-np.eye(n),delta),(-np.eye(n),poly_1.t),(-poly_1.T,p_1)],sign="=")
    constraints_list_of_tuples(model,[(np.eye(n),x),(-np.eye(n),poly_2.t),(-poly_2.T,p_2)],sign="=")
    constraints_list_of_tuples(model,[(poly_1.P.H,p_1),(-np.eye(poly_1.P.h.shape[0]),poly_1.P.h)],sign="<")
    constraints_list_of_tuples(model,[(poly_2.P.H,p_2),(-np.eye(poly_2.P.h.shape[0]),poly_2.P.h)],sign="<")
    model.setParam('OutputFlag', False)
    if norm=="infinity":
        delta_max=model.addVar(lb=0,obj=1)
        model.update()
        for i in range(n):
            model.addConstr(delta_max>=delta[i,0])
            model.addConstr(delta_max>=-delta[i,0])
        model.optimize()
#        print valuation(x).T,valuation(delta).T
        return delta_max.X
    else:
        raise NotImplementedError
        
def check_collision(poly_1,poly_2,tol=10**-6):
    return minimum_distance(poly_1,poly_2,norm="infinity")<tol
    

def Minkowski_sum(poly_1,poly_2):
    poly_1,poly_2=to_AH_polytope(poly_1),to_AH_polytope(poly_2)
    if poly_1.T.shape[0]!=poly_2.T.shape[0]:
        ValueError("The sizes of two sums do not match",poly_1.T.shape[0],poly_2.T.shape[0])    
    H=blk(poly_1.P.H,poly_2.P.H)
    h=np.vstack((poly_1.P.h,poly_2.P.h))
    P=polytope(H,h)
    T=np.hstack((poly_1.T,poly_2.T))
    t=(poly_1.t+poly_2.t)
    return AH_polytope(T,t,P)

def cartesian_product(poly_1,poly_2):
    poly_1,poly_2=to_AH_polytope(poly_1),to_AH_polytope(poly_2)
    if poly_1.T.shape[0]!=poly_2.T.shape[0]:
        ValueError("The sizes of two sums do not match",poly_1.T.shape[0],poly_2.T.shape[0])    
    H=blk(poly_1.P.H,poly_2.P.H)
    h=np.vstack((poly_1.P.h,poly_2.P.h))
    P=polytope(H,h)
    T=blk(poly_1.T,poly_2.T)
    t=np.vstack((poly_1.t,poly_2.t))
    return AH_polytope(T,t,P)  

def is_nonempty(poly):
    """
    Return if poly (AH-polytope, Zonotope, or Polytope) is non-empty
    """
    Q=to_AH_polytope(poly)
    return Q.is_nonempty()

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
