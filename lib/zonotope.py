# -*- coding: utf-8 -*-
"""
Created on Tue Nov  6 10:24:55 2018

@author: sadra
"""
import numpy as np
from gurobipy import Model,LinExpr,QuadExpr,GRB
from random import randint

from pypolycontain.lib.inclusion_encodings import subset_zonotope_both,constraints_AB_eq_CD,add_Var_matrix
from pypolycontain.utils.utils import valuation


class zonotope():
    """
    Definition of a Zonotope
    """
    def __init__(self,x,G,name=None,color=None):
        self.x=x
        self.G=G
        if name==None:
            self.name="zonotope "+str(randint(0,1000))
        else:
            self.name=name
        try:
            assert color!=None
            self.color=color
        except:
            self.color=(np.random.random(),np.random.random(),np.random.random())
        self.type="zonotope"
    
    def __repr__(self):
        return self.name


def zonotope_directed_distance(z1,z2):
    """
    A linear program for finding the minimal D such that z1 inside z2+ D*infinity_norm_ball
    """
    model=Model("Zonotope Directed")
    G1=np.hstack((z1.G,z2.x-z1.x))
    (n,N1)=G1.shape
    (n,N2)=z2.G.shape
    alpha=np.empty((N2,N1),dtype='object')
    alpha_abs=np.empty(alpha.shape,dtype='object')
    alpha=add_Var_matrix(model,alpha)
    alpha_abs=add_Var_matrix(model,alpha_abs)
    eta=np.empty((n,N1),dtype='object')
    eta_abs=np.empty(eta.shape,dtype='object')
    eta=add_Var_matrix(model,eta)
    eta_abs=add_Var_matrix(model,eta_abs)
    epsilon=model.addVar(lb=0,obj=1)
    model.update()
    absolute_value(model,alpha,alpha_abs)
    absolute_value(model,eta,eta_abs)
    infinity_norm(model,alpha_abs,1)
    infinity_norm(model,eta_abs,epsilon)
    for row in range(n):
        for column in range(N1):
            lin=LinExpr()
            for k in range(N2):
                lin.add(z2.G[row,k]*alpha[k,column])
            model.addConstr(G1[row,column]-eta[row,column]==lin)
    model.setParam('OutputFlag', False)
    model.optimize()
    return epsilon.X    
            
            
def zonotope_distance(z1,z2,eps_min=0,eps_max=10,eps=0.05):
    """
    Using binary search for finding Hausdorff distance between polytopes
    """
    d=(eps_max+eps_min)/2.0
    if (eps_max-eps_min)<eps:
        return d
    else:
        print("searching between",eps_min,"and",eps_max)
        model=Model("Zonotope Distance")
        subset_zonotope_both(model,z1.x,z1.G,z2.x,np.hstack((z2.G,d*np.eye(z2.G.shape[0]))))
        model.setParam("outputflag",False)
        model.optimize()
        if model.Status==3:
            return zonotope_distance(z1,z2,d,eps_max,eps)
        else:
            return zonotope_distance(z1,z2,eps_min,d,eps)
            
def zonotope_distance_cocenter(z1,z2):
    """
    How much z2 should be added by eps to contain z1
    """
    model=Model("zonotope co cenetr distance")
    (n,N1)=z1.G.shape
    (n,N2)=z2.G.shape
    alpha=np.empty((N2,N1),dtype='object')
    alpha_abs=np.empty(alpha.shape,dtype='object')
    alpha=add_Var_matrix(model,alpha)
    alpha_abs=add_Var_matrix(model,alpha_abs)
    eta=np.empty((n,N1),dtype='object')
    eta_abs=np.empty(eta.shape,dtype='object')
    eta=add_Var_matrix(model,eta)
    eta_abs=add_Var_matrix(model,eta_abs)
    epsilon=model.addVar(lb=0,obj=1)
    model.update()
    model.setParam('OutputFlag', False)
    absolute_value(model,alpha,alpha_abs)
    absolute_value(model,eta,eta_abs)
    infinity_norm(model,alpha_abs,1)
    infinity_norm(model,eta_abs,epsilon)
    for row in range(n):
        for column in range(N1):
            lin=LinExpr()
            for k in range(N2):
                lin.add(z2.G[row,k]*alpha[k,column])
            model.addConstr(z1.G[row,column]-eta[row,column]==lin)
    model.optimize()
#    print valuation(alpha)
    return epsilon.X

def zonotope_inside(z,x):
    """
    Arguments:
        z: a zonotope
        x: a test point
    Returns:
        Boolean: if z is inside x
    """
    model=Model("inside")
    n=z.x.shape[0]
    p=np.empty((z.G.shape[1],1),dtype='object')
    for row in range(p.shape[0]):
        p[row,0]=model.addVar(lb=-1,ub=1)
    model.update()
    constraints_AB_eq_CD(model,np.eye(n),x-z.x,z.G,p)
    model.setParam('OutputFlag', 0)
    model.optimize()
    if model.Status==3:
        return False
    else:
        return True  
    
def zonotope_distance_point(z,x):
    model=Model("inside")
    n=z.x.shape[0]
    p=np.empty((z.G.shape[1],1),dtype='object')
    e=np.empty((z.G.shape[1],1),dtype='object')
    for row in range(p.shape[0]):
        p[row,0]=model.addVar(lb=-GRB.INFINITY,ub=GRB.INFINITY)
        e[row,0]=model.addVar(lb=0,ub=GRB.INFINITY,obj=1)
    model.update()
    constraints_AB_eq_CD(model,np.eye(n),x-z.x,z.G,p)
    for row in range(p.shape[0]):
        model.addConstr(p[row,0]+e[row,0]>=-1)
        model.addConstr(p[row,0]-e[row,0]<=1)
    model.setParam('OutputFlag', 0)
    model.optimize()
    d=np.dot(z.G,np.array([e[row,0].X for row in range(p.shape[0])]).reshape(p.shape[0],1))
    return np.linalg.norm(d)
    

"""
Zonotope Order Reduction Methods
"""

def zonotope_order_reduction_outer(Z,G_r,i_max=100,delta=0.05,scale=1.05):
    """
    Zonotope Order Reduction: Inner Approximation
    """
    eps=np.empty(i_max)
    Z_r_list=[Z]*i_max
    Gamma,gamma_diag=zonotope_order_reduction_initial_Gamma0(Z,G_r)
    G_r=np.dot(G_r,np.diag(gamma_diag[:,0]))*scale
    i=0
    while i<i_max:
        print ("***** iteration",i)
        try:
            Z_r_list[i]=zonotope(Z.x,G_r)
            Gamma=zonotope_order_reduction_Gamma(Z,G_r)
            eps[i]=zonotope_distance_cocenter(Z_r_list[i],Z)
            G_r=zonotope_order_reduction_gradient(Z,G_r,Gamma,eps[i])
            i+=1
        except:
            return Z_r_list[i],Z_r_list[:i],eps[:i]
    return Z_r_list[-1],Z_r_list,eps


def zonotope_order_reduction_inner(Z,G_r,i_max=20):
    """
    Zonotope Order Reduction: Inner Approximation
    """
    Z_r_list=[Z]*i_max
    for i in range(15):
        Gamma=zonotope_order_reduction_initial_Gamma0_inner(Z,G_r)
        G_r=zonotope_order_reduction_inner_alternation(Z,G_r,Gamma)
        Z_r_list[i]=zonotope(Z.x,G_r)
    return (Z_r_list[i],Z_r_list)
 
    
def zonotope_order_reduction_gradient(z,G_r,Gamma,eps,delta=0.1):
    N=G_r.shape[1]
    model=Model("zonotope order reduction")
    G=z.G
    (n,n_G)=G.shape
    G_r_plus=np.empty((n,N),dtype='object')
    G_r_plus=add_Var_matrix(model,G_r_plus)
    Gamma_plus=np.empty((N,n_G),dtype='object')
    Gamma_plus=add_Var_matrix(model,Gamma_plus)
    Gamma_plus_abs=np.empty((N,n_G),dtype='object')
    Gamma_plus_abs=add_Var_matrix(model,Gamma_plus_abs)
    alpha=np.empty((n_G,N),dtype='object')
    alpha_abs=np.empty(alpha.shape,dtype='object')
    alpha=add_Var_matrix(model,alpha)
    alpha_abs=add_Var_matrix(model,alpha_abs)
    eta=np.empty((n,N),dtype='object')
    eta_abs=np.empty((n,N),dtype='object')
    eta=add_Var_matrix(model,eta)
    eta_abs=add_Var_matrix(model,eta_abs)
    epsilon=model.addVar(lb=0,ub=eps,obj=1)
    gamma_norm=model.addVar(lb=0,ub=1,obj=0)
    # The changes:
    delta_r=np.empty((n,N),dtype='object')
    delta_r=add_Var_matrix(model,delta_r,delta=delta)
    delta_gamma=np.empty((N,n_G),dtype='object')
    delta_gamma=add_Var_matrix(model,delta_gamma,delta=delta)
    model.update()
    for row in range(n):
        for column in range(N):
            lin=LinExpr()
            for k in range(n_G):
                lin.add(G[row,k]*alpha[k,column])
            model.addConstr(G_r_plus[row,column]-eta[row,column]==lin)
    absolute_value(model,alpha,alpha_abs)
    absolute_value(model,eta,eta_abs)
    absolute_value(model,Gamma_plus,Gamma_plus_abs)
    infinity_norm(model,alpha_abs,1)
    infinity_norm(model,eta_abs,epsilon)
    infinity_norm(model,Gamma_plus_abs,gamma_norm)
    # sums
    sum_matrix_equality(model,G_r_plus,G_r,delta_r)
    sum_matrix_equality(model,Gamma_plus,Gamma,delta_gamma)
#    constraints_AB_eq_CD(model,np.eye(n),G,G_r,Gamma)
    constraints_AB_eq_CD(model,delta_r,Gamma,-G_r,delta_gamma)
    model.setParam('OutputFlag', False)
    model.optimize()
    if model.Status!=2:
        print("model status is", model.Status)
        return 
    print("epsilon",epsilon)
    print("gamma_norm",gamma_norm)
    print(valuation(delta_r),valuation(delta_gamma))
    return valuation(G_r_plus)

def zonotope_order_reduction_initial_Gamma0(z,G_r):
    """
    Takes an initial guess of G_r, produces the best Gamma
    """
    N=G_r.shape[1]
    model=Model("zonotope order reduction")
    G=z.G
    (n,n_G)=G.shape
    Gamma=np.empty((N,n_G),dtype='object')
    Gamma=add_Var_matrix(model,Gamma)
    Gamma_abs=np.empty((N,n_G),dtype='object')
    Gamma_abs=add_Var_matrix(model,Gamma_abs)
    gamma_diag=np.empty((Gamma.shape[0],1),dtype='object')
    for row in range(gamma_diag.shape[0]):
        gamma_diag[row,0]=model.addVar(lb=0)
    model.update()
    absolute_value(model,Gamma,Gamma_abs)
    J=QuadExpr()
    for row in range(Gamma_abs.shape[0]):
        lin=LinExpr()
        for column in range(Gamma_abs.shape[1]):
            lin.add(Gamma_abs[row,column])
        model.addConstr(lin<=gamma_diag[row,0])
        J.add(gamma_diag[row,0]*gamma_diag[row,0])
    constraints_AB_eq_CD(model,np.eye(n),G,G_r,Gamma)
    model.setObjective(J)
    model.setParam('OutputFlag', False)
    model.optimize()
    return (valuation(Gamma),valuation(gamma_diag))

def zonotope_order_reduction_Gamma(z,G_r):
    N=G_r.shape[1]
    model=Model("zonotope order reduction")
    G=z.G
    (n,n_G)=G.shape
    Gamma=np.empty((N,n_G),dtype='object')
    Gamma=add_Var_matrix(model,Gamma)
    Gamma_abs=np.empty((N,n_G),dtype='object')
    Gamma_abs=add_Var_matrix(model,Gamma_abs)
    gamma_diag=np.empty((Gamma.shape[0],1),dtype='object')
    for row in range(gamma_diag.shape[0]):
        gamma_diag[row,0]=model.addVar(lb=0,ub=1)
    model.update()
    absolute_value(model,Gamma,Gamma_abs)
    J=QuadExpr()
    for row in range(Gamma_abs.shape[0]):
        lin=LinExpr()
        for column in range(Gamma_abs.shape[1]):
            lin.add(Gamma_abs[row,column])
        model.addConstr(lin<=gamma_diag[row,0])
        J.add(gamma_diag[row,0]*gamma_diag[row,0])
    constraints_AB_eq_CD(model,np.eye(n),G,G_r,Gamma)
    model.setObjective(J)
    model.setParam('OutputFlag', False)
    model.optimize()
    return valuation(Gamma)
    

def zonotope_order_reduction_initial_Gamma0_inner(z,G_r):
    """
    Takes an initial guess of G_r, produces the best Gamma that G=G_r*Gamma+Delta
    such that ||Delta||<=epsilon, and epsilon is minimzed
    """
    N=G_r.shape[1]
    model=Model("zonotope order reduction")
    G=z.G
    (n,n_G)=G.shape
    Gamma=np.empty((N,n_G),dtype='object')
    Gamma=add_Var_matrix(model,Gamma)
    Gamma_abs=np.empty((N,n_G),dtype='object')
    Gamma_abs=add_Var_matrix(model,Gamma_abs)
    eta=np.empty((n,n_G),dtype='object')
    eta_abs=np.empty(eta.shape,dtype='object')
    eta=add_Var_matrix(model,eta)
    eta_abs=add_Var_matrix(model,eta_abs)
    epsilon=model.addVar(lb=0,obj=1)
    # The changes:
    model.update()
    model.setParam('OutputFlag', False)
    for row in range(n):
        for column in range(n_G):
            lin=LinExpr()
            for k in range(N):
                lin.add(G_r[row,k]*Gamma[k,column])
            model.addConstr(G[row,column]-eta[row,column]==lin)
    absolute_value(model,eta,eta_abs)
    absolute_value(model,Gamma,Gamma_abs)
    infinity_norm(model,eta_abs,epsilon)
    infinity_norm(model,Gamma_abs,1)
    model.optimize()
    if model.Status!=2:
        print("model status is", model.Status)
        return 
    print("epsilon",epsilon)
    return valuation(Gamma)
    
def zonotope_order_reduction_inner_alternation(z,G_r,Gamma,eps=1000):
    N=G_r.shape[1]
    model=Model("zonotope order reduction")
    G=z.G
    (n,n_G)=G.shape
    G_r=np.empty((n,N),dtype='object')
    G_r=add_Var_matrix(model,G_r)
    alpha=np.empty((n_G,N),dtype='object')
    alpha_abs=np.empty(alpha.shape,dtype='object')
    alpha=add_Var_matrix(model,alpha)
    alpha_abs=add_Var_matrix(model,alpha_abs)
    eta=np.empty((n,n_G),dtype='object')
    eta_abs=np.empty(eta.shape,dtype='object')
    eta=add_Var_matrix(model,eta)
    eta_abs=add_Var_matrix(model,eta_abs)
    epsilon=model.addVar(lb=0,ub=eps,obj=1)
    # The changes:
    model.update()
    for row in range(n):
        for column in range(n_G):
            lin=LinExpr()
            for k in range(N):
                lin.add(G_r[row,k]*Gamma[k,column])
            model.addConstr(G[row,column]-eta[row,column]==lin)
    absolute_value(model,alpha,alpha_abs)
    absolute_value(model,eta,eta_abs)
    infinity_norm(model,alpha_abs,1)
    infinity_norm(model,eta_abs,epsilon)
    # sums
    constraints_AB_eq_CD(model,np.eye(n),G_r,G,alpha)
    model.setParam('OutputFlag', False)
    model.optimize()
    if model.Status!=2:
        print("model status is", model.Status)
        return 
    print(("epsilon",epsilon))
    return valuation(G_r)
 
      
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
            
def sum_matrix_equality(model,A,B,C):
    """
    add A=B+C
    """
    for row in range(A.shape[0]):
        for column in range(A.shape[1]):
            model.addConstr(A[row,column]==B[row,column]+C[row,column])