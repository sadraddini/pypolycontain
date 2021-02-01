#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 31 14:07:32 2021

@author: sadra
"""
from __future__ import print_function
import sys
import threading
from time import sleep
try:
    import thread
except ImportError:
    import _thread as thread
    
    
import warnings
import numpy as np

# Pydrake
try:
    import pydrake.solvers.mathematicalprogram as MP
    import pydrake.solvers.gurobi as Gurobi_drake
    # use Gurobi solver
    global gurobi_solver, license
    gurobi_solver=Gurobi_drake.GurobiSolver()
    license = gurobi_solver.AcquireLicense()
except:
    warnings.warn("You don't have pydrake installed properly. Methods that rely on optimization may fail.")
    
# Pypolycontain
try:
    import pypolycontain as pp
except:
    warnings.warn("You don't have pypolycontain properly installed. Can not import objects")

class affine_system:
    def __init__(self,A=None,B=None,c=None,X=None,U=None,XU=None,name='affine system'):
        self.name=name
        self.A,self.B,self.c,self.X,self.U,self.XU=A,B,c,X,U,XU
        self.build()
    
    def build(self):
        self.n,self.m= self.B.shape
        assert self.A.shape[1]==self.n
        assert self.c.shape==(self.n,1)
        
    def __repr__(self):
        return self.name
            
        
class pwa_system(affine_system):
    def __init__(self,A=None,B=None,c=None,X=None,U=None,XU=None,name='my system'):
        self.name=name
        self.modes={}
    
    def add_mode(self,affine_system):
        self.modes[affine_system.name]=affine_system
        self.build()
        
    def build(self):
        for mode in self.modes:
            S=self.modes[mode]
            self.n,self.m=S.n,S.m
        
    def __repr__(self):
        return self.name
    
    
def point_trajectory(system, start, T, goal, Q=None, R=None):
    n,m=system.n,system.m
    if type(Q)==type(None):
        Q=np.eye(n)
    if type(R)==type(None):
        R=np.eye(m)
    S=list(system.modes)    
    S_all=S+['all']
    prog=MP.MathematicalProgram()
    x={(i,t): prog.NewContinuousVariables(n,1,"x%s%d"%(i,t)) \
       for i in S_all for t in range(T+1)}
    u={(i,t): prog.NewContinuousVariables(m,1,"x%s%d"%(i,t)) \
       for i in S_all for t in range(T)}
    mu={(i,t): prog.NewBinaryVariables(m,1,"x%s%d"%(i,t)) \
       for i in S for t in range(T)}
    for i in S:
        # print( i,-system.modes[i].A, -system.modes[i].B, n, m  )
        XU=system.modes[i].XU
        _N=np.hstack((XU.H,-XU.h))
        for t in range(T):
            _w=np.vstack((x[i,t],u[i,t],mu[i,t]))
            prog.AddLinearConstraint(A=_N,ub=np.zeros((XU.h.shape)),\
                                     lb=-np.inf*np.ones(XU.h.shape),vars=_w)

    for t in range(T):            
        _M=np.hstack([np.hstack((-system.modes[i].A,-system.modes[i].B,\
                                 -system.modes[i].c)) for i in S]\
                     + [np.eye(n)])
        _v=np.vstack([np.vstack((x[i,t],u[i,t],mu[i,t])) for i in S]+[x['all',t+1]])
        prog.AddLinearEqualityConstraint(_M,np.zeros((n,1)),_v) 
    
        
    for t in range(T):
        _u=np.vstack( [u[i,t] for i in S_all])
        _uI=np.hstack( [np.eye(system.m) for i in S]+[-np.eye(m)] ) # Very non-efficient
        prog.AddLinearEqualityConstraint(_uI,np.zeros((m,1)),_u) 
        
        _mu=np.vstack( [mu[i,t] for i in S])
        prog.AddLinearEqualityConstraint(np.ones((1,len(S))),np.ones((1,1)),_mu) 
        
    for t in range(T+1):
        _x=np.vstack( [x[i,t] for i in S_all])        
        _xI=np.hstack( [np.eye(n) for i in S]+[-np.eye(n)] )
        prog.AddLinearEqualityConstraint(_xI,np.zeros((n,1)),_x) 
    
    # start
    prog.AddLinearConstraint( np.equal( x['all',0], start, dtype='object').flatten())
    
    # end
    if type(goal)==type(np.array([])):
        prog.AddLinearConstraint( np.equal( x['all',T], goal, dtype='object').flatten()) 
    else:
    # Goal set: H_polytope
        prog.AddLinearConstraint( np.less_equal( \
            np.dot(goal.H, x['all',T]), goal.h, dtype='object').flatten()) 
    
    
    # Cost function
    for t in range(T+1):
        prog.AddQuadraticCost( Q, np.zeros(n), x['all',t] )
    for t in range(T):
        prog.AddQuadraticCost( R, np.zeros(n), u['all',t] )
    # solve and result
    result=gurobi_solver.Solve(prog,None,None)
    if result.is_success():
        print('trajectory optimization succesfull')
        x_n={t: result.GetSolution(x["all",t]) for t in range(T+1)}
        u_n={t: result.GetSolution(u["all",t]) for t in range(T)}
        mu_n={(t,i):result.GetSolution(mu[i,t]).item() for i in S for t in range(T)}
        sum_mu_0= sum([mu_n[t,S[0]] for t in range(T)])
        if sum_mu_0>0 and sum_mu_0<T:
            print("contact mode change detected")
        return x_n,u_n,mu_n
    else:
        print('trajectory optimization failed')
        return 


   
def polytopic_trajectory(system, start, T, list_of_goals, q=None, Q=None, R=None):
    n,m=system.n,system.m
    if type(Q)==type(None):
        Q=np.eye(n)
    if type(R)==type(None):
        R=np.eye(m)
    if type(q)==type(None):
        q=n
    S=list(system.modes)    
    S_all=S+['all']
    prog=MP.MathematicalProgram()
    x={(i,t): prog.NewContinuousVariables(n,1,"x%s%d"%(i,t)) \
       for i in S_all for t in range(T+1)}
    u={(i,t): prog.NewContinuousVariables(m,1,"x%s%d"%(i,t)) \
       for i in S_all for t in range(T)}
    mu={(i,t): prog.NewBinaryVariables(1,1,"x%s%d"%(i,t)) \
       for i in S for t in range(T)}

    G={(i,t): prog.NewContinuousVariables(n,q,"x%s%d"%(i,t)) \
       for i in S_all for t in range(T+1)}
    theta={(i,t): prog.NewContinuousVariables(m,q,"x%s%d"%(i,t)) \
       for i in S_all for t in range(T)}    
    
    # Containment
    for i in S:
        for t in range(T):
            XU=system.modes[i].XU
            xu=np.vstack((x[i,t],u[i,t]))
            Gtheta=np.vstack((G[i,t],theta[i,t]))
            inbody=pp.zonotope(x=xu,G=Gtheta)
            circumbody=pp.H_polytope(XU.H, XU.h*mu[i,t])
            pp.subset(prog, inbody, circumbody)
            
    # Dynamics of point
    for t in range(T):            
        _M=np.hstack([np.hstack((-system.modes[i].A,-system.modes[i].B,\
                                 -system.modes[i].c)) for i in S]\
                     + [np.eye(n)])
        _v=np.vstack([np.vstack((x[i,t],u[i,t],mu[i,t])) for i in S]+[x['all',t+1]])
        prog.AddLinearEqualityConstraint(_M,np.zeros((n,1)),_v) 
    # Dynamics of polytopes
        _M=np.hstack([np.hstack((-system.modes[i].A,-system.modes[i].B)) for i in S]\
                     + [np.eye(n)])
        _v=np.vstack([np.vstack((G[i,t],theta[i,t])) for i in S]+[G['all',t+1]])
        for j in range(q):
            prog.AddLinearEqualityConstraint(_M,np.zeros((n,1)),_v[:,j])     
    
    # Summation Equation    
    for t in range(T):
        _u=np.vstack( [u[i,t] for i in S_all])
        _uI=np.hstack( [np.eye(system.m) for i in S]+[-np.eye(m)] ) # Very non-efficient
        prog.AddLinearEqualityConstraint(_uI,np.zeros((m,1)),_u) 
        
        _theta=np.vstack( [theta[i,t] for i in S_all])
        for j in range(q):
            prog.AddLinearEqualityConstraint(_uI,np.zeros((m,1)),_theta[:,j]) 
            
        
        _mu=np.vstack( [mu[i,t] for i in S])
        prog.AddLinearEqualityConstraint(np.ones((1,len(S))),np.ones((1,1)),_mu) 
        
    for t in range(T+1):
        _x=np.vstack( [x[i,t] for i in S_all])        
        _xI=np.hstack( [np.eye(n) for i in S]+[-np.eye(n)] )
        prog.AddLinearEqualityConstraint(_xI,np.zeros((n,1)),_x) 
        
        _G=np.vstack( [G[i,t] for i in S_all]) 
        for j in range(q):
            prog.AddLinearEqualityConstraint(_xI,np.zeros((n,1)),_G[:,j]) 
    
    # start
    prog.AddLinearConstraint( np.equal( x['all',0], start, dtype='object').flatten())
    # end
    mu_d,t_d,T_d=pp.add_disjunctive_subsets(prog,pp.zonotope(x=x['all',T],G=G['all',T]),list_of_goals)
    # pp.subset(prog, pp.zonotope(x=x['all',T],G=G['all',T]), goal)
    
    # Cost function
    for t in range(T+1):
        prog.AddQuadraticCost( Q, np.zeros(n), x['all',t] )
    for t in range(T):
        prog.AddQuadraticCost( R, np.zeros(n), u['all',t] )
    
    # Volume Optimization
    prog.AddLinearCost( G['all',0][0,0]*10+G['all',0][1,1]*1 )
    print("*"*10," Set up a mixed-integer optimization problem","*"*10)
    # solve and result
    result=gurobi_solver.Solve(prog,None,None)
    if result.is_success():
        print('polytopic trajectory optimization succesfull')
        x_n={t: result.GetSolution(x["all",t]) for t in range(T+1)}
        u_n={t: result.GetSolution(u["all",t]) for t in range(T)}
        mu_n={(t,i):result.GetSolution(mu[i,t]).item() for i in S for t in range(T)}
        G_n={t: result.GetSolution(G["all",t]) for t in range(T+1)}
        theta_n={t: result.GetSolution(theta["all",t]) for t in range(T)}   
        # Disjunctive Sets
        # print(mu_d,type(mu_d))
        # print({result.GetSolution(i) for i in mu_d})
        # for i in t_d:
        #     print(result.GetSolution(t_d[i]))
        #     print(result.GetSolution(T_d[i]))
        return x_n,u_n,mu_n,G_n,theta_n
    else:
        print('polytopic trajectory optimization failed')
        return 

        
    
def rci(A,B,X,U,W,eta=0.001,q=1):    
    n,m=A.shape[0],B.shape[1]
    W=pp.zonotope(x=np.zeros((n,1)),G=W.G*eta)
    program=MP.MathematicalProgram()
    q+=n
    program = MP.MathematicalProgram()
    phi=program.NewContinuousVariables(n,q,'phi')
    theta=program.NewContinuousVariables(m,q,'theta')
    alpha=program.NewContinuousVariables(1,'alpha')
    beta=program.NewContinuousVariables(2,'beta')
    program.AddBoundingBoxConstraint(0,1,alpha)
    program.AddBoundingBoxConstraint(0,10,beta)
    K=np.hstack(( (np.dot(A,phi) + np.dot(B,theta))[:,n:] , W.G ))
    program.AddLinearConstraint (  np.equal(K, phi, dtype='object').flatten() )
    inbody=pp.zonotope(x=np.zeros((n,1)), G=(np.dot(A,phi)+np.dot(B,theta))[:,0:n])
    _W=pp.to_AH_polytope(W)
    _W.P.h=_W.P.h*alpha
    _X=pp.to_AH_polytope(X)
    _X.P.h=_X.P.h*beta[0]
    _U=pp.to_AH_polytope(U)
    _U.P.h=_U.P.h*beta[1]
    pp.subset(program, inbody,circumbody=_W)
    pp.subset(program, pp.zonotope(x=np.zeros((n,1)),G=phi),circumbody=_X)
    pp.subset(program, pp.zonotope(x=np.zeros((n,1)),G=theta),circumbody=_U)
    program.AddLinearCost(beta[0])
    program.AddLinearConstraint(beta[0]<=1-alpha[0])
    program.AddLinearConstraint(beta[1]<=1-alpha[0])
    program.AddLinearConstraint(beta[0]>=beta[1])
    result=gurobi_solver.Solve(program,None,None)
    if result.is_success():
        print("sucsess")
        print("betas are",result.GetSolution(beta))
        beta_1_n=result.GetSolution(beta)[0]
        beta_2_n=result.GetSolution(beta)[1]
        alpha_n=result.GetSolution(alpha)[0]
        phi_n= result.GetSolution(phi)
        theta_n= result.GetSolution(theta)
        Omega=pp.zonotope(x=np.zeros((2,1)),G=phi_n/beta_1_n,color='red')
        pp.visualize([X,Omega],title='Robust control invariant set $\mathcal{X}$ (red) \n \
        inside safe set $\mathbb{X}$ (green)',figsize=(6,6),a=0.02)
        print("alpha was",alpha_n)
        omega_U=theta_n/beta_2_n
        print(omega_U)
        print('sum_omega=',np.sum(np.abs(omega_U),1)[0])
        return Omega
    else:
        print("failure")
        
        
def in_the_tree(x,list_of_polytopes):
    for P in list_of_polytopes:
        if np.all(np.dot(P.H,x)<=P.h):
            return True
    return False
    
    
def extend(mysystem,start,T,list_of_nodes,H_rep=True,N_H=500,tol_H=1e-5,color=None):
    if type(color)==type(None):
        color=(np.random.random(),np.random.random(),np.random.random())
    else:
        color=color
    x,u,mu,G,theta=polytopic_trajectory(mysystem, start, T, list_of_nodes)
    Z={t: pp.zonotope(x=x[t],G=G[t],color=color) for t in range(T+1)}
    funnel=[None]*T
    H_funnel=[None]*T
    for t in range(T):
        funnel[t]=pp.convex_hull(Z[t], Z[t+1])
        funnel[t].color=color
        if H_rep:
            H_funnel[t]=pp.ray_shooting_hyperplanes(funnel[t],N=N_H,tol=tol_H)
    # fig,ax=plt.subplots()
    # mu[T,'free'],mu[T,'contact']=1,1
    # ax.plot( [x[t][0] for t in range(T+1)] , [x[t][2] for t in range(T+1)] )
    # ax.plot( [x[t][0] for t in range(T+1) if mu[t,'free']==0] , \
    #          [x[t][2] for t in range(T) if mu[t,'free']==0],'o',\
    #          color='red')
    # ax.plot( [x[t][0] for t in range(T+1) if mu[t,'free']==1] , \
    #          [x[t][2] for t in range(T+1) if mu[t,'free']==1],'o',\
    #          color='blue')
    # pp.visualize([Omega]+[funnel[t] for t in range(T)],ax=ax,fig=fig,a=0.01,alpha=0.99,tuple_of_projection_dimensions=(0,2))
    return funnel,H_funnel,x,mu,G   