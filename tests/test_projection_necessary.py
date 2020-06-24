#%%
import pypolycontain as pp
import numpy as np
# Pypolycontain
from pypolycontain.objects import H_polytope,zonotope,AH_polytope,unitbox,hyperbox,V_polytope
from pypolycontain.containment import subset,extreme_rays_for_containment
from pypolycontain.operations import minkowski_sum
from pypolycontain.conversions import to_AH_polytope

import pydrake.solvers.mathematicalprogram as MP
import pydrake.solvers.gurobi as Gurobi_drake
# use Gurobi solver
global gurobi_solver, license
gurobi_solver=Gurobi_drake.GurobiSolver()
license = gurobi_solver.AcquireLicense()


solver=gurobi_solver

def find_lambda(Q,P,N=-1):
    # First solve for Lambda
    program=MP.MathematicalProgram()
    eps=program.NewContinuousVariables(1,"epsilon")
    Ball=hyperbox(N=P.n).H_polytope
    Ball.h=Ball.h*eps
    P_plus_symbolic=minkowski_sum(P,Ball)
    program.AddLinearCost(np.array([1]),np.array([0]),eps)
    subset(program,Q,P_plus_symbolic,N=-1)
    Theta,Lambda,Gamma,Beta=subset(program,P,Q,N=N)
    result=solver.Solve(program,None,None)
    if result.is_success():
        Lambda_n=result.GetSolution(Lambda)
        eps_n=result.GetSolution(eps)
        print("epsilon is",eps_n)
    else:
        raise ValueError("Not feasible")
    return Lambda_n,eps_n

def find_Hx(P,Q,Lambda_n,N=-1):
    program=MP.MathematicalProgram()
    eps=program.NewContinuousVariables(1,"epsilon")
    hx=program.NewContinuousVariables(P.h.shape[0],1,"hx")
    P.h=hx
    Q1=to_AH_polytope(P)
    Q2=to_AH_polytope(Q)
    Hx,Hy,hx,hy,X,Y,xbar,ybar=Q1.P.H,Q2.P.H,Q1.P.h,Q2.P.h,Q1.T,Q2.T,Q1.t,Q2.t
    qx,qy,nx,ny=Hx.shape[0],Hy.shape[0],X.shape[1],Y.shape[1]
    if N<0:
        Theta=np.eye(qy)
    else:
        Theta=extreme_rays_for_containment(Q,N)
    Gamma=program.NewContinuousVariables(ny,nx,'Gamma')
    beta=program.NewContinuousVariables(ny,1,'beta')
    # Constraints
    program.AddLinearConstraint(np.equal(X,np.dot(Y,Gamma),dtype='object').flatten()) #X=YGamma
    program.AddLinearConstraint(np.equal(ybar-xbar,np.dot(Y,beta),dtype='object').flatten()) 
    program.AddLinearConstraint(np.equal(np.dot(Lambda_n,Hx),np.dot(Theta.T,np.dot(Hy,Gamma)),dtype='object').flatten()) 
    program.AddLinearConstraint(np.less_equal(np.dot(Lambda_n,hx),np.dot(Theta.T,hy)+np.dot(Theta.T,np.dot(Hy,beta)),dtype='object').flatten())
    Ball=hyperbox(N=P.n).H_polytope
    Ball.h=Ball.h*eps
    P_plus_symbolic=minkowski_sum(P,Ball)
    program.AddLinearCost(np.array([1]),np.array([0]),eps)
    subset(program,Q,P_plus_symbolic,N=N)
    result=solver.Solve(program,None,None)
    if result.is_success():
        hx_n=result.GetSolution(hx)
        eps_n=result.GetSolution(eps)
        print("epsilon is",eps_n)
    else:
        raise ValueError("Not feasible")
    return hx_n,eps_n
    
#%%     
    
V=pp.objects.V_polytope([20*np.random.random((2,1))-10 for i in range(15)])
Vc=pp.conversions.to_AH_polytope(V)
Vc.color='blue'

Vc=pp.objects.zonotope(x=np.zeros((2,1)),G=4*(np.random.random((2,5))-0.5),color='green')
pp.visualize.visualize([Vc])

P0=pp.objects.hyperbox(2).H_polytope
q=10
P0.H=np.random.random((q,2))-0.5
P0.H=np.vstack((P0.H,-P0.H))
P0.h=np.random.random((q,1))*2
P0.h=np.vstack((P0.h,P0.h))
for i in range(3):
    pp.visualize.visualize([Vc,P0])
    Lambda,eps=find_lambda(Vc,P0,N=0)
    hx,eps=find_Hx(P0,Vc,Lambda,N=0)
    P0.h=hx.reshape(P0.h.shape)
