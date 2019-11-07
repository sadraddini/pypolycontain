import warnings
import numpy as np
 
# Scipy
try:
    from scipy.spatial import ConvexHull
except:
    warnings.warn("You don't have scipy package installed. You may get error while using some feautures.")
    
# Pypolycontain
try:
    from objects import H_polytope,zonotope,AH_polytope,unitbox,hyperbox,V_polytope
except:
    warnings.warn("You don't have pypolycontain properly installed. Can not execute 'import pypyplycontain'")

# pycdd
try:
    from cdd import Polyhedron,Matrix,RepType
except:
    print("WARNING: You don't have CDD package installed. Unable to visualize polytopes. You may still visualize zonotopes.")

# Pydrake
try:
    import pydrake.solvers.mathematicalprogram as MP
    import pydrake.solvers.gurobi as Gurobi_drake
    import pydrake.solvers.osqp as OSQP_drake
    # use Gurobi solver
    global gurobi_solver,OSQP_solver, license
    gurobi_solver=Gurobi_drake.GurobiSolver()
    license = gurobi_solver.AcquireLicense()
    OSQP_solver=OSQP_drake.OsqpSolver()
except:
    warnings.warn("You don't have pydrake installed properly. Methods that rely on optimization may fail.")

try:
    import pydrake.solvers.scs as scs_drake
    scs_solver=scs_drake.ScsSolver()
except:
    warnings.warn("You don't have pydrake with SCS solver.")
    

def to_AH_polytope(P):
    """
    Converts the polytopic object P into an AH-polytope
    """
    if type(P).__name__=="AH_polytope":
        return P
    elif type(P).__name__=="H_polytope":
        n=P.H.shape[1]
        return AH_polytope(np.eye(n),np.zeros((n,1)),P)
    elif type(P).__name__=="zonotope":
        q=P.G.shape[1]
        return AH_polytope(T=P.G,t=P.x,P=unitbox(N=q).H_polytope,color=P.color)
    elif type(P).__name__=="V_polytope":
        V=P.list_of_vertices
        N=len(V)
        T=np.hstack([V[i]-V[-1] for i in range(N-1)])
        t=V[-1]
        H=np.vstack((-np.eye(N-1),np.ones((1,N-1))))
        h=np.zeros((N,1))
        h[N-1,0]=1
        P=H_polytope(H,h)
        return AH_polytope(t,T,P)
    else:
        raise ValueError("object type not within my polytopic library:",P.type)

def H_to_V(P):
    r"""
    Returns the vertices of an H_polytope.
    
    Inputs:
        * P: H_polytope in :math:`\mathbb{R}^n`
    Output:
        * V: list of vertices. Each vertex is *numpy.ndarray[float[n,1]]*
    
    **Method**:
        The method is based on double description method, and is using pycddlib.
    .. warning::
        This method can be very slow or numerically unstable for polytopes in high dimensions and/or large number of hyperplanes
    """
    if type(P).__name__=="H_polytope":
        p_mat=Matrix(np.hstack((P.h,-P.H)))
        p_mat.rep_type = RepType.INEQUALITY
        poly=Polyhedron(p_mat)
        y=np.array(poly.get_generators())
        x=y[:,1:]#
        x=x[ConvexHull(x).vertices,:]
        return x
    else:
        raise ValueError(str(type(P).__name__)+" is not an H-polytope")
        

def AH_to_V(P,N=360,epsilon=1e-3,solver="Gurobi"):
    """
    Returns the V-polytope form of a 2D AH_polytope.
    The method is based on ray shooting. 
    
    Inputs:
        * P: AH-polytope
        * N ``defualt=360``: number of rays
        * solver: ``default=Gurobi``. The linear-programming optimization solver.
    Returns:
        * V: matrix 
    """
    Q=to_AH_polytope(P)
    if Q.n!=2:
        raise ValueError("Sorry, but I can only do 2D")
    v=np.empty((N,2))
    prog=MP.MathematicalProgram()
    zeta=prog.NewContinuousVariables(Q.P.H.shape[1],1,"zeta")
    prog.AddLinearConstraint(A=Q.P.H,ub=Q.P.h,lb=-np.inf*np.ones((Q.P.h.shape[0],1)),vars=zeta)
    theta=1
    c=np.array([np.cos(theta),np.sin(theta)]).reshape(2,1)
    c_T=np.dot(c.T,Q.T)
    c_T[c_T == 0] = 1e-3
    # get_nonzero_cost_vectors(c_T)
    a=prog.AddLinearCost(c_T.reshape(Q.P.n),np.zeros(1),zeta)
    if solver=="Gurobi":
        solver=gurobi_solver
    elif solver=="SCS":
        solver=scs_solver
    else:
        raise NotImplementedError
    for i in range(N):
        theta=i*N/2/np.pi+0.01
        c=np.array([np.cos(theta),np.sin(theta)]).reshape(2,1)
        c_T=np.dot(c.T,Q.T)
        e=a.evaluator()
        cost = c_T.reshape(Q.P.H.shape[1])
        cost[cost == 0] = 1e-3
        e.UpdateCoefficients(cost)
        result=solver.Solve(prog,None,None)
        assert result.is_success()
        zeta_n=result.GetSolution(zeta).reshape(zeta.shape)
        v[i,:]=(np.dot(Q.T,zeta_n)+Q.t).reshape(2)
    try:
        v=v[ConvexHull(v).vertices,:]
        return v
    except: # convexhull very small. Add some epsilon
        w=np.empty((4*N,2))
        for i in range(N):
            w[4*i,:]=v[i,:]+np.array([epsilon,epsilon])
            w[4*i+1,:]=v[i,:]+np.array([-epsilon,epsilon])
            w[4*i+2,:]=v[i,:]+np.array([-epsilon,-epsilon])
            w[4*i+3,:]=v[i,:]+np.array([epsilon,-epsilon])
        w=w[ConvexHull(w).vertices,:]
        return w
    
def zonotope_to_V(Z):
    """
    Finds the vertices of a zonotope
    """
    q=Z.G.shape[1]
    if q<13:
        v=Z.x.T+np.dot(Z.G,vcube(q).T).T
        return v[ConvexHull(v).vertices,:]
    else:
        warnings.warn('The number of generators %d is very large. Resorting to ray shooting'%q)
        return AH_to_V(to_AH_polytope(Z))
    
def to_V(P):
    r"""
    returns the vertices of polytopic object in a vertical stack form
    """
    if type(P).__name__=='zonotope':
        return zonotope_to_V(P)
    elif type(P).__name__=='AH_polytope':
        return AH_to_V(P)
    elif type(P).__name__=='H_polytope':
        return H_to_V(P)
    else:
        raise ValueError("Did not recognize the polytopic object"+str(type(P).__name__))
    

        
        
        


     
#V=[np.random.random((2,1)) for i in range(10)]
#P=V_polytope(V)
#Q=to_AH_polytope(P)
#
#B=unitbox(4)
#HB=B.H_polytope
#V=H_to_V(HB)
#Q=to_AH_polytope(V_polytope(V))
#ver=AH_to_V(Q)
        
#G=np.random.random((2,5))
#x=np.random.random((2,1))
#Z=zonotope(x,G)
#V=zonotope_to_V(Z)

def vcube(n):
    r"""
    :math:`2^n \times n` array of vectors of vertices in unit cube in :math:`\mathbb{R^n}`
    """
    from itertools import product 
    v=list(product(*list(zip([-1]*n,[1]*n))))
    return np.array(v)