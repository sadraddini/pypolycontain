import warnings
import numpy as np

# Scipy
try:
    import scipy.linalg as spa
    from scipy.spatial import ConvexHull
    from scipy.linalg import block_diag
except:
    warnings.warn("You don't have scipy package installed. You may get error while using some feautures.")

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
    

# Pypolycontain
try:
    import pypolycontain as pp
except:
    warnings.warn("You don't have pypolycontain properly installed. Can not import objects")


"""
Optimization-based Operations:
"""  


def affine_map( T, P, t=None ):
    """
    Returns the affine map of a polytope.
    """        
    if type(t)==type(None):
            t=np.zeros((P.n,1))
    if P.type=='AH_polytope':
        return pp.AH_polytope(t=t+np.dot(T,P.t),T=np.dot(T,P.T),P=P.P)
    elif P.type=='zonotope':
        return pp.zonotope(x=t+np.dot(T,P.x),G=np.dot(T,P.G))
    elif P.type=="H_polytope":
        if T.shape[0]>=T.shape[1]:
            Tinv=np.linalg.pinv(T)
            H=np.dot(P.H,Tinv)
            return pp.H_polytope(H=H,h=P.h+np.dot(H,t))
        else:
            Q=pp.to_AH_polytope(P)
            return affine_map( T, Q, t )
    else:
        return ValueError('Polytope type: ',P.type," Not recognized")

def translate( t, P ):
    """
    Shifts the polytope by t vector
    """ 
    assert t.shape[0]==P.n # Dimension match       
    if P.type=='AH_polytope':
        return pp.AH_polytope(t=t+P.t,T=P.T,P=P.P)
    elif P.type=='zonotope':
        return pp.zonotope(x=t+P.x,G=P.G)
    elif P.type=="H_polytope":
        return pp.H_polytope(H=P.H,h=P.h+np.dot(P.H,t))
    else:
        return ValueError('Polytope type: ',P.type," Not recognized")
              
def point_membership(Q,x,tol=10**-5,solver="gurobi"):
    if type(Q).__name__=="H_polytope":
        return Q.if_inside(x,tol)
    else:
        Q=pp.to_AH_polytope(Q)
        prog=MP.MathematicalProgram()
        zeta=prog.NewContinuousVariables(Q.P.H.shape[1],1,"zeta")
        prog.AddLinearConstraint(A=Q.P.H,ub=Q.P.h+tol,lb=-np.inf*np.ones((Q.P.h.shape[0],1)),vars=zeta)
        prog.AddLinearEqualityConstraint(Q.T,x-Q.t,zeta)
        if solver=="gurobi":
            result=gurobi_solver.Solve(prog,None,None)
        elif solver=="osqp":
            prog.AddQuadraticCost(np.eye(zeta.shape[0]),np.zeros(zeta.shape),zeta)
            result=OSQP_solver.Solve(prog,None,None)
        else:
            result=MP.Solve(prog)
    return result.is_success()

def check_subset(P1,P2,k=-1):
    """
    Checks if .math..`P1 \subseteq P2` 
    Inputs:
        
    """

def point_membership_fuzzy(Q,x,tol=10**-5,solver="gurobi"):
    """
    Fuzzy membership check. If x contains NaN, the entry is unconstrained
    @param Q: Polytope in R^n
    @param x: n*1 numpy array, may contain NaNs
    @param tol:
    @param solver: solver to use
    @return: boolean of whether x is in Q
    """
    Q=pp.to_AH_polytope(Q)
    prog=MP.MathematicalProgram()
    zeta=prog.NewContinuousVariables(Q.P.H.shape[1],1,"zeta")
    prog.AddLinearConstraint(A=Q.P.H,ub=Q.P.h+tol,lb=-np.inf*np.ones((Q.P.h.shape[0],1)),vars=zeta)
    assert(x.shape[1]==1)
    for i, xi in enumerate(x):
        if not np.isnan(xi):
            prog.AddLinearEqualityConstraint(np.atleast_2d(Q.T[i,:]),(x[i]-Q.t[i]).reshape([-1,1]),zeta)
    if solver=="gurobi":
        result=gurobi_solver.Solve(prog,None,None)
    elif solver=="osqp":
        prog.AddQuadraticCost(np.eye(zeta.shape[0]),np.zeros(zeta.shape),zeta)
        result=OSQP_solver.Solve(prog,None,None)
    else:
        result=MP.Solve(prog)
    return result.is_success()

def check_non_empty(Q,tol=10**-5,solver="gurobi"):
    Q=pp.to_AH_polytope(Q)
    prog=MP.MathematicalProgram()
    zeta=prog.NewContinuousVariables(Q.P.H.shape[1],1,"zeta")
    prog.AddLinearConstraint(A=Q.P.H,ub=Q.P.h+tol,lb=-np.inf*np.ones((Q.P.h.shape[0],1)),vars=zeta)
    if solver=="gurobi":
            result=gurobi_solver.Solve(prog,None,None)
    elif solver=="osqp":
        prog.AddQuadraticCost(np.eye(zeta.shape[0]),np.zeros(zeta.shape),zeta)
        result=OSQP_solver.Solve(prog,None,None)
    else:
        result=MP.Solve(prog)
    return result.is_success()

def directed_Hausdorff_distance(Q1,Q2,ball="infinty_norm",solver="gurobi"):
    r"""
    Computes the directed Hausdorff distance of Q_1 and Q_2 (AH_polytopes)
    ***************************************************************************
    The optimization problem is:
                        Minimize    epsilon  
                        such that   Q1 \subset Q2+epsilon(Ball)
                        
    It is zero if and only if Q1 subset Q2. The method is based on 
                
                    Sadraddini&Tedrake, 2019, CDC (available on ArXiv)
                    
    We solve the following problem:
        D*ball+Q1 subset Q2
    We solve the following linear program:
    ..math::
        \min     D
        s.t.    Lambda_1 H_1=H_2 Gamma_1
                Lambda_2 H_1=H_ball Gamma_2
                Lambda_1 h_1<=h_2 + H_2 beta_1
                Lambda_2 h_2<=D h_ball + H_ball beta_2
                x_2 - X_2 beta_1 - beta_2 = x_1
                X_2 Gamma_1 + Gamma_2 = X_1
    ***************************************************************************
    """
    Q1,Q2=pp.to_AH_polytope(Q1),pp.to_AH_polytope(Q2)
    n=Q1.t.shape[0]
    if ball=="infinty_norm":
        HB=np.vstack((np.eye(n),-np.eye(n)))
        hB=np.vstack((np.ones((n,1)),np.ones((n,1))))
    elif ball=="l1":
        HB,hb=make_ball(ball)
    prog=MP.MathematicalProgram()
    # Variables
    D=prog.NewContinuousVariables(1,1,"D")
    Lambda_1=prog.NewContinuousVariables(Q2.P.H.shape[0],Q1.P.H.shape[0],"Lambda_1")
    Lambda_2=prog.NewContinuousVariables(HB.shape[0],Q1.P.H.shape[0],"Lambda2")
    Gamma_1=prog.NewContinuousVariables(Q2.P.H.shape[1],Q1.P.H.shape[1],"Gamma1")
    Gamma_2=prog.NewContinuousVariables(HB.shape[1],Q1.P.H.shape[1],"Gamma1")
    beta_1=prog.NewContinuousVariables(Q2.P.H.shape[1],1,"beta1")
    beta_2=prog.NewContinuousVariables(HB.shape[1],1,"beta1")
    # Constraints
    # Lambda_1 and Lambda_2 positive
    prog.AddBoundingBoxConstraint(0,np.inf,Lambda_1)
    prog.AddBoundingBoxConstraint(0,np.inf,Lambda_2)
    # Lambda_1 H_1
    Lambda_H_Gamma(prog,Lambda_1,Q1.P.H,Q2.P.H,Gamma_1)
    # Lambda_2 H_1
    Lambda_H_Gamma(prog,Lambda_2,Q1.P.H,HB,Gamma_2)
    # Lambda_1 h_1
    Lambda_h_Inequality(prog,Lambda_1,beta_1,Q2.P.H,Q1.P.h,Q2.P.h)
    # Lambda_2 h_1
    Lambda_h_Inequality_D(prog,Lambda_2,beta_2,HB,Q1.P.h,hB,D)
    # X2 beta_1   
    prog.AddLinearEqualityConstraint(-np.hstack((Q2.T,np.eye(n))),Q1.t-Q2.t,np.vstack((beta_1,beta_2)))
    # X2 Gamma_1
    Aeq=np.hstack((Q2.T,np.eye(Q2.T.shape[0])))
    for i in range(Gamma_1.shape[1]):
        beq=Q1.T[:,i]
        var=np.hstack((Gamma_1[:,i],Gamma_2[:,i]))
        prog.AddLinearEqualityConstraint(Aeq,beq,var)
    # Cost
    # Optimize
    if solver=="gurobi":
            prog.AddLinearCost(D[0,0])
            result=gurobi_solver.Solve(prog,None,None)
    elif solver=="osqp":
        prog.AddQuadraticCost(D[0,0]*D[0,0])
        result=OSQP_solver.Solve(prog,None,None)
    else:
        result=MP.Solve(prog)
    if result.is_success():
        return np.asscalar(result.GetSolution(D))

def Hausdorff_distance(Q1,Q2,ball="infinty_norm",solver="gurobi"):
    return max(directed_Hausdorff_distance(Q1,Q2,ball,solver),directed_Hausdorff_distance(Q2,Q1,ball,solver))
    
def distance_polytopes(Q1,Q2,ball="infinity",solver="gurobi"):
    Q1,Q2=pp.to_AH_polytope(Q1),pp.to_AH_polytope(Q2)
    n=Q1.n
    prog=MP.MathematicalProgram()
    zeta1=prog.NewContinuousVariables(Q1.P.H.shape[1],1,"zeta1")
    zeta2=prog.NewContinuousVariables(Q2.P.H.shape[1],1,"zeta2")
    delta=prog.NewContinuousVariables(n,1,"delta")
    prog.AddLinearConstraint(A=Q1.P.H,ub=Q1.P.h,lb=-np.inf*np.ones((Q1.P.h.shape[0],1)),vars=zeta1)
    prog.AddLinearConstraint(A=Q2.P.H,ub=Q2.P.h,lb=-np.inf*np.ones((Q2.P.h.shape[0],1)),vars=zeta2)
    prog.AddLinearEqualityConstraint( np.hstack((Q1.T,-Q2.T,np.eye(n))),Q2.t-Q1.t,np.vstack((zeta1,zeta2,delta)) )
    if ball=="infinity":
        delta_abs=prog.NewContinuousVariables(1,1,"delta_abs")
        prog.AddBoundingBoxConstraint(0,np.inf,delta_abs)
        prog.AddLinearConstraint(np.greater_equal( np.dot(np.ones((n,1)),delta_abs),delta,dtype='object' ))
        prog.AddLinearConstraint(np.greater_equal( np.dot(np.ones((n,1)),delta_abs),-delta,dtype='object' ))
        cost=delta_abs
    elif ball=="l1":
        delta_abs=prog.NewContinuousVariables(n,1,"delta_abs")
        prog.AddBoundingBoxConstraint(0,np.inf,delta_abs)
        prog.AddLinearConstraint(np.greater_equal( delta_abs,delta,dtype='object' ))
        prog.AddLinearConstraint(np.greater_equal( delta_abs,-delta,dtype='object' ))
        cost=np.dot(np.ones((1,n)),delta_abs)
    else:
        raise NotImplementedError
    if solver=="gurobi":
        prog.AddLinearCost(cost[0,0])
        result=gurobi_solver.Solve(prog,None,None)
    elif solver=="osqp":
        prog.AddQuadraticCost(cost[0,0]*cost[0,0])
        result=OSQP_solver.Solve(prog,None,None)
    else:
        prog.AddLinearCost(cost[0,0])
        result=MP.Solve(prog)
    if result.is_success():
        return np.sum(result.GetSolution(delta_abs)),\
            np.dot(Q1.T,result.GetSolution(zeta1).reshape(zeta1.shape[0],1))+Q1.t,\
            np.dot(Q2.T,result.GetSolution(zeta2).reshape(zeta2.shape[0],1))+Q2.t

def _setup_program_distance_point(P,ball="infinity",solver="Gurobi"):
    """
    Initilize the mathematial program
    Choice of balls:
        infinity: L-infinity norm
        l1: l1 norm (Manhattan Distance)
        l2: l2 norm (Euclidean Distance)
    """
    if P.distance_program is None:
        prog=MP.MathematicalProgram()
        Q=pp.to_AH_polytope(P)
        n=Q.n
        x=np.zeros((n,1))
        P.zeta=prog.NewContinuousVariables(Q.P.H.shape[1],1,"zeta")
        delta=prog.NewContinuousVariables(n,1,"delta")
        prog.AddLinearConstraint(A=Q.P.H,ub=Q.P.h,lb=-np.inf*np.ones((Q.P.h.shape[0],1)),vars=P.zeta)
        P.distance_constraint=prog.AddLinearEqualityConstraint( np.hstack((Q.T,-np.eye(n))),x-Q.t,np.vstack((P.zeta,delta)) )
        if ball=="infinity":
            delta_abs=prog.NewContinuousVariables(1,1,"delta_abs")
            prog.AddBoundingBoxConstraint(0,np.inf,delta_abs)
            prog.AddLinearConstraint(np.greater_equal( np.dot(np.ones((n,1)),delta_abs),delta,dtype='object' ))
            prog.AddLinearConstraint(np.greater_equal( np.dot(np.ones((n,1)),delta_abs),-delta,dtype='object' ))
            prog.AddLinearCost(delta_abs[0,0])
        elif ball=="l1":
            delta_abs=prog.NewContinuousVariables(n,1,"delta_abs")
            prog.AddBoundingBoxConstraint(0,np.inf,delta_abs)
            prog.AddLinearConstraint(np.greater_equal( delta_abs,delta,dtype='object' ))
            prog.AddLinearConstraint(np.greater_equal( delta_abs,-delta,dtype='object' ))
            cost=np.dot(np.ones((1,n)),delta_abs)
            prog.AddLinearCost(cost[0,0])
        elif ball=="l2":
            prog.AddQuadraticCost(np.eye(n),np.zeros(n),delta)
        else:
            print(("Not a valid choice of norm",str(ball)))
            raise NotImplementedError
        P.distance_program=prog
        return 
    else:
        return
            
        
def distance_point_polytope(P, x, ball="infinity", solver="Gurobi"):
    """
    Computes the distance of point x from AH-polytope Q 
    """
    x_vector = np.atleast_2d(x) #in case x is not n*1 vector
    P = pp.to_AH_polytope(P)
    _setup_program_distance_point(P,ball,solver)
    prog=P.distance_program
    Q=pp.to_AH_polytope(P)
    a=P.distance_constraint.evaluator()
    x_vector=x_vector.reshape(max(x_vector.shape),1)
    a.UpdateCoefficients(np.hstack((Q.T,-np.eye(Q.n))), x_vector - Q.t)
    if solver=="Gurobi":
        result=gurobi_solver.Solve(prog,None,None)
    elif solver=="osqp":
        result=OSQP_solver.Solve(prog,None,None)
    else:
        result=MP.Solve(prog)
    if result.is_success():
        zeta_num=result.GetSolution(P.zeta).reshape(P.zeta.shape[0],1)
        x_nearest=np.dot(Q.T,zeta_num)+Q.t
        delta=(x_vector - x_nearest).reshape(Q.n)
        if ball=="infinity":
            d=np.linalg.norm(delta,ord=np.inf)
        elif ball=="l1":
            d=np.linalg.norm(delta,ord=1)
        elif ball=="l2":
            d=np.linalg.norm(delta,ord=2)  
        return d,x_nearest
    
def bounding_box(Q,solver="Gurobi"):
    r"""
    Computes the bounding box of a polytope by solving :math:`2n` linear programs. 
    Each linear program is in the form:
        
    .. math:: 
        \begin{array}{lll}
        l_i= & \min & e_i^T x \\
        & \text{subject to} & x \in \mathbb{P}
        \end{array}
        
    and
    
    .. math:: 
        \begin{array}{lll}
        u_i= & \max & e_i^T x \\
        & \text{subject to} & x \in \mathbb{P}
        \end{array}  
        
    where :math:`l,u` define the lower and upper corners of the bounding box.
    """
    Q=pp.to_AH_polytope(Q)
    prog=MP.MathematicalProgram()
    zeta=prog.NewContinuousVariables(Q.P.H.shape[1],1,"zeta")
    x=prog.NewContinuousVariables(Q.n,1,"x")
    prog.AddLinearConstraint(A=Q.P.H,ub=Q.P.h,lb=-np.inf*np.ones((Q.P.h.shape[0],1)),vars=zeta)
    prog.AddLinearEqualityConstraint(np.hstack((-Q.T,np.eye(Q.n))),Q.t,np.vstack((zeta,x)))
    lower_corner=np.zeros((Q.n,1))
    upper_corner=np.zeros((Q.n,1))
    c=prog.AddLinearCost(np.ones(Q.n),0,x)
    if solver=="Gurobi":
        solver=gurobi_solver
    else:
        raise NotImplementedError
    a=np.zeros((Q.n,1))
    # Lower Corners
    for i in range(Q.n):
        e=c.evaluator()
        a[i,0]=1
        e.UpdateCoefficients(a.reshape(Q.n))
        result=solver.Solve(prog,None,None)
        assert result.is_success()
        lower_corner[i,0]=result.GetSolution(x)[i]
        a[i,0]=0
    # Upper Corners
    for i in range(Q.n):
        e=c.evaluator()
        a[i,0]=-1
        e.UpdateCoefficients(a.reshape(Q.n))
        result=solver.Solve(prog,None,None)
        assert result.is_success()
        upper_corner[i,0]=result.GetSolution(x)[i]
        a[i,0]=0
    return pp.hyperbox(corners=(lower_corner,upper_corner))
        
        
def directed_Hausdorff_hyperbox(b1,b2):
    """
    The directed Hausdorff hyperbox 
    min epsilon such that b1 \in b2+epsilon
    """       
    return max(0,np.max(np.hstack((b1.u-b2.u,b2.l-b1.l))))           
    
def distance_hyperbox(b1,b2):
    """
    The distance between boxes
    """
    return max(0,np.max(np.hstack((b1.l-b2.u,b2.l-b1.u))))      
    

def make_ball(n,norm):
    # TO COMPLETE
    if norm=="l1":
        pass
    elif norm=="infinity":
        pass
    return 

def get_nonzero_cost_vectors(cost):
     cost[cost == 0] = 1e-3

def AH_polytope_vertices(P,N=200,epsilon=0.001,solver="Gurobi"):
    """
    Returns N*2 matrix of vertices
    """
    try:
        P.vertices_2D
        if type(P.vertices_2D) == type(None):
            raise Exception
    except:
        Q=pp.to_AH_polytope(P)
        v=np.empty((N,2))
        prog=MP.MathematicalProgram()
        zeta=prog.NewContinuousVariables(Q.P.H.shape[1],1,"zeta")
        prog.AddLinearConstraint(A=Q.P.H,ub=Q.P.h,lb=-np.inf*np.ones((Q.P.h.shape[0],1)),vars=zeta)
        theta=1
        c=np.array([np.cos(theta),np.sin(theta)]).reshape(2,1)
        c_T=np.dot(c.T,Q.T)
        get_nonzero_cost_vectors(c_T)
        # get_nonzero_cost_vectors(c_T)
        a=prog.AddLinearCost(np.dot(c_T,zeta)[0,0])
        if solver=="Gurobi":
            solver=gurobi_solver
        else:
            raise NotImplementedError
        for i in range(N):
            theta=i*N/2/np.pi+0.01
            c=np.array([np.cos(theta),np.sin(theta)]).reshape(2,1)
            c_T=np.dot(c.T,Q.T)
            e=a.evaluator()
            cost = c_T.reshape(Q.P.H.shape[1])
            get_nonzero_cost_vectors(cost)
            e.UpdateCoefficients(cost)
            result=solver.Solve(prog,None,None)
            assert result.is_success()
            zeta_n=result.GetSolution(zeta).reshape(zeta.shape)
            v[i,:]=(np.dot(Q.T,zeta_n)+Q.t).reshape(2)
        w=np.empty((4*N,2))
        for i in range(N):
            w[4*i,:]=v[i,:]+np.array([epsilon,epsilon])
            w[4*i+1,:]=v[i,:]+np.array([-epsilon,epsilon])
            w[4*i+2,:]=v[i,:]+np.array([-epsilon,-epsilon])
            w[4*i+3,:]=v[i,:]+np.array([epsilon,-epsilon])
        P.vertices_2D=v,w
        return v,w
    else:
        return P.vertices_2D
    
def convex_hull_of_point_and_polytope(x, Q):
    r"""
    Inputs:
        x: numpy n*1 array
        Q: AH-polytope in R^n
    Returns:
        AH-polytope representing convexhull(x,Q)
    
    .. math::
        \text{conv}(x,Q):=\{y | y= \lambda q + (1-\lambda) x, q \in Q\}.
    """
    Q=pp.to_AH_polytope(Q)
    q=Q.P.H.shape[1]
    new_T=np.hstack((Q.T,Q.t-x))
    new_t=x
    new_H_1=np.hstack((Q.P.H,-Q.P.h))
    new_H_2=np.zeros((2,q+1))
    new_H_2[0,q],new_H_2[1,q]=1,-1
    new_H=np.vstack((new_H_1,new_H_2))
    new_h=np.zeros((Q.P.h.shape[0]+2,1))
    new_h[Q.P.h.shape[0],0],new_h[Q.P.h.shape[0]+1,0]=1,0
    new_P=pp.H_polytope(new_H,new_h)
    return pp.AH_polytope(new_T,new_t,new_P)
    

def minkowski_sum(P1,P2):
    r"""
    Inputs: 
        P1, P2: AH_polytopes
    Returns:
        returns the Mkinkowski sum :math:`P_1 \oplus P_2` as an AH-polytope.
        
    **Background**: The Minkowski sum of two sets is defined as:
        
    .. math::
        A \oplus B = \{ a + b \big | a \in A, b \in B\}.
    """
    Q1,Q2=pp.to_AH_polytope(P1),pp.to_AH_polytope(P2)
    T=np.hstack((Q1.T,Q2.T))
    t=Q1.t+Q2.t
    H=spa.block_diag(*[Q1.P.H,Q2.P.H])
    h=np.vstack((Q1.P.h,Q2.P.h))
    new_P=pp.H_polytope(H,h)
    return pp.AH_polytope(t=t,T=T,P=new_P)  
    
    
def intersection_old(P1,P2):
    """
    Inputs: 
        P1, P2: AH_polytopes :math:`\mathbb{P}_1,\mathbb{P}_2`. Converted to AH-polytopes
    Output:
        returns :math:`\mathbb{P}_1 \cap \mathbb{P}_2` as an AH-polytope
    """
    Q1,Q2=pp.to_AH_polytope(P1),pp.to_AH_polytope(P2)
    T=np.hstack((Q1.T,Q2.T*0))
    t=Q1.t
    H_1=spa.block_diag(*[Q1.P.H,Q2.P.H])
    H_2=np.hstack((Q1.T,-Q2.T))
    H=np.vstack((H_1,H_2,-H_2))
    h=np.vstack((Q1.P.h,Q2.P.h,Q2.t-Q1.t,Q1.t-Q2.t))
    new_P=pp.H_polytope(H,h)
    return pp.AH_polytope(T=T,t=t,P=new_P)

def convex_hull(P1,P2):
    """
    Inputs:
        P1, P2: AH_polytopes
    Output:
        returns :math:`\text{ConvexHull}(\mathbb{P}_1,\mathbb{P}_2)` as an AH-polytope
    """
    Q1,Q2=pp.to_AH_polytope(P1),pp.to_AH_polytope(P2)
    T=np.hstack((Q1.T,Q2.T,Q1.t-Q2.t ))
    H_1 = np.hstack((Q1.P.H,np.zeros((Q1.P.H.shape[0],Q2.P.n)), -Q1.P.h ))
    H_2 = np.hstack((np.zeros((Q2.P.H.shape[0],Q1.P.n)), Q2.P.H, Q2.P.h ))
    H_3 = np.zeros(( 2, Q1.P.n + Q2.P.n +1 ))
    H_3[:,-1:]=np.array([1,-1]).reshape(2,1)
    H=np.vstack((H_1,H_2,H_3))
    h=np.vstack((Q1.P.h*0,Q2.P.h,1,0))
    new_P=pp.H_polytope(H,h)
    return pp.AH_polytope(T=T,t=Q2.t,P=new_P)

def intersection(P1,P2):
    """
    Inputs:
        P1, P2: polytopic objects
    Output:
        returns :math:`\mathbb{P}_1 \cap \mathbb{P}_2` as an AH-polytope
        
        If both objects are H-polytopes, return H-polytope
    """
    if P1.type=="H_polytope" and P2.type=="H_polytope":
        H=np.vstack((P1.H,P2.H))
        h=np.vstack((P1.h,P2.h))
        return pp.H_polytope(H,h)
    else:
        X,Y=pp.to_AH_polytope(P1),pp.to_AH_polytope(P2)
        T=np.hstack(( X.T, np.zeros((X.T.shape[0]  ,  Y.T.shape[1] ))  ))
        H_1 = np.hstack((X.P.H,   np.zeros((X.P.H.shape[0],Y.P.H.shape[1])) ))
        Ty_inv=np.linalg.pinv(Y.T)
        H_2 = np.hstack(( np.linalg.multi_dot([Y.P.H,Ty_inv,X.T]),\
                         np.dot(Y.P.H,np.eye(Y.T.shape[1])-np.dot(Ty_inv,Y.T))  ))
        H=np.vstack((H_1,H_2))
        h=np.vstack((X.P.h,Y.P.h-np.linalg.multi_dot([Y.P.H,Ty_inv,X.t-Y.t])))
        new_P=pp.H_polytope(H,h)
        return pp.AH_polytope(T=T,t=X.t,P=new_P)

    
"""
Pydrake Mathematical Program Helper: Matrix based Constraints
"""
def AddMatrixInequalityConstraint_classical(mathematical_program,A,X,B):
    raise NotImplementedError    
    
def Lambda_h_Inequality(mathematical_program,Lambda,beta,H,h_1,h_2):
    """
    Adds Lambda H-1 \le h_2 + H beta to the Mathematical Program
    """
    M1=np.kron(np.eye(h_2.shape[0]),h_1.T)
    M=np.hstack((M1,-H))
    v1=Lambda.reshape((Lambda.shape[0]*Lambda.shape[1],1))
    v=np.vstack((v1,beta))
    mathematical_program.AddLinearConstraint(A=M,ub=h_2,lb=-np.inf*np.ones(h_2.shape),vars=v)
    
def Lambda_h_Inequality_D(mathematical_program,Lambda,beta,H,h_1,h_2,D):
    """
    Adds Lambda H-1 \le h_2 D + H beta to the Mathematical Program
    """
    M1=np.kron(np.eye(h_2.shape[0]),h_1.T)
    M=np.hstack((M1,-H,-h_2))
    v1=Lambda.reshape((Lambda.shape[0]*Lambda.shape[1],1))
    v=np.vstack((v1,beta,D))
    mathematical_program.AddLinearConstraint(A=M,ub=np.zeros(h_2.shape),lb=-np.inf*np.ones(h_2.shape),vars=v)
    
def positive_matrix(mathematical_program,Lambda):
    """
    All elements are non-negative 
    """
#    q=Lambda.shape[0]
    mathematical_program.AddBoundingBoxConstraint(0,np.inf,Lambda)
#    [mathematical_program.AddLinearConstraint(A=np.eye(q),vars=Lambda[:,i],
#                                              ub=np.inf*np.ones((q,1)),lb=np.zeros((q,1)))
#                                                for i in range(Lambda.shape[1])]
#    q=Lambda.shape[0]*Lambda.shape[1]
#    mathematical_program.AddLinearConstraint(A=np.eye(q),vars=Lambda.reshape(q),ub=np.inf*np.ones((q)),lb=np.zeros((q)))
        
def Lambda_H_Gamma(mathematical_program,Lambda,H_1,H_2,Gamma):
    """
    Lambda H_1 = H_2 Gamma
    """
#    v_1=Lambda.reshape((Lambda.shape[0]*Lambda.shape[1],1))
#    for j in range(Gamma.shape[1]):
#        M1=np.kron(np.eye(H_2.shape[0]),H_1[:,j].reshape(1,H_1.shape[0]))
#        M2=-H_2
#        M=np.hstack((M1,M2))
#        v=np.vstack((v_1,Gamma[:,j].reshape(Gamma.shape[0],1)))
#        mathematical_program.AddLinearEqualityConstraint(M,np.zeros((M.shape[0],1)),v)
    for i in range(Lambda.shape[0]):
        for j in range(Gamma.shape[1]):
            M=np.hstack((H_1[:,j],-H_2[i,:]))
            v=np.hstack((Lambda[i,:],Gamma[:,j]))
            mathematical_program.AddLinearEqualityConstraint(M.reshape(1,M.shape[0]),np.zeros(1),v)




def sorting_generator(G,desired_numberofcolumns):
    """
    The goal is deviding the generator into to parts.
    One part that is used for zonotope order reduction methods.
    And the other part which is used to enforce the reduced zonotope to have the desire order.
    """
    dimension = len(G)

    if dimension == desired_numberofcolumns:
        return G , None
    indexlist = np.argsort(np.linalg.norm( G ,axis=0 ,ord = None))
    sortedG = (G)[:,indexlist]
    G_untouched = sortedG[: , - (desired_numberofcolumns - dimension ): ]
    G_reduced = sortedG[: , : -(desired_numberofcolumns - dimension )  ]
    
    return G_reduced , G_untouched


def boxing_order_reduction(zonotope,desired_order=1):
    """
    boxing method for zonotope order reduction
    inputs: input zonotope , order of the output zonotope
    output: zonotope

    Based on Kopetzki, Anna-Kathrin, Bastian Schürmann, and Matthias Althoff.
    "Methods for order reduction of zonotopes."
    2017 IEEE 56th Annual Conference on Decision and Control (CDC). IEEE, 2017.
    """
    assert(type(zonotope) == pp.zonotope),"TypeError: The first argument need to be from \"zonotope\" class in pypolycontain package "
    assert(type(desired_order) == int or type(desired_order) == float) , "TypeError: The second argument need to be a number greater than 1. \n \
                                                                            It is the order of the reduced zonotpe which equal to the number of \
                                                                                columns over the space dimension."
    assert(desired_order >= 1), "desired order of the outcome zonotope needs to be greater or equal to 1"
    
    x = np.array(zonotope.x)
    G = np.array(zonotope.G)
    
    dimension , numberofcolumns =  G.shape
    desired_numberofcolumns = round(desired_order * dimension)
    
    if numberofcolumns <= desired_numberofcolumns:
        return zonotope
    
    elif dimension == desired_numberofcolumns:
        G_box = np.diag(np.sum(abs( G ) ,axis=1 ))
        return pp.zonotope( G_box , x)  
    
    else:
        G_reduced , G_untouched = sorting_generator( G , desired_numberofcolumns )
        G_box = np.concatenate(  ( np.diag(np.sum(abs( G_reduced ) ,axis=1 )) , G_untouched   ), axis=1  )
        return pp.zonotope( G_box , x)  


def pca_order_reduction(zonotope,desired_order=1):

    """
    PCA method for zonotope order reduction
    inputs: input zonotope , order of the output zonotope
    output: zonotope

    Based on Kopetzki, Anna-Kathrin, Bastian Schürmann, and Matthias Althoff.
    "Methods for order reduction of zonotopes."
    2017 IEEE 56th Annual Conference on Decision and Control (CDC). IEEE, 2017.
    """
    assert(type(zonotope) == pp.zonotope),"TypeError: The first argument need to be from \"zonotope\" class in pypolycontain package "
    assert(type(desired_order) == int or type(desired_order) == float) , "TypeError: The second argument need to be a number greater than 1. \n \
                                                                            It is the order of the reduced zonotpe which equal to the number of \
                                                                                columns over the space dimension."
    assert(desired_order >= 1), "desired order of the outcome zonotope needs to be greater or equal to 1"
    
    x = np.array(zonotope.x)
    G = np.array(zonotope.G)
    
    dimension , numberofcolumns =  G.shape
    desired_numberofcolumns = round(desired_order * dimension)

    if numberofcolumns <= desired_numberofcolumns:
        return zonotope

    G_reduced , G_untouched = sorting_generator( G , desired_numberofcolumns )
    X = np.concatenate(  (G_reduced ,  - G_reduced), axis=1 ).T
    covariance = np.dot(X.T , X)
    U,_,_ =np.linalg.svd( covariance )
    interval_hull = boxing_order_reduction( pp.zonotope( np.dot(U.T , G_reduced) , x ) ).G 
    G_pca = np.dot(U , interval_hull)
    
    if type(G_untouched)== np.ndarray:
        return pp.zonotope( np.concatenate((  G_pca , G_untouched   ), axis=1  ) , x) 
    elif G_untouched == None:
        return pp.zonotope(G_pca , x)


def decompose(zonotope,dimensions):
    """
    @author: kasra
    Decompising a given set into bunch of fewer dimensional sets such that 
    the Cartesian product of those sets is a subset of the given set.
    """

    assert(sum(dimensions) == len(zonotope.G)), "ValueError: sum of the given dimensions has \
                                                to be equal to the dimension of the input set"
    #number_of_sets = len(dimensions)

    prog=MP.MathematicalProgram()

    #defining varibales
    x_i = [prog.NewContinuousVariables(i) for i in dimensions]
    G_i = [prog.NewContinuousVariables(i,i) for i in dimensions]                #non-symmetric G_i
    #G_i = [prog.NewSymmetricContinuousVariables(i) for i in dimensions ]               #symmentric G_i

    #inbody_zonotope
    inbody_x = np.concatenate(x_i)
    inbody_G = block_diag(*G_i)
    inbody_zonotope = pp.zonotope(inbody_G , inbody_x)

    #Defining Constraints
    prog.AddPositiveSemidefiniteConstraint(inbody_G)
    pp.subset(prog, inbody_zonotope , zonotope)
    # ASSUMPTION for PSD of inbody_G
    # ASSUMPTION for SYMMENTRICITY of G_i

    #Defining the cost function
    prog.AddMaximizeLogDeterminantSymmetricMatrixCost( inbody_G )

    #Solving
    result = MP.Solve(prog)
    
    print(f"Is optimization successful? {result.is_success()}")
    print(f"Solution to x_i: {result.GetSolution(x_i[0])}")
    print(f"Solution to G_i: {result.GetSolution(G_i[0])}")
    print(f"optimal cost: {result.get_optimal_cost()}")
    print('solver is: ', result.get_solver_id().name())

    x_i_result = [result.GetSolution(x_i[i]) for i in range(len(dimensions))]
    G_i_result = [result.GetSolution(G_i[i]) for i in range(len(dimensions))]
    list_zon = [pp.zonotope(G = G_i_result[i],x = x_i_result[i]) for i in range(len(dimensions))]

    return list_zon