import warnings
import numpy as np

# Scipy
try:
    import scipy.linalg as spa
except:
    warnings.warn("You don't have scipy package installed. You may get error while using some feautures.")

    
#pycdd    
try:
    from cdd import Polyhedron,Matrix,RepType
except:
    warnings.warn("You don't have CDD package installed. Unable to run cone ray generation.")

# Internal Imports
try:
    import pypolycontain as pp
except:
    warnings.warn("You don't have pypolycontain not properly installed.")
#try:
#    from pypolycontain.conversions import to_AH_polytope
#except:
#    pass
#    warnings.warn("You don't have pypolycontain not properly installed.")

# try:
#    import pydrake.solvers.mathematicalprogram as MP
#    import pydrake.solvers.gurobi as Gurobi_drake
#    import pydrake.solvers.osqp as OSQP_drake
#    # use Gurobi solver
#    global gurobi_solver,OSQP_solver, license
#    gurobi_solver=Gurobi_drake.GurobiSolver()
#    license = gurobi_solver.AcquireLicense()
#    OSQP_solver=OSQP_drake.OsqpSolver()
# except:
#    warnings.warn("You don't have pydrake installed properly. Methods that rely on optimization may fail.")
        

"""
***********
Description
***********

The necessary and sufficient conditions for encoding was provided here.

"""


def extreme_rays_for_containment(circumbody,k=0,i=0):
    """
    This is from the Section 4 of the paper [Sadraddini and Tedrake, 2020].
    
    Inputs:
    
        * AH-polytope: ``circumbody``
        * int ``N``: the number of columns added to the subspace of ``U``. 
        For more information, look at the paper.
        *Default* is 0.
    
    Outputs:
    
        * 2D numpy array ``Theta``, which is defined in the paper
    """
    circumbody_AH=pp.to_AH_polytope(circumbody)
    H_y=circumbody_AH.P.H
    q_y=H_y.shape[0]
    if k<0:
        return np.eye(q_y)
    Y=circumbody_AH.T
    # First identify K=[H_y'^Y'+  ker(H_y')]
    S_1=np.dot(np.linalg.pinv(H_y.T),Y.T)
    S_2=spa.null_space(H_y.T)
    if S_2.shape[1]>0:
        S=np.hstack((S_1,S_2))
    else:
        S=S_1
    # phiw>=0. Now with the augmentation 
    S_complement=spa.null_space(S.T)   
    circumbody.dim_complement=S_complement.shape[1]
    number_of_columns=min(k,S_complement.shape[1])
    psi=np.hstack((S,S_complement[:,i:number_of_columns+i]))
    p_mat=Matrix(np.hstack((np.zeros((psi.shape[0],1)),psi)))
    p_mat.rep_type = RepType.INEQUALITY
    poly=Polyhedron(p_mat)
    R=np.array(poly.get_generators())
    r=R[:,1:].T
    Theta=np.dot(psi,r)
    return Theta


def subset(program,inbody,circumbody,k=-1,Theta=None,i=0):
    """
    Adds containment property Q1 subset Q2
    
    Inputs:
        * program: a `pydrake` mathematical program
        * inbody: a polytopic object
        * circumbody: a polytopic object
        * N: 
            * **Default**: :math:`-1``. Sufficient as in [Sadraddini and Tedrake, 2020]
            * pick `0` for necessary and sufficient encoding (may be too slow) (2019b)
            * pick any positive number. As the number is smaller, 
            the condition becomes closer to necessity. However, this may be too slow.
    
    Output:
        * No direct output, adds :\math:`inbody \subseteq circumbody` to the model
    """
    if type(inbody).__name__=="zonotope" and type(circumbody).__name__=="zonotope":
        """
        For the case when both inbody and circumbody sets are zonotope:
        """
        from itertools import product
        #Defining Variables
        Gamma=program.NewContinuousVariables( circumbody.G.shape[1], inbody.G.shape[1], 'Gamma')
        Lambda=program.NewContinuousVariables( circumbody.G.shape[1],'Lambda')
        
        #Defining Constraints
        program.AddLinearConstraint(np.equal(inbody.G,np.dot(circumbody.G,Gamma),dtype='object').flatten()) #inbody_G = circumbody_G * Gamma
        program.AddLinearConstraint(np.equal(circumbody.x - inbody.x ,np.dot(circumbody.G,Lambda),dtype='object').flatten())

        Gamma_Lambda = np.concatenate((Gamma,Lambda.reshape(circumbody.G.shape[1],1)),axis=1)
        
        comb = np.array(list(product([-1, 1], repeat= Gamma_Lambda.shape[1]))).reshape(-1, Gamma_Lambda.shape[1])
        
        for i in range(Gamma_Lambda.shape[0]):
            program.AddLinearConstraint(
                A=comb,
                lb= -np.inf * np.ones(comb.shape[0]),
                ub=np.ones(comb.shape[0]),
                vars=Gamma_Lambda[i,:]
            ) 
            
        #from pydrake.symbolic import abs as exp_abs
        # Gamma_abs = np.array([exp_abs(Gamma[i,j]) for i in range(circumbody.G.shape[1]) for j in range(inbody.G.shape[1])]).reshape(*Gamma.shape)
        # Lambda_abs = np.array([exp_abs(Lambda[i]) for i in range(circumbody.G.shape[1])]).reshape(circumbody.G.shape[1],1)
        # Gamma_lambda_abs = np.concatenate((Gamma_abs,Lambda_abs),axis=1)
        # infnorm_Gamma_lambda_abs = np.sum(Gamma_lambda_abs,axis=1)

        #[program.AddConstraint( infnorm_Gamma_lambda_abs[i] <= 1) for i in range(circumbody.G.shape[1])]
        
        #program.AddBoundingBoxConstraint(-np.inf * np.ones(circumbody.G.shape[1]) , np.ones(circumbody.G.shape[1]), infnorm_Gamma_lambda_abs)
        # program.AddLinearConstraint(
        #     A=np.eye(circumbody.G.shape[1]),
        #     lb= -np.inf * np.ones(circumbody.G.shape[1]),
        #     ub=np.ones(circumbody.G.shape[1]),
        #     vars=infnorm_Gamma_lambda_abs
        # )

        return Lambda, Gamma



    Q1=pp.to_AH_polytope(inbody)
    Q2=pp.to_AH_polytope(circumbody)
    Hx,Hy,hx,hy,X,Y,xbar,ybar=Q1.P.H,Q2.P.H,Q1.P.h,Q2.P.h,Q1.T,Q2.T,Q1.t,Q2.t
    qx,qy,nx,ny=Hx.shape[0],Hy.shape[0],X.shape[1],Y.shape[1]
    if k<0:
        Theta=np.eye(qy)
    else:
        Theta=extreme_rays_for_containment(circumbody,k,i)
    Lambda=program.NewContinuousVariables(Theta.shape[1],qx,'Lambda')
    Gamma=program.NewContinuousVariables(ny,nx,'Gamma')
    gamma=program.NewContinuousVariables(ny,1,'gamma')
    # Constraints
    program.AddBoundingBoxConstraint(0,np.inf,Lambda) # Lambda Non-Negative
    program.AddLinearConstraint(np.equal(X,np.dot(Y,Gamma),dtype='object').flatten()) #X=YGamma
    program.AddLinearConstraint(np.equal(ybar-xbar,np.dot(Y,gamma),dtype='object').flatten()) 
    program.AddLinearConstraint(np.equal(np.dot(Lambda,Hx),np.dot(Theta.T,np.dot(Hy,Gamma)),dtype='object').flatten()) 
    program.AddLinearConstraint(np.less_equal(np.dot(Lambda,hx),\
          np.dot(Theta.T,hy)+np.dot(Theta.T,np.dot(Hy,gamma)),dtype='object').flatten())
    return Theta,Lambda,Gamma,gamma

    