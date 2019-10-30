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

#pycdd    
try:
    from pypolycontain.lib.operations import to_AH_polytope
except:
    warnings.warn("You don't have pyplycontain not properly installed.")
    

"""
***********
Description
***********

The necessary and sufficient conditions for encoding was provided here.
"""


def extreme_rays_for_containment(circumbody,N=0):
    """
    This is from the paper of Sadraddini and Tedrake (2019).
    
    Inputs:
    
        * AH-polytope: ``circumbody``
        * int ``N``: the number of columns added to the subspace of ``U``. 
        For more information, look at the paper.
        *Default* is 0.
    
    Outputs:
    
        * 2D numpy array ``Theta``, which is defined in the paper
    """
    circumbody_AH=to_AH_polytope(circumbody)
    H_y=circumbody_AH.P.H
    Y=circumbody_AH.T
    # First identify K=[H_y'^+  ker(H_y')]
    K_1=np.dot(np.linalg.pinv(H_y.T),Y.T)
    K_2=spa.null_space(H_y.T)
    if K_2.shape[1]>0:
        phi=np.hstack((K_1,K_2))
    else:
        phi=K_1
    # phiw>=0. Now with the augmentation 
    phi_complement=spa.null_space(phi.T)   
    number_of_columns=min(N,phi_complement.shape[1])
    psi=np.hstack((phi,phi_complement[:,0:number_of_columns]))
    p_mat=Matrix(np.hstack((np.zeros((psi.shape[0],1)),psi)))
    p_mat.rep_type = RepType.INEQUALITY
    poly=Polyhedron(p_mat)
    R=np.array(poly.get_generators())
    r=R[:,1:].T
    Theta=np.dot(psi,r)
    return Theta


def subset(program,inbody,circombody,N=-1):
    """
    Adds containment property Q1 subset Q2
    
    Inputs:
        * program: a `pydrake` mathematical program
        * inbody: a polytopic object
        * circumbody: a polytopic object
        * N: 
            * **Default**: :math:`-1``. Sufficient as in Sadraddini and Tedrake (2019a)
            * pick `0` for necessary and sufficient encoding (may be too slow) (2019b)
            * pick any positive number. As the number is smaller, the condition becomes closer to necessity. However, this may be too slow.
    
    Output:
        * No direct output, adds :\math:`inbody \subseteq circumbody` to the model
    """
    Q1=to_AH_polytope(inbody)
    Q2=to_AH_polytope(circombody)
    Hx,Hy,hx,hy,X,Y,xbar,ybar=Q1.P.H,Q2.P.H,Q1.P.h,Q2.P.h,Q1.T,Q2.T,Q1.t,Q2.t
    qx,qy,nx,ny=Hx.shape[0],Hy.shape[0],X.shape[1],Y.shape[1]
    if N<0:
        Theta=np.eye(qy)
    else:
        Theta=extreme_rays_for_containment(circombody,N)
    Lambda=program.NewContinuousVariables(Theta.shape[1],qx,'Lambda')
    Gamma=program.NewContinuousVariables(ny,nx,'Gamma')
    beta=program.NewContinuousVariables(ny,1,'beta')
    # Constraints
    program.AddBoundingBoxConstraint(0,np.inf,Lambda) # Lambda Non-Negative
    program.AddLinearConstraint(np.equal(X,np.dot(Y,Gamma),dtype='object').flatten()) #X=YGamma
    program.AddLinearConstraint(np.equal(ybar-xbar,np.dot(Y,beta),dtype='object').flatten()) 
    program.AddLinearConstraint(np.equal(np.dot(Lambda,Hx),np.dot(Theta.T,np.dot(Hy,Gamma)),dtype='object').flatten()) 
    program.AddLinearConstraint(np.less_equal(np.dot(Lambda,hx),np.dot(Theta.T,hy)+np.dot(Theta.T,np.dot(Hy,beta)),dtype='object').flatten())