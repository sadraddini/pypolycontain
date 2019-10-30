import warnings
import numpy as np
# Scipy
try:
    import scipy.linalg as spa
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