import warnings
import time
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

try:
    import pydrake.solvers.mathematicalprogram as MP
    import pydrake.solvers.gurobi as Gurobi_drake
    import pydrake.solvers.osqp as OSQP_drake
    # use Gurobi solver
    global gurobi_solver,OSQP_solver, license
    gurobi_solver=Gurobi_drake.GurobiSolver()
    license = gurobi_solver.AcquireLicense()
except:
    warnings.warn("You don't have pydrake installed properly. Methods that rely on optimization may fail.")
        

"""
***********
Description
***********

The necessary and sufficient conditions for encoding was provided here.

"""


def theta_k(circumbody,k=0,i=0):
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


def be_in_set(program,point,zonotope):
    """
    It forces a point to be a member of a zonotope
    """
    dimension = zonotope.G.shape[1]
    b=program.NewContinuousVariables( dimension,'b')
    program.AddBoundingBoxConstraint(-1,1,b)
    program.AddLinearConstraint( np.equal(point, zonotope.x+np.dot(zonotope.G , b) ,dtype='object').flatten() )
    return b

def member_in_set(program,x,P):
    """
    Inputs:
        @ program: mathematical program
        # x: point
        # P: my polytope
    Adds the constraint ..math::`x in P` to the mathematical program
    """
    P=pp.to_AH_polytope(P)
    x.reshape(len(x),1)
    zeta=program.NewContinuousVariables( P.P.n,1, 'zeta')
    program.AddLinearConstraint(A=P.P.H,lb=-np.inf*np.ones((P.P.h.shape[0],1)),ub=P.P.h,vars=zeta)
    _f=np.equal(np.dot(P.T,zeta),x-P.t,dtype='object')
    program.AddLinearConstraint(_f.flatten())
    
    
    
def subset(program,inbody,circumbody,k=-1,Theta=None,i=0,alpha=None,verbose=False):
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
        program.AddLinearConstraint(np.equal(inbody.G,np.dot(circumbody.G,Gamma),dtype='object').flatten())             #inbody_G = circumbody_G * Gamma
        program.AddLinearConstraint(np.equal(circumbody.x - inbody.x ,np.dot(circumbody.G,Lambda),dtype='object').flatten())    #circumbody_x - inbody_x = circumbody_G * Lambda

        Gamma_Lambda = np.concatenate((Gamma,Lambda.reshape(circumbody.G.shape[1],1)),axis=1)
        comb = np.array( list(product([-1, 1], repeat= Gamma_Lambda.shape[1])) ).reshape(-1, Gamma_Lambda.shape[1])
        if alpha=='scalar' or alpha== 'vector':
            comb= np.concatenate( (comb,-1*np.ones((comb.shape[0],1))) , axis=1)
        
        # Managing alpha
        if alpha==None:
            variable = Gamma_Lambda
        elif alpha=='scalar':
            alfa = program.NewContinuousVariables(1,'alpha')
        elif alpha=='vector':
            alfa=program.NewContinuousVariables( circumbody.G.shape[1],'alpha')
            variable = np.concatenate((Gamma_Lambda, alfa.reshape(-1,1)),axis=1)
        else:
            raise ValueError('alpha needs to be \'None\', \'scalaer\', or \'vector\'')

        # infinity norm of matrxi [Gamma,Lambda] <= alfa
        for j in range(Gamma_Lambda.shape[0]):
            program.AddLinearConstraint(
                A= comb,
                lb= -np.inf * np.ones(comb.shape[0]),
                ub= np.ones(comb.shape[0]) if alpha==None else np.zeros(comb.shape[0]),
                vars= variable[j,:] if alpha!='scalar' else np.concatenate((Gamma_Lambda[j,:], alfa ))
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
        if alpha==None:
            return Lambda, Gamma 
        else:
            return Lambda, Gamma , alfa



    Q1=pp.to_AH_polytope(inbody)
    Q2=pp.to_AH_polytope(circumbody)
    Hx,Hy,hx,hy,X,Y,xbar,ybar=Q1.P.H,Q2.P.H,Q1.P.h,Q2.P.h,Q1.T,Q2.T,Q1.t,Q2.t
    qx,qy,nx,ny=Hx.shape[0],Hy.shape[0],X.shape[1],Y.shape[1]
    if type(Theta)!=type(None):
        if verbose:
            print("theta already computed")
        pass
    elif k<0:
        Theta=np.eye(qy)
        if verbose:
            print("Using Positive Orthant")
    else:
        if verbose:
            print("Computing theta with k=%d"%k)
        Theta=theta_k(circumbody,k,i)
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

def necessity_gap(X,Y,Theta):
    """
    The necessity gap for the encoding using cone defined by Theta
    
    .. warning:
        To be added later
    """
    alpha_0=alpha(X,Y,theta_k(Y,k=0))
    my_alpha= alpha(X,Y,Theta)
    return 1-my_alpha/alpha_0

def necessity_gap_k(X,Y,plot=True):
    print("\n\n","="*75,"\n","="*75,"\n\t\tComputing Necessity Gaps ")
    print("="*75,"\n","="*75)
    print("k\tTheta.shape \t delta(X,Y,C) \t alpha \t Computation Time")
    print("-"*75)
    alpha_0=alpha(X,Y,theta_k(Y,k=0))
    Table={}        
    Z=pp.to_AH_polytope(Y)
    kappa=int(Z.T.shape[1]-Z.T.shape[0])
    for k in range(kappa+1):
        t_0=time.time()
        Theta=theta_k(Y,k)
#        print(k,"theta shape",Theta.shape,Theta)
        a=alpha(X,Y,Theta)
        delta=1-a/alpha_0
        cpu_time=time.time()-t_0
        print(k, "\t",Theta.shape,"\t %0.03f"%delta,"\t\t %0.03f"%a,"\t\t %0.003f"%cpu_time)
        Table[k]=np.array([Theta.shape[1],delta,a,cpu_time])
    if plot:
        import matplotlib.pyplot as plt
        plt.figure()
        fig, ax1 = plt.subplots()
        color = (0.7,0,0)
        ax1.set_xlabel(r'$k$')
        ax1.set_ylabel(r'#Rows in $\Theta_k$', color=color,FontSize=20)
        ax1.plot(range(kappa+1),[Table[k][0] for k in range(kappa+1)],'-',color=color)
        ax1.plot(range(kappa+1),[Table[k][0] for k in range(kappa+1)],'o',color=color)
        ax1.tick_params(axis='y', labelcolor=color)
        
        ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
        
        color = (0,0,0.7)
        ax2.set_xlabel(r'$k$')
        ax2.set_ylabel(r'$\delta(\mathbb{X},\mathbb{Y},\mathbb{C}_k)$', color=color,FontSize=20)
        ax2.plot(range(kappa+1),[Table[k][1] for k in range(kappa+1)],'-',color=color)
        ax2.plot(range(kappa+1),[Table[k][1] for k in range(kappa+1)],'o',color=color)
        ax2.tick_params(axis='y', labelcolor=color)
        
        ax2.set_title(r'$k^*=%d, n=%d$'%(kappa,Z.n),FontSize=20)
        ax1.grid()
        fig.tight_layout()  # otherwise the right y-label is slightly clipped
    
    return Table

    
def alpha(X,Y,Theta):
    X2=pp.to_AH_polytope(X)
    Y2=pp.to_AH_polytope(Y)
    prog=MP.MathematicalProgram()
    alpha=prog.NewContinuousVariables(1,"alpha")
    t=prog.NewContinuousVariables(X2.n,1,"t")
    Y2.P.h=Y2.P.h*alpha
    X2.t=X2.t+t
    subset(prog,X2,Y2,Theta=Theta)
    prog.AddLinearCost(np.eye(1),np.zeros((1)),alpha)
    result=gurobi_solver.Solve(prog,None,None)
    if result.is_success():
#        print("alpha test successful")
#        print(result.GetSolution(alpha),result.GetSolution(t))
        return 1/result.GetSolution(alpha)[0]
    else:
        print("not a subset") 