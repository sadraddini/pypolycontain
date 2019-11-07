# Numpy
import numpy as np
# Pypolycontain
#from ..utils.utils import unique_rows


class H_polytope():
    r"""
    An H-polytope is a set defined as follows:

    .. math:: 
        \mathbb{P}=\{x \in \mathbb{R}^n  | H x \le h \},
        
    where 
    
    Attributes:
        * :math:`H \in \mathbb{R}^{q \times n}`:  `numpy.ndarray[float[[q,n]]`
        * :math:`h\in \mathbb{R}^q` : `numpy.ndarray[float[[q,1]]`
    
    define the hyperplanes. The inequality is 
    interpreted element-wise and q is the number of hyperplanes. 
    """
    def __init__(self,H,h,symbolic=False,color='red'):
        # Sanity Checks (H should be q*n matrix and h should be q*1 vector)
        h=np.atleast_2d(h)
        if h.shape[0]==1:
            h=h.T
        if h.shape[1]!=1:
            ValueError("Error: not appropriate h size, it is",h.shape)
        if H.shape[0]!=h.shape[0]:
            ValueError("Error: not consistent dimension of H: %d and h: %d"%(H.shape[0],h.shape[0]))
        # End of sanity checks
        self.H,self.h=H,h
        self.type="H_polytope"
#        if type(h[0,0]):
#            self.H,self.h=unique_rows(H,h)
        self.n=H.shape[1]
        self.hash_value = None
        self.distance_program=None
        self.color=color

    def __repr__(self):
        return ("H_polytope in R^%d"%self.n)

    def __hash__(self):
        if self.hash_value is None:
            self.hash_value = hash(str(np.hstack([self.H, self.h])))
        return self.hash_value
    
    def __add__(self,another_one):
        return minkowski_sum(self,another_one)
 
class zonotope():
    r"""
    A Zonotope is a set defined as follows:

    .. math:: 
        \mathbb{Z}=\langle x,G \rangle = \{x + G p | p \in [-1,1]^q \},
        
    where 
    
    Attributes:
        * :math:`G \in \mathbb{R}^{n \times q}`: `numpy.ndarray[float[[n,q]]` is the zonotope generator. 
        * :math:`x\in \mathbb{R}^n`: `numpy.ndarray[float[[n,1]]` is the zonotope centroid. 
    
    The order
    of the zonotope is defined as :math:`\frac{q}{n}`. 
    """
    def __init__(self,x,G,name=None,color="green"):
        self.x=x
        self.G=G
        if name==None:
            self.name="zonotope"
        else:
            self.name=name
        self.color=color
        self.hash_value = None
        self.distance_program=None
#        self.color="red"

    def __repr__(self):
        return self.name

    def __hash__(self):
        if self.hash_value is None:
            self.hash_value = hash(str(np.hstack([self.G, self.x])))  # FIXME: better hashing implementation
        return self.hash_value

   
class AH_polytope():
    r"""
    An AH_polytope is an affine transformation of an H-polytope and is defined as:
        
    .. math::
        \mathbb{Q}=\{t+Tx  | p \in \mathbb{R}^p, H p \le h \}
    
    Attributes:
        * P: The underlying H-polytope :math:`P:\{x \in \mathbb{R}^p | Hx \le h\}`
        * T: :math:`\mathbb{R}^{n \times p}` matrix: linear transformation
        * t: :math:`\mathbb{R}^{n}` vector: translation 
    """
    def __init__(self,t,T,P,color='blue'):
        """
        Initilization: T,t,P. X=TP+t
        """
        self.T=T # Matrix n*n_p
        self.t=t # vector n*1
        self.P=P # Polytope in n_p dimensions
        self.n=T.shape[0]
        if T.shape[1]!=P.H.shape[1]:
            ValueError("Error: not appropriate T size, it is",T.shape[1],P.n)
        self.method="Gurobi"
        self.hash_value = None
        self.distance_program=None
        self.vertices_2D=None
        self.color=color

    def __repr__(self):
        return "AH_polytope from R^%d to R^%d"%(self.P.n,self.n)

    def __hash__(self):
        if self.hash_value is None:
            self.hash_value = hash(self.P) + hash(str(np.hstack([self.T, self.t])))  # FIXME: better hashing implementation
        return self.hash_value

class V_polytope():
    r"""
    V-polytopes are a convex hull of vertices. 
    .. math:: 
        \mathbb{V}= \{ x \in \mathbb{R}^n | x = \sum_{i=1}^N \lambda_i v_i,  \sum_{i=1}^N \lambda_i=1, \lambda_i \ge 0, i=1,\cdots,N \}
    where each :math:`v_i, i=1,\cdots,N` is a point (some or all are effectively vertices).  
    
    Attributes:
        * list_of_vertices= `list` of `numpy.ndarray[float[n,1]]`.  
    """
    def __init__(self,list_of_vertices):
        self.list_of_vertices=list_of_vertices
        
        
class hyperbox():
    def __init__(self,N=None,corners=[],d=1):
        """
        returns N-dimensional Box 
        corners=typle of 2 numpy arrays: lower_corner and upper_corner
        """
        if len(corners)==0:
            H=np.vstack((np.eye(N),-np.eye(N)))
            h=d*np.ones((2*N,1))
            self.l=-np.ones((N,1))*d
            self.u=-self.l
        else:
            l,u=corners[0:2]
            assert all(u>=l)
            N=l.shape[0]
            assert N==u.shape[0]
            H=np.vstack((np.eye(N),-np.eye(N)))
            if not all(u>=l):
                raise ValueError("Upper-right corner not totally ordering lower-left corner")
            h=np.vstack((u,-l))
            self.l=l
            self.u=u
        self.H_polytope=H_polytope(H,h)
        self.zonotope=zonotope((self.l+self.u)/2,np.diagflat((self.u-self.l)/2))
        self.n=N

class unitbox():
    r"""
    A unitbox in :math:`\mathbb{R}^n` is :math:`[-1,1]^n`.
    """    
    def __init__(self,N):
        H=np.vstack((np.eye(N),-np.eye(N)))
        h=np.ones((2*N,1))
        self.H_polytope=H_polytope(H,h)
        self.zonotope=zonotope(x=np.zeros((N,1)),G=np.eye(N))
        
        
        
        
def box(N=None,d=1,corners=[]):
    """
    returns N-dimensional Box 
    corners=typle of 2 numpy arrays: lower_corner and upper_corner
    """
    if len(corners)==0:
        H=np.vstack((np.eye(N),-np.eye(N)))
        h=d*np.ones((2*N,1))
    else:
        l,u=corners[0:2]
        N=l.shape[0]
        assert N==u.shape[0]
        H=np.vstack((np.eye(N),-np.eye(N)))
        if not all(u>=l):
            raise ValueError("Upper-right corner not totally ordering lower-left corner")
        h=np.vstack((u,-l))
    return H_polytope(H,h)