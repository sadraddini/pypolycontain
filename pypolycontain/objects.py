# Numpy
import numpy as np
from itertools import combinations
from scipy.linalg import block_diag

# Pypolycontain
import pypolycontain as pp

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
#        h=np.atleast_2d(h) # I want the vector to be n*1
        if h.shape[0]==1 or len(h.shape)==1:
            h=h.reshape(len(h),1)
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

    def __add__(self,other):
        return pp.minkowski_sum(self,other) # It is returning the minkowski sum in the form of a AH_Polytope!

    def __rmul__(self,scalar):
        """
        Scaling zonotopes by a scalar. The scalar needs to be an integer or a float.
        """
        return pp.H_polytope(self.H,self.h*scalar)
 
class zonotope():
    r"""
    A Zonotope is a set defined as follows:

    .. math:: 
        \mathbb{Z}=\langle x,G \rangle = \{x + G p | p \in [-1,1]^q \},
        
    where 
    
    Attributes:
        * :math:`G \in \mathbb{R}^{n \times q}`: `numpy.ndarray[float[[n,q]]` is the zonotope generator. 
        * :math:`x\in \mathbb{R}^n`: `numpy.ndarray[float[[n,1]]` is the zonotope centroid. Default is zero vector.
    
    The order
    of the zonotope is defined as :math:`\frac{q}{n}`. 
    """
    def __init__(self,G,x=None,name=None,color="green"):
        self.G=np.array(G)        
        if type(x)==type(None):
            self.x=np.zeros((G.shape[0],1))
        else:
            self.x=np.array(x)
        if name==None:
            self.name="zonotope"
        else:
            self.name=name
        self.color=color
        self.type='zonotope'
        self.hash_value = None
        self.distance_program=None
#        self.color="red"

    def __repr__(self):
        return self.name

    def __hash__(self):
        if self.hash_value is None:
            self.hash_value = hash(str(np.hstack([self.G, self.x])))  # FIXME(Kasraghasemi) : better hashing implementation
        return self.hash_value
    
    def __add__(self,other):
        x = self.x + other.x
        G = np.concatenate([self.G,other.G],axis=1)
        return zonotope(G=G,x=x)

    def __pow__(self,other):
        x = np.concatenate([self.x , other.x])
        G = block_diag(self.G,other.G)
        return zonotope(G=G,x=x)

    def __rmul__(self,other):
        """
        Scaling zonotopes by a scalar. The scalar needs to be an integer or a float.
        """
        assert( type(other)==int or type(other)==float or other.shape == (1,) ), "A zonotope can be multiplied by an int or float"
        x = other * self.x
        G = other * self.G
        return zonotope(G=G,x=x)

    def __mul__(self,other):
        return self.__rmul__(other)

    def volume(self):
        r"""
        Computes the volume of the zonotope in :math:`\mathbb{n}` dimensions. 
        The formula is based on the paper in [Gover2002]_
            
            
            
        .. [Gover2002] 
            Gover, Eugene, and Nishan Krikorian. 
            "Determinants and the volumes of parallelotopes and zonotopes." 
            Linear Algebra and its Applications 433, no. 1 (2010): 28-40.
                
        """
        G=self.G
        S=combinations(range(G.shape[1]),G.shape[0])
        V=0
        for s in S:
            Gs=np.hstack([G[:,i:i+1] for i in s])
            V+=abs(np.linalg.det(Gs))
        return (2**G.shape[0]) *V
    
    def volume_gradient(self):
        G=self.G
        S=combinations(range(G.shape[1]),G.shape[0])
        V_dot=np.zeros(G.shape)
        for s in S:
            Gs=np.hstack([G[:,i:i+1] for i in s])
            D=np.linalg.det(Gs)
            if D!=0:
                adj_Gs=D*np.linalg.inv(Gs)
                X=adj_Gs.T*np.sign(np.linalg.det(Gs))
            else:
#                print("Warning: we have a singular matrix")
                e=np.eye(G.shape[0])*10**(-5)
                adj_Gs=np.linalg.det(Gs+e)*np.linalg.inv(Gs+e)
                X=adj_Gs.T*np.sign(np.linalg.det(Gs+e))
#                X=adj_Gs.T*(np.linalg.det(Gs+e))
            for i in range(len(s)):
                V_dot[:,s[i]]+=X[:,i]
        return V_dot

   
class AH_polytope():
    r"""
    An AH_polytope is an affine transformation of an H-polytope and is defined as:
        
    .. math::
        \mathbb{Q}=\{t+Tx  | x \in \mathbb{R}^p, H x \le h \}
    
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
        self.t=np.atleast_2d(t) # vector n*1
        self.P=P # Polytope in n_p dimensions
        self.type='AH_polytope'
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
    
    def copy(self):
        return AH_polytope(self.t,self.T,H_polytope(self.P.H,self.P.h))

    def __add__(self,other):
        from pypolycontain.operations import minkowski_sum
        return minkowski_sum(self,other)
    
    def __rmul__(self,scalar):
        return AH_polytope(self.t,self.T*scalar,self.P)


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
        self.H_polytope=H_polytope(H,h,color='cyan')
        self.zonotope=zonotope(x=(self.l+self.u)/2,G=np.diagflat((self.u-self.l)/2),color='cyan')
        self.n=N
        self.color='cyan'



class unitbox():
    r"""
    A unitbox in :math:`\mathbb{R}^n` is :math:`[-1,1]^n`.
    """    
    def __init__(self,N):
        H=np.vstack((np.eye(N),-np.eye(N)))
        h=np.ones((2*N,1))
        self.H_polytope=H_polytope(H,h)
        self.zonotope=zonotope(x=np.zeros((N,1)),G=np.eye(N))
        
             
def box(N=None,d=1,corners=[],color='cyan'):
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
    return H_polytope(H,h,color=color)