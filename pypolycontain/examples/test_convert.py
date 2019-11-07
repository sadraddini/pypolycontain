import warnings
import numpy as np
 

# Pypolycontain
try:
    import pypolycontain
    from pypolycontain.conversions import to_AH_polytope
    from pypolycontain.objects import zonotope,unitbox
    from pypolycontain.visualize import visualize
except:
    warnings.warn("You don't have pypolycontain properly installed. Can not execute 'import pypyplycontain'")
    
G=np.random.random((2,20))-0.8
x=np.random.random((2,1))
Z=zonotope(x,G)
Z2=zonotope(x,G*0.5)

B=unitbox(2)
HB=B.H_polytope

visualize([Z,HB,to_AH_polytope(Z2)],alpha=0.7)