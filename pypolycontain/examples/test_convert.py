import warnings
import numpy as np
 

# Pypolycontain
try:
    from . import *
    import objects,operations,visualize,conversions
except:
    warnings.warn("You don't have pypolycontain properly installed. Can not execute 'import pypyplycontain'")
    
G=np.random.random((2,20))-0.8
x=np.random.random((2,1))
Z=objects.zonotope(x,G)
Z2=objects.zonotope(x,G*0.5)

B=objects.unitbox(2)
HB=B.H_polytope

visualize.visualize([Z,HB,conversions.to_AH_polytope(Z2)],alpha=0.7)