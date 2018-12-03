# pypolycontain: A Python Package for Polytope Containment Problems

The polytope containment problem is deciding whether a polytope is a contained within another polytope.
This problem has strong relevance in computational convexity, and arises in applications such as verification and control of linear/hybrid systems.
The complexity of the problem heavily depends on how the polytopes are represented. 
We focus on representations given as affine transformations of polytopes that are described by their hyperplanes (H-polytopes). 
Affine transformations of H-polytopes are able to capture a broad range of polytopes for which hyperplane representations are not efficient,
 such as Zonotopes, convex hulls and Minkowski sums of multiple H-polytopes.
 We provide sufficient conditions, or necessary and sufficient conditions for some cases,
 to efficiently cast the polytope containment problem as a set of linear constraints.
 These efficient encodings enable us to designate certain components of  polytopes as decision variables,
 which may be added to a given set of (mixed-integer) linear constraints. 
We present examples on applications to Zonotope algebra and formal controller verification for a hybrid system.

## Dependencies:
* Gurobi 7.0 or later [Gurobi](https://gurobi.com) (Free Academic License)
* (For visualization purposes) [pycdd](https://pycddlib.readthedocs.io/en/latest/index.html)



## Polytopic Containment Encodings
Our package provides linear encodings for the following problems.
* A-Polytope inside H-Polytope (necessary and sufficient conditions)
* A-Polytope inside A-Polytope (sufficient conditions)
* Zonotope inside Zonotope (sufficient conditions)
* Minkowski Sum of A-Polytopes inside H-Polytope (necessary and sufficient conditions)
* A-Polytope inside Minkowski Sum of A-Polytopes (sufficient conditions)
* A-Polytope inside Convex Hull of A-Polytopes (sufficient conditions)
* Disjunctive A-Polytope inside H-Polytope (necessary and sufficient conditions)

## Polytopic Distance Functions:
Available Functions:
* Hausdorff Distance Between Polytopes
* Hausdorff Distance Between Zonotopes

## Zonotope Order Reduction:
Available Functions:
* Optimal Zonotope Order Reduction: Outer Approximation
* Optimal Zonotope Order Reduction: Inner Approximation

## Examples:

### Finding the distance between two zonotopes
```python
import numpy as np
from pypolycontain.lib.zonotope import zonotope,zonotope_directed_distance
from pypolycontain.visualization.visualize_2D import visualize_2D_zonotopes as visZ

G_1=np.array([[1,0,0,1,1],[0,1,0,-1,-3]])
G_2=np.array([[1,0,1,1,2,4,-4],[0,1,1,-1,3,3,1]])
x_1=np.array([0,1]).reshape(2,1)
x_2=np.array([1,0]).reshape(2,1)

z1=zonotope(x_1,G_1,color="red")
z2=zonotope(x_2,G_2,color="green")
visZ([z2,z1],title="Zonotopes")
D12=zonotope_directed_distance(z1,z2)
D21=zonotope_directed_distance(z2,z1) 
```

We get D12=0, D21=11. The underlying norms are infinity-norm. D12=0 implies the red zonotope is contained within green zonotope. 
The Hausdorff distance is max(11,0)=11.
![](https://imgur.com/jSO5DaM.png "Zonotope red is contained within zonotope Green") {:height="50%" width="50%"}




### Zonotope Containment: Sufficient Conditions (Empirically shown to be very close to necessary)
![image|20%](https://imgur.com/bG5ykUa.png "Zonotope Containment")
![image | 20%](https://imgur.com/bIHKoUI.png "Zonotope Containment")

### Zonotope Containment: Sufficient Conditions for Zonotope in Convex Hull of Zonotopes
![](https://imgur.com/GHQo7nf.ong "Zonotope Containment")


### Zonotope Order Reduction: An Optimization-based Approach
![](tests/figures/zonotope_reduction_outer.gif "Order Reduction - Outer-Approximation")
![](tests/figures/zonotope_reduction_inner.gif "Order Reduction - Inner-Approximation")
