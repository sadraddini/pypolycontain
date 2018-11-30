# pypolycontain: A Python Package for Polytope Containment Problems

The polytope containment problem is deciding whether a polytope is a contained within another polytope. This problem has strong relevance in computational convexity, and arises in applications such as verification and control of linear/hybrid systems. The complexity of the problem heavily depends on how the polytopes are represented. We focus on representations given as affine transformations of polytopes that are described by their hyperplanes (H-polytopes). Affine transformations of H-polytopes are able to capture a broad range of polytopes for which hyperplane representations are not efficient, such as Zonotopes, convex hulls and Minkowski sums of multiple H-polytopes.  We provide sufficient conditions, or necessary and sufficient conditions for some cases, to efficiently cast the polytope containment problem as a set of linear constraints. These efficient encodings enable us to designate certain components of  polytopes as decision variables, which may be added to a given set of (mixed-integer) linear constraints. We present examples on applications to Zonotope algebra and formal controller verification for a hybrid system.

## Dependencies:
* Gurobi 7.0 or later
* (For visualization purposes) [pycdd](https://pycddlib.readthedocs.io/en/latest/index.html)



## Description: Under Construction
Our package provides linear encodings for the following problems.
* A-Polytope inside H-Polytope (necessary and sufficient conditions)
* A-Polytope inside A-Polytope (sufficient conditions)
* Zonotope inside Zonotope (sufficient conditions)
* Minkowski Sum of A-Polytopes inside H-Polytope (necessary and sufficient conditions)
* A-Polytope inside Minkowski Sum of A-Polytopes (sufficient conditions)
* A-Polytope inside Convex Hull of A-Polytopes (sufficient conditions)
* Disjunctive A-Polytope inside H-Polytope (necessary and sufficient conditions)



### Zonotope Order Reduction: An Optimization Approach
![](tests/figures/zonotope_reduction_outer.gif "Order Reduction") [tests/figures/zonotope_reduction_outer.gif]
