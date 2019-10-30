pypolycontain
=============
pypolycontain is a python package for polytopic objects, operations, and polytope containment Problems


Documentation (NEW - Under Construction):
*****************************************
We are working on the `documentation <https://pypolycontain.readthedocs.io/en/latest/>`_ of this package. Check for updates. 


Publications
***************
    * A preprint of the paper on the theoretical background of this package is available `here <https://arxiv.org/pdf/1903.05214.pdf>`_.
    * A version of this work is going to be presented in `CDC2019 <https://cdc2019.ieeecss.org/>`_


Dependencies:
*************
* numpy


Optional dependencies (for some features)
-----------------------------------------
    * `Drake <https://drake.mit.edu/>`_
    * `pycdd <https://pycddlib.readthedocs.io/en/latest/index.html>`_ 
    * Gurobi 8.0.1 or later `Gurobi Website <https://gurobi.com>`_ *Free Academic License*
    * Scipy




Setup
***************************
As of right now, pip is not available for the latest release. Please clone the repository and add it to your python path.

.. code-block:: bash
    $ git clone https://github.com/sadraddini/pypolycontain.git
    
    $ cd pypolcontain
    
    $ echo "export CLASSPATH=\".:$PWD/:$CLASSPATH\"" >> ~/.bashrc


Polytopic Containment Encodings
###############################

Our package provides linear encodings for the following problems.
    * A-Polytope inside H-Polytope (necessary and sufficient conditions)
    * A-Polytope inside A-Polytope (sufficient conditions)
    * Zonotope inside Zonotope (sufficient conditions)
    * Minkowski Sum of A-Polytopes inside H-Polytope (necessary and sufficient conditions)
    * A-Polytope inside Minkowski Sum of A-Polytopes (sufficient conditions)
    * A-Polytope inside Convex Hull of A-Polytopes (sufficient conditions)
    * Disjunctive A-Polytope inside H-Polytope (necessary and sufficient conditions)

Polytopic Distance Functions:
*****************************
    Available Functions:
    * Hausdorff Distance Between Polytopes
    * Hausdorff Distance Between Zonotopes

Zonotope Order Reduction:
*****************************
    Available Functions:
    * Optimal Zonotope Order Reduction: Outer Approximation
    * Optimal Zonotope Order Reduction: Inner Approximation

Orthogonal Projection:
*****************************
Available Function:
* Orthogonal Projection of Polytopes

.. Examples:
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
    
    We get D12=0, D21=11. The underlying norms are infinity-norm. D12=0 implies the red zonotope is contained within the green zonotope. The Hausdorff distance is max(11,0)=11.
    
    ![](https://imgur.com/jSO5DaM.png "Zonotope red is contained within zonotope Green")
    
    
    
    
    ### Zonotope Containment: Sufficient Conditions (Empirically shown to be very close to necessary)
    ![image|20%](https://imgur.com/bG5ykUa.png "Zonotope Containment")
    ![image | 20%](https://imgur.com/bIHKoUI.png "Zonotope Containment")
    
    ### Zonotope Containment: Sufficient Conditions for Zonotope in the Convex Hull of Zonotopes
    ![](https://imgur.com/GHQo7nf.png "Zonotope Containment")
    
    ### Zonotope Order Reduction: An Optimization-based Approach
    ![](tests/figures/zonotope_reduction_outer.gif "Order Reduction - Outer-Approximation")
    ![](tests/figures/zonotope_reduction_inner.gif "Order Reduction - Inner-Approximation")
    
    ### Orthogonal Projections for Feasible Sets of Model Predictive Controllers
    ![](figures/projection.gif "Orthogonal Projection")