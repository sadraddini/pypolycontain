.. pypolycontain documentation master file, created by
   sphinx-quickstart on Mon Oct 28 17:45:43 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to pypolycontain's documentation!
=========================================

.. figure:: ../../pypolycontain.png
    :scale: 120 %


pypolycontain is a python package for polytopic objects, operations, and polytope containment problems. It is written as part of a 
project for verification and control of hybrid systems.

Setup
***************************
As of right now, pip is not available for the latest release. Please clone the repository and add it to your python path.

.. code-block:: python
    $ git clone https://github.com/sadraddini/pypolycontain.git
    $ cd pypolcontain   
    $ echo "export CLASSPATH=\".:$PWD/:$CLASSPATH\"" >> ~/.bashrc
    
    
Dependencies:
*************
* `numpy <https://numpy.org/>`_ 

Optional dependencies (for some features)
-----------------------------------------
    * `Drake <https://drake.mit.edu/>`_
    * `pycdd <https://pycddlib.readthedocs.io/en/latest/index.html>`_ 
    * Gurobi 8.0.1 or later `Gurobi Website <https://gurobi.com>`_ *Free Academic License*
    * Scipy





    
.. toctree::
   :maxdepth: 4
   :caption: Table of Contents:
    
   objects
   operations
   containment
   visualization


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


Developers
==================
* `Sadra Sadraddini <http://www.sadraddini.com/>`_
* `Albert Wu <https://github.com/wualbert/>`_