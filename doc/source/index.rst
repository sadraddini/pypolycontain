.. pypolycontain documentation master file, created by
   sphinx-quickstart on Mon Oct 28 17:45:43 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to pypolycontain's documentation!
=========================================

.. figure:: ../../pypolycontain.png
   :height: 400px
   :width: 1000 px
   :scale: 50 %
   :alt: alternate text
   :align: center


pypolycontain is a python package for polytopic objects, operations, and polytope containment problems. It is written as part of a 
project for verification and control of hybrid systems.

Setup
***************************
Installation is now easy::

    ls -lsa .
    make file    
    
Dependencies:
*************
* `numpy <https://numpy.org/>`_ (Use latest version)

Optional dependencies (for some features)
-----------------------------------------
    * `Drake <https://drake.mit.edu/>`_ (Use latest version)
    * `pycdd <https://pycddlib.readthedocs.io/en/latest/index.html>`_ (Use latest version)
    * `Gurobi <https://gurobi.com>`_  (Version 8.0.1 or later) *Free Academic License*
    * `scipy <https://scipy.org//>`_ (Use latest version)


    
.. toctree::
   :maxdepth: 3
   :caption: Contents:
   :glob:
    
   examples
   objects
   conversions
   operations
   visualization
   


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


Main Developer
==================
* `Sadra Sadraddini <http://www.sadraddini.com/>`_

Contributors
==================
* `Albert Wu <https://github.com/wualbert/>`_
* `Kasra Ghasemi <https://github.com/Kasraghasemi>`_