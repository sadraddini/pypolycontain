"""
@author: kasra
Testing Boxing and PCA order reduction methods by illustration in 2D dimension
"""
import numpy as np
import matplotlib.pyplot as plt
import pypolycontain as pp

numberofexamples = 10
for i in range(numberofexamples):
    x=np.random.rand( 2,1  ) # offset
    numofcolumns = np.random.randint(2,30)
    G=np.random.rand( 2, numofcolumns )
    my_zonotope=pp.zonotope(x=x,G=G)
    my_zonotope.color = 'red'

    order = np.random.uniform( 1, numofcolumns/2 ) if numofcolumns != 2 else 1
    zon_boxing = pp.boxing_order_reduction(my_zonotope , order )
    zon_boxing.color = 'yellow'
    zon_pca = pp.pca_order_reduction(my_zonotope , order)
    zon_pca.color = 'green'

    pp.visualize(
        [zon_boxing , zon_pca , my_zonotope] , 
        title= 'red : original zonotope \n number of columns={0} \n \
            yellow: Boxing \n green: PCA \n order ={1}'.format(numofcolumns,order)
        )
plt.show()






