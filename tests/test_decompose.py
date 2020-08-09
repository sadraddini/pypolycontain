"""
@author: kasra
Testing decompose and decompose_cost_functiontion functions 
"""
import numpy as np
import pypolycontain as pp

#Testing decompose_cost_functiontion():
n = 10 
a = [np.eye(2)*i  for i in range(1,n)]
assert(-pp.decompose_cost_function(a) ==  sum([np.log(i**4) for i in range(1,n)] )), "decompose_cost_functiontion() has failed"


print("=================================================================================================")

#Testing decompose()
circumbody = pp.zonotope(10* np.ones((4,4)) , np.ones(4))
dimensions = [2,2]
# circumbody = pp.zonotope(np.array([[1,0],[0,1]]) , [0,0])
# dimensions = [1,1]
x,G = pp.decompose(circumbody,dimensions)


print('x=',x)
print('G=',G)



