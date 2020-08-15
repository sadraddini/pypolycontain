"""
@author: kasra
Testing decompose and decompose_cost_functiontion functions 
"""
import numpy as np
import pypolycontain as pp
import matplotlib.pyplot as plt

circumbody = pp.zonotope( 4 * np.random.rand(2,5) , np.ones(2) , color='green')
dimensions = [1,1]
x,G = pp.decompose(circumbody,dimensions)
#x,G = pp.decompose(circumbody,dimensions, obj_coef= [1,1])

print('x=',x)
print('G=',G)
inbody = pp.zonotope([[G[0],0],[0,G[1]]],x , color= 'red')
pp.visualize([circumbody,inbody])
plt.show()


print("=================================================================================================")

#Testing decompose()
circumbody_G = 10 * np.random.rand(6,12)
#circumbody_G = np.random.uniform(-10,10,(6,12))
#circumbody_G = 10 * np.eye(6)
circumbody_x = np.ones(6)
circumbody = pp.zonotope(circumbody_G,circumbody_x)
circumbody_1 = pp.zonotope( circumbody_G[:2,:] , circumbody_x[:2] ,color='green')
circumbody_2 = pp.zonotope( circumbody_G[2:4,:] , circumbody_x[2:4] ,color='green')
circumbody_3 = pp.zonotope( circumbody_G[4:6,:] , circumbody_x[4:6] ,color='green')

dimensions = [2,2,2]
x,G = pp.decompose(circumbody,dimensions)
#x,G = pp.decompose(circumbody,dimensions, obj_coef= [1,1,1])
#G = 1.1* np.array(G)
print('x=',x)
print('G=',G)

inbody_1 = pp.zonotope(G[0],x[0], color= 'red')
inbody_2 = pp.zonotope(G[1],x[1], color= 'red')
inbody_3 = pp.zonotope(G[2],x[2], color= 'red')

inbody = inbody_1 ** inbody_2 **inbody_3

print('directed_Hausdorff_distance = ',pp.directed_Hausdorff_distance(circumbody, inbody))
print('directed_Hausdorff_distance = ',pp.directed_Hausdorff_distance(inbody,circumbody))

fig, axs = plt.subplots(3)
pp.visualize([circumbody_1, inbody_1], ax = axs[0],fig=fig, title='')
axs[0].axis('equal')
pp.visualize([circumbody_2, inbody_2], ax = axs[1],fig=fig, title='')
axs[1].axis('equal')
pp.visualize([circumbody_3, inbody_3], ax = axs[2],fig=fig, title='')
axs[2].axis('equal')
plt.show()
