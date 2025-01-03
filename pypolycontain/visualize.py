import numpy as np
import pypolycontain as pp
#import warnings

from matplotlib.patches import Polygon
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection

# Pypolycontain
#try:
#    from pypolycontain.conversions import to_V    
#except:
#    pass
#    warnings.warn("You don't have pypolycontain properly installed. Can not import conversions")
    

def visualize(list_of_objects,fig=None,ax=None,a=0.5,alpha=0.8,tuple_of_projection_dimensions=[0,1],\
              title=r'pypolycontain visualization',\
              show_vertices=False,TitleSize=15,fontsize=15,equal_axis=False,grid=True,\
              N_points=1000,figsize=(8,8)):
    r"""
    Visualization.
    
    inputs: 
        * list_of_objects:
        * fig:
        * tuple_of_projection_dimensions: 
        * 
    """
    if type(ax)==type(None):
        fig,ax=plt.subplots()
    fig.set_size_inches(figsize[0],figsize[1])
    p_list,x_all=[],np.empty((0,2))  
    for p in list_of_objects:
        if p.n>2:
            print('projection on ',tuple_of_projection_dimensions[0],\
                  ' and ',tuple_of_projection_dimensions[1], 'dimensions')
            p=_projection(p,tuple_of_projection_dimensions)
        x=pp.to_V(p,N=N_points)
        mypolygon=Polygon(x)
        p_list.append(mypolygon) 
        x_all=np.vstack((x_all,x))
        if show_vertices:
            ax.plot(x[:,0],x[:,1],'*',color=p.color)
    p_patch = PatchCollection(p_list,color=[p.color for p in list_of_objects], alpha=alpha)
    ax.add_collection(p_patch)
    ax.set_xlim([np.min(x_all[:,0])-a,a+np.max(x_all[:,0])])
    ax.set_ylim([np.min(x_all[:,1])-a,a+np.max(x_all[:,1])])
    if grid:
        ax.grid(color=(0,0,0), linestyle='--', linewidth=0.3)
    ax.set_title(title,fontsize=TitleSize)
    ax.set_xlabel(r"$x_{%d}$"%(tuple_of_projection_dimensions[0]+1),fontsize=fontsize)
    ax.set_ylabel(r"$x_{%d}$"%(tuple_of_projection_dimensions[1]+1),fontsize=fontsize)
    if equal_axis:
        ax.axis('equal')
        
def _projection(P,tuple_of_projection_dimensions):
    p_matrix=np.zeros((2,P.n))
    p_matrix[0,tuple_of_projection_dimensions[0]]=1
    p_matrix[1,tuple_of_projection_dimensions[1]]=1
    return pp.affine_map( p_matrix, P)