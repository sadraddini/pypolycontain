import numpy as np
import warnings

from matplotlib.patches import Polygon
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection

# Pypolycontain
try:
    from pypolycontain.conversions import to_V    
except:
    warnings.warn("You don't have pypolycontain properly installed. Can not import conversions")
    


def visualize(list_of_objects,ax=None,a=0.5,alpha=0.5,list_of_dimensions=[0,1],title=r'pypolycontain visualization'):
    r"""
    Visualization.
    
    inputs: 
        * list_of_objects
    """
    if type(ax)==type(None):
        fig,ax=plt.subplots()
    p_list,x_all=[],np.empty((0,2))  
    for p in list_of_objects:
        x=to_V(p)
        p=Polygon(x)
        p_list.append(p) 
        x_all=np.vstack((x_all,x))
    p_patch = PatchCollection(p_list,color=[p.color for p in list_of_objects], alpha=alpha)
    ax.add_collection(p_patch)
    ax.set_xlim([np.min(x_all[:,0])-a,a+np.max(x_all[:,0])])
    ax.set_ylim([np.min(x_all[:,1])-a,a+np.max(x_all[:,1])])
    ax.grid(color=(0,0,0), linestyle='--', linewidth=0.3)
    ax.set_title(title,FontSize=10)