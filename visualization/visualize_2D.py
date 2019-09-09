# -*- coding: utf-8 -*-
"""
Created on Wed Oct 24 16:35:05 2018

@author: sadra

This part is only for visualization of 2D Polytopes
"""

from matplotlib.patches import Polygon
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from scipy.spatial import ConvexHull
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

from pypolycontain.lib.operations import AH_polytope_vertices    


import numpy as np
try:
    from cdd import Polyhedron,Matrix,RepType
except:
    print "WARNING: You don't have CDD package installed. Unable to visualize polytopes. You may still visualize zonotopes."

def visualize_2D(list_of_polytopes,a=1.5,title="polytopes"):
    """
    Given a polytope in its H-representation, plot it
    """ 
    ana_color=[(1,0,0),(0,1,0),(0,0,1),(1,1,0),(1,0,1),(0,0,1),(0,0,1)]
    p_list=[]
    x_all=np.empty((0,2))
    for polytope in list_of_polytopes:
        p_mat=Matrix(np.hstack((polytope.h,-polytope.H)))
        p_mat.rep_type = RepType.INEQUALITY
        poly=Polyhedron(p_mat)
        y=np.array(poly.get_generators())
        x=y[:,1:3]#/y[:,2].reshape(y.shape[0],1)
        x=x[ConvexHull(x).vertices,:]
        x_all=np.vstack((x_all,x))
        p=Polygon(x)
        p_list.append(p)
#    p_patch = PatchCollection(p_list, color=[(np.random.random(),np.random.random(),np.tanh(np.random.random())) \
#        for polytope in list_of_polytopes],alpha=0.7)
    p_patch = PatchCollection(p_list,color=ana_color[0:len(list_of_polytopes)], alpha=0.5)
    fig, ax = plt.subplots()
    ax.add_collection(p_patch)
    ax.set_xlim([np.min(x_all[:,0])-a,a+np.max(x_all[:,0])])
    ax.set_ylim([np.min(x_all[:,1])-a,a+np.max(x_all[:,1])])
    ax.grid(color=(0,0,0), linestyle='--', linewidth=0.3)
    ax.set_title(title)
    return fig


    
def visualize_2D_zonotopes(list_of_zonotopes,a=1.5,list_of_dimensions=None,title="zonotopes",axis_limit=[0],\
                           alpha = 0.975,fig=None, ax = None, color = None):
    """
    Given a list of zonotopes, draw them. The zonotopes already have colors.
    """       
    if type(list_of_dimensions)==type(None):
        list_of_dimensions=[0,1]
    p_list=[]
    x_all=np.empty((0,2))
    for zono in list_of_zonotopes:
        y=zono.x.T+np.dot(zono.G,vcube(zono.G.shape[1]).T).T
        x=y[:,list_of_dimensions]#/y[:,2].reshape(y.shape[0],1)
        x=x[ConvexHull(x).vertices,:]
        x_all=np.vstack((x_all,x))
        p=Polygon(x)
        p_list.append(p)
    if color is None:
        p_patch = PatchCollection(p_list, color=[Z.color for Z in list_of_zonotopes],alpha=alpha)
    else:
        p_patch = PatchCollection(p_list, color=color,alpha=alpha)
#    p_patch = PatchCollection(p_list, color=[(1-zono.x[0,0]>=1,0,zono.x[0,0]>=1) \
#        for zono in list_of_zonotopes],alpha=0.75)
    if fig is None or ax is None:
        fig, ax = plt.subplots()
    ax.add_collection(p_patch)
#    print(axis_limit)
    if len(axis_limit)==1:
        ax.set_xlim([np.min(x_all[:,0])-a,a+np.max(x_all[:,0])])
        ax.set_ylim([np.min(x_all[:,1])-a,a+np.max(x_all[:,1])])
    else:
        ax.set_xlim([axis_limit[0],axis_limit[1]])
        ax.set_ylim([axis_limit[2],axis_limit[3]])
    ax.grid(color=(0,0,0), linestyle='--', linewidth=0.3)
    ax.set_title(title)
    return fig,ax

    
def visualize_2D_zonotopes_convexhull(fig,ax,list_of_zonotopes,a=1.5,list_of_dimensions=None,title="zonotopes",axis_limit=[0]):
    if type(list_of_dimensions)==type(None):
        list_of_dimensions=[0,1]
    p_list=[]
    x_all=np.empty((0,2))
    for zono in list_of_zonotopes:
        y=zono.x.T+np.dot(zono.G,vcube(zono.G.shape[1]).T).T
        x=y[:,list_of_dimensions]#/y[:,2].reshape(y.shape[0],1)
        x_all=np.vstack((x_all,x))
    x_all=x_all[ConvexHull(x_all).vertices,:]
    p=Polygon(x_all)
    p_list.append(p)
    p_patch = PatchCollection(p_list, color=(0.5,0.5,0.5),alpha=0.75)
#    p_patch = PatchCollection(p_list, color=[(1-zono.x[0,0]>=1,0,zono.x[0,0]>=1) \
#        for zono in list_of_zonotopes],alpha=0.75)
    ax.add_collection(p_patch)
    if len(axis_limit)==1:
        ax.set_xlim([np.min(x_all[:,0])-a,a+np.max(x_all[:,0])])
        ax.set_ylim([np.min(x_all[:,1])-a,a+np.max(x_all[:,1])])
    else:
        ax.set_xlim([axis_limit[0],axis_limit[1]])
        ax.set_ylim([axis_limit[2],axis_limit[3]])
    ax.grid(color=(0,0,0), linestyle='--', linewidth=0.3)
    ax.set_title(title)

def visualize_3D_zonotopes(list_of_zonotopes,a=1.5,list_of_dimensions=None):
    """
    Given a polytope in its H-representation, plot it
    """ 
    if type(list_of_dimensions)==type(None):
        list_of_dimensions=[0,1,2]
    p_list=[]
    x_all=np.empty((0,3))
    for zono in list_of_zonotopes:
        y=np.dot(zono.G,vcube(zono.G.shape[1]).T).T
        x=y[:,list_of_dimensions]#/y[:,2].reshape(y.shape[0],1)
        x=x[ConvexHull(x).vertices,:]
        p_mat=Matrix(x)
        p_mat.rep_type = RepType.GENERATOR
        x_all=np.vstack((x_all,x))
        p=Poly3DCollection([x])
        p_list.append(p)
#    p_patch = PatchCollection(p_list, color=[(np.random.random(),np.random.random(),np.tanh(np.random.random())) \
#        for zono in list_of_zonotopes],alpha=0.6)
    fig = plt.figure()
    ax = Axes3D(fig)
    ax.set_xlim3d([np.min(x_all[:,0])-a,a+np.max(x_all[:,0])])
    ax.set_ylim3d([np.min(x_all[:,1])-a,a+np.max(x_all[:,1])])
    ax.set_zlim3d([np.min(x_all[:,2])-a,a+np.max(x_all[:,2])])
    ax.add_collection3d(p)
    ax.grid3d(color=(0,0,0), linestyle='--', linewidth=0.3)
    
def visualize_2D_AH_polytope(list_of_AH_polytopes,a=1.5,color=None,alpha=0.5,fig=None,ax=None,axis_limit=[0],title=r"AH-Polytopes",N=20,epsilon=0.001):
    p_list=[]
    v_all=np.empty((0,2))
    for Q in list_of_AH_polytopes:
        v,w=AH_polytope_vertices(Q,N=N,epsilon=epsilon)
        # Minkowski sum with epsilon ball
        try:
            v=v[ConvexHull(v).vertices,:]
        except:
            print "WARNING:", Q,": was degenerate or very thin to plot. Adding a tube of", epsilon,"N:", v.shape
            v=w[ConvexHull(w).vertices,:]
        v_all=np.vstack((v_all,v))
        p=Polygon(v)
        p_list.append(p)
    if color is None:
        p_patch = PatchCollection(p_list, color=[Z.color for Z in list_of_AH_polytopes],alpha=alpha)
    else:
        p_patch = PatchCollection(p_list, color=color,alpha=alpha)
#    p_patch = PatchCollection(p_list, color=[(1-zono.x[0,0]>=1,0,zono.x[0,0]>=1) \
#        for zono in list_of_zonotopes],alpha=0.75)
    if fig is None or ax is None:
        fig, ax = plt.subplots()
    ax.add_collection(p_patch)
#    print(axis_limit)
    if len(axis_limit)==1:
        ax.set_xlim([np.min(v_all[:,0])-a,a+np.max(v_all[:,0])])
        ax.set_ylim([np.min(v_all[:,1])-a,a+np.max(v_all[:,1])])
    else:
        ax.set_xlim([axis_limit[0],axis_limit[1]])
        ax.set_ylim([axis_limit[2],axis_limit[3]])
    ax.grid(color=(0,0,0), linestyle='--', linewidth=0.3)
    ax.set_title(title)
    return fig,ax



"""
The following functions involve the ax object, or the plot, as one of the arguments
"""

def visualize_2D_ax(ax,list_of_polytopes,a=1.5,title="polytopes",color=[True],alpha=0.5):
    """
    Given a polytope in its H-representation, plot it
    """ 
    if any(color)==True:
        ana_color=[(0,1,0),(0,0,1),(1,1,0),(1,0,1),(0,0,1),(0,0,1),(1,0,0),
                   (0.5,0.5,0),(0,0.5,0.5),(0.5,0,0.5),(0.5,0.5,0.5),
                   (1,1,0.5),(1,0.5,1),(0.5,1,1),
                   (0.5,0.5,1),(0.5,1,0.5),(1,0.5,0.5)]*len(list_of_polytopes)
    else:
        ana_color=[color]*len(list_of_polytopes)
    p_list=[]
    x_all=np.empty((0,2))
    for polytope in list_of_polytopes:
        p_mat=Matrix(np.hstack((polytope.h,-polytope.H)))
        p_mat.rep_type = RepType.INEQUALITY
        poly=Polyhedron(p_mat)
        y=np.array(poly.get_generators())
        x=y[:,1:3]#/y[:,2].reshape(y.shape[0],1)
        x=x[ConvexHull(x).vertices,:]
        x_all=np.vstack((x_all,x))
        p=Polygon(x)
        p_list.append(p)
    p_patch = PatchCollection(p_list,color=ana_color[0:len(list_of_polytopes)], alpha=alpha)
    ax.add_collection(p_patch)
    ax.set_xlim([np.min(x_all[:,0])-a,a+np.max(x_all[:,0])])
    ax.set_ylim([np.min(x_all[:,1])-a,a+np.max(x_all[:,1])])
    ax.grid(color=(0,0,0), linestyle='--', linewidth=0.3)
    ax.set_title(title)


def visualize_2D_zonotopes_ax(ax,list_of_zonotopes,a=1.5,list_of_dimensions=None,title="zonotopes",axis_limit=[0],alpha=0.5):
    """
    Given a plot, add zonotopes
    """ 
    print("*"*30,"Getting a plot of your zonotopes, be patient!")
    print(ax)
    if type(list_of_dimensions)==type(None):
        list_of_dimensions=[0,1]
    p_list=[]
    x_all=np.empty((0,2))
    for zono in list_of_zonotopes:
        y=zono.x.T+np.dot(zono.G,vcube(zono.G.shape[1]).T).T
        x=y[:,list_of_dimensions]#/y[:,2].reshape(y.shape[0],1)
        x=x[ConvexHull(x).vertices,:]
        x_all=np.vstack((x_all,x))
        p=Polygon(x)
        p_list.append(p)
    p_patch = PatchCollection(p_list, color=[Z.color for Z in list_of_zonotopes],alpha=alpha)
#    p_patch = PatchCollection(p_list, color=[(1-zono.x[0,0]>=1,0,zono.x[0,0]>=1) \
#        for zono in list_of_zonotopes],alpha=0.75)
    ax.add_collection(p_patch)
    if len(axis_limit)==1:
        ax.set_xlim([np.min(x_all[:,0])-a,a+np.max(x_all[:,0])])
        ax.set_ylim([np.min(x_all[:,1])-a,a+np.max(x_all[:,1])])
    else:
        ax.set_xlim([axis_limit[0],axis_limit[1]])
        ax.set_ylim([axis_limit[2],axis_limit[3]])
    ax.grid(color=(0,0,0), linestyle='--', linewidth=0.3)
    ax.set_title(title)


"""
Auxilary functions
"""

    
def vcube(T):
    """
    Description: 2**n * n array of vectors of vertices in unit cube in R^n
    """
    from itertools import product 
    v=list(product(*zip([-1]*T,[1]*T)))
    return np.array(v)