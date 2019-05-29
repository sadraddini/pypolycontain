import numpy as np
from pypolycontain.lib.zonotope import zonotope
from pypolycontain.lib.polytope import polytope
from pypolycontain.lib.AH_polytope import AH_polytope, to_AH_polytope
from time import time
from multiprocessing import Pool


def get_single_random_zonotope(argument):
    dim, generator_range, centroid_range, return_type, color,seed = argument
    np.random.seed(seed)
    m = dim + 1  # np.random.random_integers(dim, 2*dim)
    G = (np.random.rand(dim, m) - 0.5) * generator_range * 1
    x = 2 * (np.random.rand(dim, 1) - 0.5) * centroid_range
    if return_type == 'AH_polytope':
        return to_AH_polytope(zonotope(x, G))
    elif return_type == 'zonotope':
        return zonotope(x, G, color=color)
    else:
        raise NotImplementedError

def get_uniform_random_zonotopes(zonotope_count, dim, centroid_range=None, generator_range = 10, return_type = 'AH_polytope', seed = None, color=None, process_count= 8, as_nparray=True):
    if centroid_range is None:
        centroid_range = zonotope_count*5
    p = Pool(process_count)
    arguments = [[dim, generator_range, centroid_range, return_type, color]]*zonotope_count
    arguments = np.asarray(arguments)
    if seed is None:
        seed = int(time())
    seeds = np.atleast_2d(np.full(zonotope_count, seed)+np.arange(zonotope_count))
    arguments = np.hstack([arguments,seeds.T])
    print('Generating uniform random zonotopes...')
    start_time = time()
    polytopes = p.map(get_single_random_zonotope, arguments)
    print('Completed random zonotope generation in %f seconds' %(time()-start_time))
    if as_nparray:
        return np.asarray(polytopes)
    else:
        return polytopes

def get_uniform_density_random_polytopes(density, dim, centroid_range=None, generator_range = 10, return_type = 'AH_polytope', seed = None, color=None, process_count = 8, as_nparray=True):
    zonotope_count = density**dim
    return get_uniform_random_zonotopes(zonotope_count, dim, centroid_range, generator_range, return_type, seed, color, process_count, as_nparray)


def get_line_random_zonotopes(zonotope_count, dim, centroid_range=None, line_width = None, generator_range = 10, return_type = 'AH_polytope', seed = None, color=None, as_nparray=True):
    if seed is not None:
        np.random.seed(seed)
    polytopes = []
    if centroid_range is None:
        centroid_range = zonotope_count*5
    if line_width is None:
        line_width = 1
    for i in range(zonotope_count):
        m = np.random.random_integers(dim, 2*dim)
        G = (np.random.rand(dim, m) - 0.5) * generator_range * 1
        x = 2*(np.random.rand(dim-1,1) - 0.5)*line_width
        x = np.vstack([2*(np.random.rand(1,1) - 0.5) * centroid_range,x])
        if return_type == 'AH_polytope':
            polytopes.append(to_AH_polytope(zonotope(x, G)))
        elif return_type == 'zonotope':
            polytopes.append(zonotope(x, G, color=color))
        else:
            raise NotImplementedError
    if as_nparray:
        return np.asarray(polytopes)
    return polytopes
