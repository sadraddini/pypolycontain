from pypolycontain.utils.random_polytope_generator import *
from pypolycontain.lib.AH_polytope import distance_point

def test_get_k_random_vertices_in_zonotope():
    zonotope = get_uniform_random_zonotopes(1,4, return_type='zonotope')[0]
    for p in get_k_random_edge_points_in_zonotope(zonotope, 3):
        assert(distance_point(zonotope,p)[0]<1e-3)
    return

if __name__ == '__main__':
    test_get_k_random_vertices_in_zonotope()
