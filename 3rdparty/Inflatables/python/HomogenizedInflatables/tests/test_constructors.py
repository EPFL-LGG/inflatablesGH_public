import sys
sys.path.append('../..')
sys.path.append('..')
import inflation, numpy as np, importlib, fd_validation, visualization, parametric_pillows, wall_generation
from numpy.linalg import norm
import MeshFEM, parallelism, benchmark, utils
import periodic_unit_helper
import numpy.linalg as la

def construct_square_pillow_five_vertices():
    triArea = 1
    n_vx, n_edge = periodic_unit_helper.getBox(triArea)
    m, fuseMarkers, fuseSegments = wall_generation.triangulate_channel_walls(n_vx, n_edge, triArea, flags="Y")
    ipu = inflation.InflatablePeriodicUnit(m, np.array(fuseMarkers) != 0, epsilon = 1e-5)
    return ipu

def construct_two_disconnected_square():
    n_vx = [[0, 0], [0, 1], [0, 2],
        [1, 0], [1, 1], [1, 2],
        [2, 0], [2, 1], [2, 2]]
    n_edge = [(0, 1), (1, 2), 
            (3, 4), (4, 5),
            (6, 7), (7, 8),
            (0, 3), (3, 6),
            (1, 4), (4, 7),
            (2, 5), (5, 8)]
    triArea = 0.5
    m, fuseMarkers, fuseSegments = wall_generation.triangulate_channel_walls(n_vx, n_edge, triArea, flags="Y")
    fuseMarkers = [0] * 9
    ipu = inflation.InflatablePeriodicUnit(m, np.array(fuseMarkers) != 0, epsilon = 1e-5)
    return ipu

def construct_two_squares_connected_at_center():
    n_vx = [[0, 0], [0, 1], [0, 2],
        [1, 0], [1, 1], [1, 2],
        [2, 0], [2, 1], [2, 2]]
    n_edge = [(0, 1), (1, 2), 
            (3, 4), (4, 5),
            (6, 7), (7, 8),
            (0, 3), (3, 6),
            (1, 4), (4, 7),
            (2, 5), (5, 8)]
    triArea = 0.5
    m, fuseMarkers, fuseSegments = wall_generation.triangulate_channel_walls(n_vx, n_edge, triArea, flags="Y")
    fuseMarkers = [0] * 9
    fuseMarkers[4] = 1
    ipu = inflation.InflatablePeriodicUnit(m, np.array(fuseMarkers) != 0, epsilon = 1e-5)
    return ipu
    
def construct_rectangle_connected_at_one_side():
    n_vx = [[0, 0], [0, 1], [0, 2], [0, 3],
        [1, 0], [1, 1], [1, 2], [1, 3]]
    n_edge = [(0, 1), (1, 2), (2, 3), 
            (4, 5), (5, 6), (6, 7),
            (0, 4), (1, 5), (2, 6), (3, 7),
            (1, 4), (2, 5), (3, 6)]
    triArea = 1
    m, fuseMarkers, fuseSegments = wall_generation.triangulate_channel_walls(n_vx, n_edge, triArea, flags="Y")
    fuseMarkers = [0] * 8
    fuseMarkers[0] = 1
    fuseMarkers[1] = 1
    fuseMarkers[3] = 1
    fuseMarkers[4] = 1
    fuseMarkers[5] = 1
    fuseMarkers[7] = 1
    ipu = inflation.InflatablePeriodicUnit(m, np.array(fuseMarkers) == 1)
    return ipu

def construct_dashline():
    h_big = 5
    w_big = 5
    h_small = 3
    w_small = 3
    avg_len = 0.3
    ipu, points, segment_edges, m, marker= periodic_unit_helper.get_boundary_aligned_dashline(w_small, w_big, h_small, h_big, avg_len)
    return ipu

def test_ipu_constructors():
    construct_square_pillow_five_vertices()
    construct_two_disconnected_square()
    construct_two_squares_connected_at_center()
    construct_dashline()


if __name__ == "__main__":
    construct_square_pillow_five_vertices()
    construct_two_disconnected_square()
    construct_two_squares_connected_at_center()