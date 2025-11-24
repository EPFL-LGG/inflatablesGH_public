import experiment_helper
import igl
from periodic_simulation_setup import *
import json
from periodic_unit_helper import get_scaled_box, point_in_box, point_in_line_segment
allowBending = False

name = 'resolution_regular_grid'
time_stamp = time.strftime("%Y_%m_%d_%H_%M")
result_folder = 'output/{}/{}'.format(name, time_stamp)
if not os.path.exists(result_folder):
    os.makedirs(result_folder)  



disableFusedRegionTFT = False

use_square = True

if use_square:
    stiffness_pressure = 1
    scale_factor_pressure = 1
else:
    stiffness_pressure = 0.4
    scale_factor_pressure = 0.01

for factor in ((2**np.arange(10)[2:] + 1) if use_square else (np.arange(48)[4:])):
    if (use_square):
        resolution = int(5 / 0.7) * factor + 1
    else:
        resolution = factor * 2 - 1
    print(resolution)
    x = np.linspace(0, 5, resolution)
    y = np.linspace(0, 5, resolution)
    pts = np.transpose([np.tile(x, len(y)), np.repeat(y, len(x))])
    edges = []
    for j in range(resolution):
        for i in range(resolution-1):
            edges.append([i + resolution * j, i + 1 + resolution * j])
            
    for j in range(resolution ):
        for i in range(resolution - 1):
            edges.append([j + resolution * i, j + resolution * (i + 1)])
            
    for j in range(resolution-1):
        for i in range(resolution - 1):
            edges.append([j + resolution * i, j + resolution * (i + 1) + 1])
    triArea = 100
    m, fuseMarkers, fuseSegments = wall_generation.triangulate_channel_walls(pts, edges, triArea, flags="Y")
    use_boundary_aligned = False
    box = [[1, 1], [4, 4]]
    wall_box_1 = get_scaled_box(np.array(box), 1)
    wall_box_2 = get_scaled_box(np.array(box) - [1, 1], 1)
    wall_box_3 = get_scaled_box(np.array(box) - [1, -4], 1)
    wall_box_4 = get_scaled_box(np.array(box) - [-4, 1], 1)
    wall_box_5 = get_scaled_box(np.array(box) - [-4, -4], 1)

    marker1 = [point_in_box(pt, wall_box_1) for pt in m.vertices()]

    marker2 = [point_in_box(pt, wall_box_2) or point_in_box(pt, wall_box_3) or point_in_box(pt, wall_box_4) or point_in_box(pt, wall_box_5) for pt in m.vertices()]

    marker3 = [point_in_line_segment(pt[:2], np.array([[2.5, 1], [2.5, 4]])) for pt in m.vertices()]
    if (use_square):
        finalMarkers = np.where(np.array(marker2 if use_boundary_aligned else marker1) == 1)[0]
    else:
        finalMarkers = np.where(np.array(marker3) == 1)[0]

    if (use_square):
        m, finalMarkers = periodic_unit_helper.shift_and_merge_2D_periodic_mesh(m, finalMarkers)
        m, finalMarkers = periodic_unit_helper.shift_and_merge_2D_periodic_mesh(m, finalMarkers, axis = 1)

    fusedVtx = get_fusedVtx_using_markers(len(m.vertices()), finalMarkers)

    ipu = inflation.InflatablePeriodicUnit(m, fusedVtx = fusedVtx, epsilon = 1e-9)

    experiment_helper.run_experiment(ipu, m, fusedVtx, stiffness_pressure, scale_factor_pressure, name, factor, allowBending, result_folder)
