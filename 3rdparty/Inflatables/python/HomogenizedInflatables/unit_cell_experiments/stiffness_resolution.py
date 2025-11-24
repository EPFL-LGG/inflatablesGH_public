import experiment_helper
import igl
from periodic_simulation_setup import *
import json

allowBending = False

name = 'stiffness_resolution'
time_stamp = time.strftime("%Y_%m_%d_%H_%M")
result_folder = 'output/{}/{}'.format(name, time_stamp)
if not os.path.exists(result_folder):
    os.makedirs(result_folder)  

pressure = 0.8
pressure = 2

disableFusedRegionTFT = False
# # boundary aligned dashline
# for res in np.linspace(0.03, 0.2, 18):
#     h_big = 5
#     w_big = 5

#     h_small = 3
#     w_small = 3

#     avg_len = res
#     ipu, points, segment_edges, m, marker= periodic_unit_helper.get_boundary_aligned_dashline(w_small, w_big, h_small, h_big, avg_len)

#     finalMarkers = np.where(np.array(marker) == 1)[0]
#     m, finalMarkers = periodic_unit_helper.shift_and_merge_2D_periodic_mesh(m, finalMarkers)
#     m, finalMarkers = periodic_unit_helper.shift_and_merge_2D_periodic_mesh(m, finalMarkers, axis = 1)

#     fusedVtx = get_fusedVtx_using_markers(len(m.vertices()), finalMarkers)
#     ipu = inflation.InflatablePeriodicUnit(m, fusedVtx = fusedVtx, epsilon = 1e-9)

#     experiment_helper.run_experiment(ipu, m, fusedVtx, pressure, name, res, allowBending, result_folder)

# two dash shifted and opposite angle
for res in np.linspace(0.03, 0.2, 18):
    h = 2
    w = 0.5
    avg_len = res

    shift = np.array([1.2, 0.])

    ipu, points, segment_edges, m, marker = periodic_unit_helper.get_shifted_dashline(h, w, avg_len, shift, angle = 37 / 180 * np.pi, two_dash = True, opposite_angle = True)

    finalMarkers = np.where(np.array(marker) == 1)[0]
    m, finalMarkers = periodic_unit_helper.shift_and_merge_2D_periodic_mesh(m, finalMarkers)
    m, finalMarkers = periodic_unit_helper.shift_and_merge_2D_periodic_mesh(m, finalMarkers, axis = 1)

    fusedVtx = get_fusedVtx_using_markers(len(m.vertices()), finalMarkers)
    ipu = inflation.InflatablePeriodicUnit(m, fusedVtx = fusedVtx, epsilon = 1e-9)

    experiment_helper.run_experiment(ipu, m, fusedVtx, pressure, name, res, allowBending, result_folder)
