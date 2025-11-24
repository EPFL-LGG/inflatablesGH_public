import experiment_helper
import igl
from periodic_simulation_setup import *
import json

allowBending = False

name = 'dash_line'
time_stamp = time.strftime("%Y_%m_%d_%H_%M")
result_folder = 'output/{}/{}'.format(name, time_stamp)
if not os.path.exists(result_folder):
    os.makedirs(result_folder)  

pressure = 0.8
pressure = 0.3

labels = np.linspace(0.1, 0.05, 21)
for label in labels:
    h = 4
    w = 0.5
    res = 100
    triArea = h * w / res
    avg_len = label

    shift = np.array([0., 0.])
    ipu, points, segment_edges, m, marker= periodic_unit_helper.get_shifted_dashline(h, w, avg_len, shift, False, angle = 0)

    finalMarkers = np.where(np.array(marker) == 1)[0]
    # m, finalMarkers = periodic_unit_helper.shift_and_merge_2D_periodic_mesh(m, finalMarkers)
    # m, finalMarkers = periodic_unit_helper.shift_and_merge_2D_periodic_mesh(m, finalMarkers, axis = 1)

    fusedVtx = get_fusedVtx_using_markers(len(m.vertices()), finalMarkers)
    ipu = inflation.InflatablePeriodicUnit(m, fusedVtx = fusedVtx, epsilon = 1e-9)

    experiment_helper.run_experiment(ipu, m, fusedVtx, pressure, name, '%.4f'%label, allowBending, result_folder)
