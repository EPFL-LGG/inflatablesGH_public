import experiment_helper
import igl
from periodic_simulation_setup import *

h = 3
w = 1
avg_len = 0.15
shift = [0, 0]
allowBending = True

name = 'dashline_angle'
time_stamp = time.strftime("%Y_%m_%d_%H_%M")
result_folder = 'output/{}/{}'.format(name, time_stamp)
if not os.path.exists(result_folder):
    os.makedirs(result_folder)  

pressure = 2

for angle in np.linspace(30, 45, 31):
    ipu, points, segment_edges, m, markers= periodic_unit_helper.get_shifted_dashline(h, w, avg_len, shift, angle = angle / 180.0 * np.pi, two_dash = False)

    finalMarkers = np.where(np.array(markers) == 1)[0]
    m, finalMarkers = periodic_unit_helper.shift_and_merge_2D_periodic_mesh(m, finalMarkers)
    m, finalMarkers = periodic_unit_helper.shift_and_merge_2D_periodic_mesh(m, finalMarkers, axis = 1)

    fusedVtx = get_fusedVtx_using_markers(len(m.vertices()), finalMarkers)
    ipu = inflation.InflatablePeriodicUnit(m, fusedVtx = fusedVtx, epsilon = 1e-9)

    experiment_helper.run_experiment(ipu, m, fusedVtx, pressure, name, angle, allowBending, result_folder)
