import experiment_helper
import igl
from periodic_simulation_setup import *
import json

allowBending = False

name = 'parallel_tubes'
time_stamp = time.strftime("%Y_%m_%d_%H_%M")
result_folder = 'output/{}/{}'.format(name, time_stamp)
if not os.path.exists(result_folder):
    os.makedirs(result_folder)  

pressure = 0.005

for width in np.linspace(0.1, 1, 10):
    h = 5
    w = width
    print("width", w)

    ipu, m, marker = periodic_unit_helper.get_parallel_tube_periodic(h, w, 0.05)

    finalMarkers = np.where(np.array(marker) == 1)[0]

    fusedVtx = get_fusedVtx_using_markers(len(m.vertices()), finalMarkers)
    ipu = inflation.InflatablePeriodicUnit(m, fusedVtx = fusedVtx, epsilon = 1e-9)

    experiment_helper.run_experiment(ipu, m, fusedVtx, pressure, name, width, allowBending, result_folder)
