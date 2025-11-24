import experiment_helper
import igl
from periodic_simulation_setup import *
import json

allowBending = False
useTFT = True

name = "grid_dots"
time_stamp = time.strftime("%Y_%m_%d_%H_%M")
result_folder = 'output/{}/{}'.format(name, time_stamp)
if not os.path.exists(result_folder):
    os.makedirs(result_folder)  

# pressure = 0.8
stiffness_pressure = 0.4
scale_factor_pressure = 0.01


num_experiments = 100
for i in range(num_experiments):
    ipu, m, marker, dots = pattern_generator_using_gmsh.get_random_grid_dots(0.1, 0.1, 10)
    np.save("{}/dots_{}_{}.npy".format(result_folder, name, i), dots)
    
    finalMarkers = np.where(np.array(marker) == 1)[0]
    m, finalMarkers = periodic_unit_helper.shift_and_merge_2D_periodic_mesh(m, finalMarkers, flip_orientation= 0)
    m, finalMarkers = periodic_unit_helper.shift_and_merge_2D_periodic_mesh(m, finalMarkers, axis = 1, flip_orientation= 1)

    fusedVtx = get_fusedVtx_using_markers(len(m.vertices()), finalMarkers)

    ipu = inflation.InflatablePeriodicUnit(m, fusedVtx = fusedVtx, epsilon = 1e-9)

    experiment_helper.run_experiment(ipu, m, fusedVtx, stiffness_pressure, scale_factor_pressure, name, i, allowBending, result_folder, useTFT = useTFT, use_low_pressure=False)