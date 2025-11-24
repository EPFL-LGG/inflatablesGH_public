import experiment_helper
import igl
from periodic_simulation_setup import *
import json

allowBending = False
useTFT = True

name = 'zero_area_zigzag_dash'
time_stamp = time.strftime("%Y_%m_%d_%H_%M")
result_folder = 'output/{}/{}'.format(name, time_stamp)
if not os.path.exists(result_folder):
    os.makedirs(result_folder)  

stiffness_pressure = 0.5
scale_factor_pressure = 0.01


# for dash_height in np.linspace(2.5, 4.5, 41):
#     h = 10
#     w = 5
#     print("dash_height", dash_height)
#     dash_point = np.array([-2, dash_height])

#     m, marker, n_vx, n_edge = periodic_unit_helper.get_zero_area_dashline(h, w, 0.16, dash_point)

#     finalMarkers = np.where(np.array(marker) == 1)[0]
#     m, finalMarkers = periodic_unit_helper.shift_and_merge_2D_periodic_mesh(m, finalMarkers)
#     m, finalMarkers = periodic_unit_helper.shift_and_merge_2D_periodic_mesh(m, finalMarkers, axis = 1)


#     fusedVtx = get_fusedVtx_using_markers(len(m.vertices()), finalMarkers)
#     ipu = inflation.InflatablePeriodicUnit(m, fusedVtx = fusedVtx, epsilon = 1e-9)

#     experiment_helper.run_experiment(ipu, m, fusedVtx, pressure, name, '%.2f'%dash_height, allowBending, result_folder)


for angle in np.linspace(0, 90, 46):
    h = 10
    w = 5
    print("angle", angle)
    dash_point = np.array([np.cos(angle / 180 * np.pi), np.sin(angle / 180 * np.pi)]) * 2.4 + np.array([0, 2.5])

    m, marker, n_vx, n_edge = periodic_unit_helper.get_zero_area_dashline(h, w, 0.1, dash_point)

    finalMarkers = np.where(np.array(marker) == 1)[0]
    m, finalMarkers = periodic_unit_helper.shift_and_merge_2D_periodic_mesh(m, finalMarkers)
    # m, finalMarkers = periodic_unit_helper.shift_and_merge_2D_periodic_mesh(m, finalMarkers, axis = 1)


    fusedVtx = get_fusedVtx_using_markers(len(m.vertices()), finalMarkers)
    ipu = inflation.InflatablePeriodicUnit(m, fusedVtx = fusedVtx, epsilon = 1e-9)

    experiment_helper.run_experiment(ipu, m, fusedVtx, stiffness_pressure, scale_factor_pressure, name, '%.2f'%angle, allowBending, result_folder, useTFT = useTFT)
