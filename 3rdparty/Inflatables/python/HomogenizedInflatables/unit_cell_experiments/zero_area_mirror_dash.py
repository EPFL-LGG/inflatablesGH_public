
import experiment_helper
import igl
from periodic_simulation_setup import *
import json

allowBending = False
useTFT = True

name = "mirror_dash"
time_stamp = time.strftime("%Y_%m_%d_%H_%M")
result_folder = 'output/{}/{}'.format(name, time_stamp)
if not os.path.exists(result_folder):
    os.makedirs(result_folder)  

# pressure = 0.8
stiffness_pressure = 0.4
scale_factor_pressure = 0.01

angles = np.linspace(0, 90, 46)[34:]
radiuses = np.linspace(1, 2.4, 15)

def run_experiment_on_angle_radius(angle, r):
    h = 5
    w = 5
    print("angle", angle, "radius: ", r)
    dash_point = np.array([np.cos(angle / 180 * np.pi), np.sin(angle / 180 * np.pi)]) * r + np.array([0, 0])

    m, marker, n_vx, n_edge = periodic_unit_helper.get_zero_area_dashline(h, w, 0.099, dash_point)

    finalMarkers = np.where(np.array(marker) == 1)[0]
    m, finalMarkers = periodic_unit_helper.shift_and_merge_2D_periodic_mesh(m, finalMarkers, flip_orientation=0)
    m, finalMarkers = periodic_unit_helper.shift_and_merge_2D_periodic_mesh(m, finalMarkers, axis = 1, flip_orientation = 1)


    fusedVtx = get_fusedVtx_using_markers(len(m.vertices()), finalMarkers)
    ipu = inflation.InflatablePeriodicUnit(m, fusedVtx = fusedVtx, epsilon = 1e-9)

    experiment_helper.run_experiment(ipu, m, fusedVtx, stiffness_pressure, scale_factor_pressure, name, "{}_{:.2f}".format(angle, r), allowBending, result_folder, useTFT = useTFT)

# for angle in angles:
#     for r in radiuses:
angle = 42
r = 1.9
h = 5
w = 5
run_experiment_on_angle_radius(angle, r)