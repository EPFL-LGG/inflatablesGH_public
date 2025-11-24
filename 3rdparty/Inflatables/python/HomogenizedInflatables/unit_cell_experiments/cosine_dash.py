import experiment_helper
import igl
from periodic_simulation_setup import *
import json

allowBending = False
useTFT = True
useMirror = True

name = 'cosine_dash_{}'.format("mirror" if useMirror else "no_mirror")
time_stamp = time.strftime("%Y_%m_%d_%H_%M")
result_folder = 'output/{}/{}'.format(name, time_stamp)
if not os.path.exists(result_folder):
    os.makedirs(result_folder)  

# pressure = 0.8
stiffness_pressure = 0.4
scale_factor_pressure = 0.01

amplitudes = np.linspace(-1, 1, 21)
radius = np.linspace(1, 1.5, 6)
angles = np.linspace(0, 45, 16)

for amp in amplitudes:
    for r in radius:
        for angle in angles:
            h = 5
            avg_len = 0.1
            dash_point = np.array([np.cos(angle / 180 * np.pi), np.sin(angle / 180 * np.pi)]) * r + np.array([0, 0])
            ipu, m, marker = pattern_generator_using_gmsh.get_cosine_dash(h, avg_len, avg_len, amplitude=amp, dash_point = dash_point)
            
            finalMarkers = np.where(np.array(marker) == 1)[0]
            m, finalMarkers = periodic_unit_helper.shift_and_merge_2D_periodic_mesh(m, finalMarkers, flip_orientation= 0)
            m, finalMarkers = periodic_unit_helper.shift_and_merge_2D_periodic_mesh(m, finalMarkers, axis = 1, flip_orientation= 1)

            fusedVtx = get_fusedVtx_using_markers(len(m.vertices()), finalMarkers)

            ipu = inflation.InflatablePeriodicUnit(m, fusedVtx = fusedVtx, epsilon = 1e-9)

            experiment_helper.run_experiment(ipu, m, fusedVtx, stiffness_pressure, scale_factor_pressure, name, "{:.1f}_{:.1f}_{:.2f}".format(amp, r, angle), allowBending, result_folder, useTFT = useTFT)