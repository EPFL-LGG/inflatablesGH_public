import sys
sys.path.append('../../3rdparty/elastic_rods/3rdparty/MeshFEM/python')  # path to MeshFEM to fetch parallelism.cpython

import parallelism, multiprocessing, itertools, setproctitle
import os, time, numpy as np


def function(arg1, arg2):
    "To be run in parallel"

    # Setting the process name allows to identify it in htop
    script_name = os.path.basename(__file__)
    process_name = '{}-{}-{}'.format(script_name, arg1, arg2)
    setproctitle.setproctitle(process_name)

    time.sleep(1)

    return


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

angles = np.linspace(0, 90, 46)
radiuses = np.linspace(1, 2.4, 15)


def run_mirror_dash(angle, r):
    h = 5
    w = 5
    print("angle", angle, "radius: ", r)
    dash_point = np.array([np.cos(angle / 180 * np.pi), np.sin(angle / 180 * np.pi)]) * r + np.array([0, 0])

    m, marker, n_vx, n_edge = periodic_unit_helper.get_zero_area_dashline(h, w, 0.05, dash_point)

    finalMarkers = np.where(np.array(marker) == 1)[0]
    m, finalMarkers = periodic_unit_helper.shift_and_merge_2D_periodic_mesh(m, finalMarkers, flip_orientation=0)
    m, finalMarkers = periodic_unit_helper.shift_and_merge_2D_periodic_mesh(m, finalMarkers, axis = 1, flip_orientation = 1)


    fusedVtx = get_fusedVtx_using_markers(len(m.vertices()), finalMarkers)
    ipu = inflation.InflatablePeriodicUnit(m, fusedVtx = fusedVtx, epsilon = 1e-9)

    experiment_helper.run_experiment(ipu, m, fusedVtx, stiffness_pressure, scale_factor_pressure, name, "{}_{:.2f}".format(angle, r), allowBending, result_folder, useTFT = useTFT)


def mute():
    "Suppress stdout output"
    sys.stdout = open(os.devnull, 'w')
    

# For embarassingly parallel runs, using single thread should be more efficient
parallelism.set_max_num_tbb_threads(1)

# Function arguments
arg1 = angles
arg2 = radiuses
n1 = len(arg1)
n2 = len(arg2)
args = itertools.product(arg1, arg2)

# Set num logical cores
logical_core_count = multiprocessing.cpu_count()
logical_cores_used = min(n1*n2, logical_core_count)
print("Using {} logical cores".format(logical_cores_used))

# Run in parallel
print("Running {}*{}={} instances of {}".format(n1, n2, n1*n2, function.__name__))
with multiprocessing.Pool(logical_cores_used, initializer=mute) as p:
    p.starmap(run_mirror_dash, args)
