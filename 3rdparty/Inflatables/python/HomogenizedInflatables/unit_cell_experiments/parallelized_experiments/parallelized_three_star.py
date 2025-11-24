

import sys; sys.path.append('../experiments')

import experiment_helper
import igl
from periodic_simulation_setup import *
import json

import parallelism, multiprocessing, itertools, setproctitle
import os, time, numpy as np

allowBending = False
useTFT = True
use_mirror = True 

name = "three_star"
time_stamp = time.strftime("%Y_%m_%d_%H_%M")
result_folder = 'output/{}/{}'.format(name, time_stamp)
if not os.path.exists(result_folder):
    os.makedirs(result_folder)  

# pressure = 0.8
stiffness_pressure = 0.4
scale_factor_pressure = 0.01

radius = np.linspace(1, 2.4, 15)
angles = np.linspace(0, 15, 16)


def run_three_star(r, angle):
    h = 5
    avg_len = 0.1
    dash_point = np.array([np.cos(angle / 180 * np.pi), np.sin(angle / 180 * np.pi)]) * r + np.array([0, 0])
    ipu, m, marker = pattern_generator_using_gmsh.get_three_star(h, avg_len, avg_len, dash_point = dash_point)        

    finalMarkers = np.where(np.array(marker) == 1)[0]
    m, finalMarkers = periodic_unit_helper.shift_and_merge_2D_periodic_mesh(m, finalMarkers, flip_orientation= 0)
    m, finalMarkers = periodic_unit_helper.shift_and_merge_2D_periodic_mesh(m, finalMarkers, axis = 1, flip_orientation= 1)

    fusedVtx = get_fusedVtx_using_markers(len(m.vertices()), finalMarkers)

    ipu = inflation.InflatablePeriodicUnit(m, fusedVtx = fusedVtx, epsilon = 1e-9)

    experiment_helper.run_experiment(ipu, m, fusedVtx, stiffness_pressure, scale_factor_pressure, name, "{:.1f}_{:.2f}".format(r, angle), allowBending, result_folder, useTFT = useTFT, use_low_pressure=False)


def mute():
    "Suppress stdout output"
    sys.stdout = open(os.devnull, 'w')
    

num_thread = 1
# For embarassingly parallel runs, using single thread should be more efficient
parallelism.set_max_num_tbb_threads(num_thread)

# Function arguments
arg1 = radius
arg2 = angles

args = itertools.product(arg1, arg2)

# Set num logical cores
logical_core_count = int(multiprocessing.cpu_count() * 0.9)

print("Using {} logical cores".format(logical_core_count))


import time
st = time.time()

with multiprocessing.Pool(logical_core_count, initializer=mute) as p:
    p.starmap(run_mirror_cosine_dash, args)
        
# get the end time
et = time.time()

# get the execution time
elapsed_time = et - st
print('Execution time:', elapsed_time, 'seconds')