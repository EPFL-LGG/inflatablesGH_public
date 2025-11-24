import os
import sys; sys.path.append(os.path.join(os.path.dirname(__file__), '..'));

import experiment_helper
from periodic_simulation_setup import *

import parallelism, multiprocessing, itertools, setproctitle
import os, time, numpy as np

import json

name = "square_with_ellipse_hole"
time_stamp = time.strftime("%Y_%m_%d_%H_%M")
time_stamp = "2024_01_13_21_03"
base_folder = 'output/{}/{}'.format(name, time_stamp)
if not os.path.exists(base_folder):
    os.makedirs(base_folder)  

useTFT = True
stiffness_pressure = 0.4
scale_factor_pressure = 0.3
avg_len = 0.08

def get_label(radius, angle):
    return "{:.2f}_{:.2f}".format(radius, angle)

def run_square_with_ellipse_holes(radius, angle):
    iw = np.cos(angle / 180 * np.pi) * radius
    ih = np.sin(angle / 180 * np.pi) * radius
   
    result_folder = "{}/{}".format(base_folder, get_label(radius, angle))
    if not os.path.exists(result_folder):
        os.makedirs(result_folder)  
    if (os.path.isfile('{}/experiment_result.json'.format(result_folder))):
        data = json.load(open('{}/experiment_result.json'.format(result_folder)))
        if data["Planar equilibrium"] == 1:
            return
    print("running {:.2f}_{:.2f}".format(radius, angle))

    h = 5
    m, marker = pattern_generator_using_gmsh.get_pattern_with_elliptic_holes(h, iw, ih, avg_len)

    finalMarkers = np.where(np.array(marker) == 1)[0]
    m, finalMarkers = periodic_unit_helper.shift_and_merge_2D_periodic_mesh(m, finalMarkers, flip_orientation= 0)
    m, finalMarkers = periodic_unit_helper.shift_and_merge_2D_periodic_mesh(m, finalMarkers, axis = 1, flip_orientation= 1)

    fusedVtx = get_fusedVtx_using_markers(len(m.vertices()), finalMarkers)

    ipu = inflation.InflatablePeriodicUnit(m, fusedVtx = fusedVtx, epsilon = 1e-9)

    experiment_log = experiment_helper.run_experiment(ipu, m, fusedVtx, stiffness_pressure, scale_factor_pressure, name, get_label(radius, angle), result_folder, useTFT = useTFT, use_low_pressure=False)

    with open('{}/experiment_result.json'.format(result_folder), 'w') as fp:
        json.dump(experiment_log, fp, indent=4)


if __name__ == '__main__':
    angles = np.linspace(47, 59, 7)
    radius = np.linspace(0.6, 2, 8)[1: -1]

    num_thread = 1
    # For embarassingly parallel runs, using single thread should be more efficient
    parallelism.set_max_num_tbb_threads(num_thread)

    # Function arguments
    arg1 = radius
    arg2 = angles

    args = itertools.product(arg1, arg2)

    # Set num logical cores
    logical_core_count = int(multiprocessing.cpu_count() * 0.4)

    print("Using {} logical cores".format(logical_core_count))

    import time
    st = time.time()

    with multiprocessing.Pool(logical_core_count) as p:
        p.starmap(run_square_with_ellipse_holes, args)
            
    # get the end time
    et = time.time()

    # get the execution time
    elapsed_time = et - st
    print('Execution time:', elapsed_time, 'seconds')

    experiment_logs = []
    for a in angles:
        for r in radius:
            result_folder = "{}/{}".format(base_folder, get_label(r, a))
            # Load the experiment log json and add it to the total experiment logs
            with open('{}/experiment_result.json'.format(result_folder), 'r') as fp:
                curr_json = json.load(fp)
                curr_json['name'] = get_label(r, a)
                experiment_logs.append(curr_json)


    import json
    with open('{}/experiment_result.json'.format(base_folder), 'w') as fp:
        data = {}
        data['pattern_parameters'] = [
            {
                'name': 'radius',
                'values': radius.tolist(),
            },
            {
                'name': 'angles',
                'values': angles.tolist(),
            },
        ]
        data['configuration'] = {
            'name': name,
            'base_folder': os.path.join(os.path.dirname(os.path.abspath(__file__)), base_folder),
            'time_stamp': time_stamp,
            'useTFT': useTFT,
            'stiffness_pressure': stiffness_pressure,
            'scale_factor_pressure': scale_factor_pressure,
            'avg_len': avg_len,
            'label_format': "{:.2f}_{:.2f}"
        }
        data['data'] = experiment_logs
        data['elapsed_time'] = elapsed_time
        json.dump(data, fp, indent=4)


    # visualize the experiment log json file with histogram plots on its fields
