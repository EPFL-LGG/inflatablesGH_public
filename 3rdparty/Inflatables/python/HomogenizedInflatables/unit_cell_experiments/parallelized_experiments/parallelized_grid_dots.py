import os
import sys; sys.path.append(os.path.join(os.path.dirname(__file__), '..'));

import experiment_helper
from periodic_simulation_setup import *

import parallelism, multiprocessing, itertools, setproctitle
import os, time, numpy as np

import json

grid_size = 10
num_dots = 20

name = "grid_dots_{}_{}".format(grid_size, num_dots)
time_stamp = time.strftime("%Y_%m_%d_%H_%M")
# time_stamp = "2023_11_07_13_58"
base_folder = 'output/{}/{}'.format(name, time_stamp)
if not os.path.exists(base_folder):
    os.makedirs(base_folder)  

useTFT = True
stiffness_pressure = 0.4
scale_factor_pressure = 0.01
avg_len = 0.1

num_experiment = 1000

def run_grid_dots(i):
    result_folder = "{}/{}".format(base_folder, i)
    if not os.path.exists(result_folder):
        os.makedirs(result_folder)  
    if (os.path.isfile('{}/experiment_result.json'.format(result_folder))):
        print("has {}".format(i))
        return
    print("running {}".format(i))

    ipu, m, marker, dots = pattern_generator_using_gmsh.get_random_grid_dots(avg_len, avg_len, grid_size, num_dots)
    np.save("{}/dots_{}_{}.npy".format(result_folder, name, i), dots)
    
    finalMarkers = np.where(np.array(marker) == 1)[0]
    m, finalMarkers = periodic_unit_helper.shift_and_merge_2D_periodic_mesh(m, finalMarkers, flip_orientation= 0)
    m, finalMarkers = periodic_unit_helper.shift_and_merge_2D_periodic_mesh(m, finalMarkers, axis = 1, flip_orientation= 1)

    fusedVtx = get_fusedVtx_using_markers(len(m.vertices()), finalMarkers)

    ipu = inflation.InflatablePeriodicUnit(m, fusedVtx = fusedVtx, epsilon = 1e-9)

    experiment_log = experiment_helper.run_experiment(ipu, m, fusedVtx, stiffness_pressure, scale_factor_pressure, name, i, result_folder, useTFT = useTFT, use_low_pressure=False)

    experiment_log['dots'] = dots.tolist()

    with open('{}/experiment_result.json'.format(result_folder), 'w') as fp:
        json.dump(experiment_log, fp, indent=4)

if __name__ == '__main__':
    num_thread = 1
    # For embarassingly parallel runs, using single thread should be more efficient
    parallelism.set_max_num_tbb_threads(num_thread)

    # Set num logical cores
    logical_core_count = int(multiprocessing.cpu_count() * 0.9)

    print("Using {} logical cores".format(logical_core_count))

    import time
    st = time.time()

    with multiprocessing.Pool(logical_core_count) as p:
        p.map(run_grid_dots, iter(range(num_experiment)))

    # get the end time
    et = time.time()

    # get the execution time
    elapsed_time = et - st
    print('Execution time:', elapsed_time, 'seconds')

    experiment_logs = []
    for i in range(num_experiment):
                result_folder = "{}/{}".format(base_folder, i)
                # Load the experiment log json and add it to the total experiment logs
                with open('{}/experiment_result.json'.format(result_folder), 'r') as fp:
                    curr_json = json.load(fp)
                    curr_json['name'] = i
                    experiment_logs.append(curr_json)


    import json
    with open('{}/experiment_result.json'.format(base_folder), 'w') as fp:
        data = {}
        data['pattern_parameters'] = [
            {
                'name': 'random_count',
                'values': list(range(num_experiment)),
            },
            {
                'name': 'grid_size',
                'values': [grid_size],
            },
            {
                'name': 'num_dots',
                'values': [num_dots],
            }
        ]
        data['configuration'] = {
            'name': name,
            'base_folder': os.path.join(os.path.dirname(os.path.abspath(__file__)), base_folder),
            'time_stamp': time_stamp,
            'useTFT': useTFT,
            'stiffness_pressure': stiffness_pressure,
            'scale_factor_pressure': scale_factor_pressure,
            'avg_len': avg_len,
            'label_format': "{}"
        }
        data['data'] = experiment_logs
        json.dump(data, fp, indent=4)


    # visualize the experiment log json file with histogram plots on its fields