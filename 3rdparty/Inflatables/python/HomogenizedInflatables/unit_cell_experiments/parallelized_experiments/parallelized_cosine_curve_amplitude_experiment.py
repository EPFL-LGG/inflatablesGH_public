import os
import sys; sys.path.append(os.path.join(os.path.dirname(__file__), '..'));

import experiment_helper
from periodic_simulation_setup import *
import json

import parallelism, multiprocessing, itertools, setproctitle
import os, time, numpy as np

allowBending = False
useTFT = True

use_half_period = False

name = 'cosine_curve_amplitude_{}'.format('half_period' if use_half_period else 'full_period')
time_stamp = time.strftime("%Y_%m_%d_%H_%M")
# time_stamp = '2023_12_04_10_34'
base_folder = 'output/{}/{}'.format(name, time_stamp)
if not os.path.exists(base_folder):
    os.makedirs(base_folder)  

# pressure = 0.8
stiffness_pressure = 0.3
scale_factor_pressure = 0.1
avg_len = 0.07

num_experiment = 50
h = 5
amplitudes = np.linspace(0, 4.5 * (0.5 if use_half_period else 1) / h, num_experiment)

def run_cosine_curve_experiment(label):
    result_folder = "{}/{}".format(base_folder, label)
    if not os.path.exists(result_folder):
        os.makedirs(result_folder)  

    if (os.path.isfile('{}/experiment_result.json'.format(result_folder))):
        print("has {}".format(label))
        return
    print("running {}".format(label))

    m, marker = pattern_generator_using_gmsh.get_cosine_curve(h, avg_len, avg_len, amplitude=amplitudes[label], end_threshold = 0.0, use_half_period=use_half_period)
    if (m is None):
        return
    finalMarkers = np.where(np.array(marker) == 1)[0]
    m, finalMarkers = periodic_unit_helper.shift_and_merge_2D_periodic_mesh(m, finalMarkers, flip_orientation = 0)
    m, finalMarkers = periodic_unit_helper.shift_and_merge_2D_periodic_mesh(m, finalMarkers, axis = 1, flip_orientation = 1)

    fusedVtx = get_fusedVtx_using_markers(len(m.vertices()), finalMarkers)

    ipu = inflation.InflatablePeriodicUnit(m, fusedVtx = fusedVtx, epsilon = 1e-9)

    experiment_log = experiment_helper.run_experiment(ipu, m, fusedVtx, stiffness_pressure, scale_factor_pressure, name, label, result_folder, useTFT = useTFT, use_low_pressure=False, allow_bending=False)

    with open('{}/experiment_result.json'.format(result_folder), 'w') as fp:
        json.dump(experiment_log, fp, indent=4)



if __name__ == '__main__':
    labels = np.arange(num_experiment)

    num_thread = 1
    # For embarassingly parallel runs, using single thread should be more efficient
    parallelism.set_max_num_tbb_threads(num_thread)

    # Set num logical cores
    logical_core_count = int(multiprocessing.cpu_count() * 0.5)

    print("Using {} logical cores".format(logical_core_count))

    import time
    st = time.time()

    with multiprocessing.Pool(logical_core_count) as p:
        p.map(run_cosine_curve_experiment, iter(labels))
            
    # get the end time
    et = time.time()

    # get the execution time
    elapsed_time = et - st
    print('Execution time:', elapsed_time, 'seconds')

    experiment_logs = []
    for label in labels:
        result_folder = "{}/{}".format(base_folder, label)
        # Load the experiment log json and add it to the total experiment logs
        with open('{}/experiment_result.json'.format(result_folder), 'r') as fp:
            curr_json = json.load(fp)
            curr_json['name'] = label
            experiment_logs.append(curr_json)


    import json
    with open('{}/experiment_result.json'.format(base_folder), 'w') as fp:
        data = {}
        data['pattern_parameters'] = [
            {
                'name': 'amplitudes',
                'values': amplitudes.tolist(),
            }
        ]
        data['configuration'] = {
            'name': name,
            'base_folder': os.path.join(os.path.dirname(__file__), base_folder),
            'time_stamp': time_stamp,
            'useTFT': useTFT,
            'stiffness_pressure': stiffness_pressure,
            'scale_factor_pressure': scale_factor_pressure,
        }
        data['data'] = experiment_logs
        data['duration'] = elapsed_time
        json.dump(data, fp, indent=4, cls = experiment_helper.NpEncoder)


    # visualize the experiment log json file with histogram plots on its fields
