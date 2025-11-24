import os
import sys; sys.path.append(os.path.join(os.path.dirname(__file__), '..'));

import experiment_helper
from periodic_simulation_setup import *

import parallelism, multiprocessing, itertools, setproctitle
import os, time, numpy as np

import json

num_points = 5

name = "voronoi_{}".format(num_points)
time_stamp = time.strftime("%Y_%m_%d_%H_%M")
# time_stamp = "2024_01_07_22_56"
base_folder = 'output/{}/{}'.format(name, time_stamp)
if not os.path.exists(base_folder):
    os.makedirs(base_folder)  

useTFT = True
stiffness_pressure = 0.4
scale_factor_pressure = 0.01
avg_len = 0.1

num_experiment = 200
all_points = np.random.uniform(low=-2.5, high=2.5, size=(num_points * num_experiment, 2))


def run_voronoi(i):
    result_folder = "{}/{}".format(base_folder, i)
    if not os.path.exists(result_folder):
        os.makedirs(result_folder)  
    if (os.path.isfile('{}/experiment_result.json'.format(result_folder))):
        data = json.load(open('{}/experiment_result.json'.format(result_folder)))
        # print("has {:.2f}_{:.2f}_{:.2f}".format(amp, r, angle))
        return
    print("running {}".format(i))
    points  = all_points[i * num_points : (i + 1) * num_points]

    m, marker, vertices, edges = pattern_generator_using_gmsh.generate_voronoi_mesh(points, avg_len, avg_len, shrink_percentage=0.66)
    np.save("{}/vertices_{}_{}.npy".format(result_folder, name, i), vertices)
    np.save("{}/edges_{}_{}.npy".format(result_folder, name, i), edges)
    
    finalMarkers = np.where(np.array(marker) == 1)[0]
    m, finalMarkers = periodic_unit_helper.shift_and_merge_2D_periodic_mesh(m, finalMarkers, flip_orientation= 0)
    m, finalMarkers = periodic_unit_helper.shift_and_merge_2D_periodic_mesh(m, finalMarkers, axis = 1, flip_orientation= 1)

    fusedVtx = get_fusedVtx_using_markers(len(m.vertices()), finalMarkers)

    ipu = inflation.InflatablePeriodicUnit(m, fusedVtx = fusedVtx, epsilon = 1e-9)

    experiment_log = experiment_helper.run_experiment(ipu, m, fusedVtx, stiffness_pressure, scale_factor_pressure, name, i, result_folder, useTFT = useTFT, use_low_pressure=False)

    with open('{}/experiment_result.json'.format(result_folder), 'w') as fp:
        json.dump(experiment_log, fp, indent=4)

if __name__ == '__main__':
    num_thread = 1
    # For embarassingly parallel runs, using single thread should be more efficient
    parallelism.set_max_num_tbb_threads(num_thread)

    # Set num logical cores
    logical_core_count = int(multiprocessing.cpu_count() * 0.7)

    print("Using {} logical cores".format(logical_core_count))

    import time
    st = time.time()

    with multiprocessing.Pool(logical_core_count) as p:
        p.map(run_voronoi, iter(range(num_experiment)))

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
                'name': 'num_points',
                'values': [num_points],
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
        data['elapsed_time'] = elapsed_time
        json.dump(data, fp, indent=4)


    # visualize the experiment log json file with histogram plots on its fields