import experiment_helper
from periodic_simulation_setup import *
import json
import os

import parallelism, multiprocessing, itertools, setproctitle
import os, time, numpy as np

allowBending = False
useTFT = True


name = 'scale_zigzag_line_10_to_90'
time_stamp = time.strftime("%Y_%m_%d_%H_%M")
# time_stamp = '2023_12_08_11_46'
base_folder = 'output/{}/{}'.format(name, time_stamp)
if not os.path.exists(base_folder):
    os.makedirs(base_folder)  

# pressure = 0.8
stiffness_pressure = 0.3
scale_factor_pressure = 0.01
scales = np.concatenate([np.flip(1 / np.linspace(1, 5, 17))[:-1], np.linspace(1, 5, 17)], axis = 0)


def run_zigzag_experiment(index):
    label = 33
    result_folder = "{}/{}".format(base_folder, index)
    if not os.path.exists(result_folder):
        os.makedirs(result_folder)  

    if (os.path.isfile('{}/experiment_result.json'.format(result_folder))):
        print("has {}".format(index))
        return
    print("running {}".format(index))

    with open('../data/ZigZag_10_90_deg/ZigZag_10_90_deg{}.json'.format(label), 'r') as f:
        data = json.load(f)
    fusedVertices = data['FusedVertices']
    fusedVertices = [True if vx == 1 else False for vx in fusedVertices]
    V = data['Vertices']
    V = np.array(V) * scales[index]
    F = data['Faces']
    m = MeshFEM.Mesh(V, F)
    ipu = inflation.InflatablePeriodicUnit(m, fusedVertices)

    finalMarkers = np.where(np.array(fusedVertices) == 1)[0]
    m, finalMarkers = periodic_unit_helper.shift_and_merge_2D_periodic_mesh(m, finalMarkers)
    m, finalMarkers = periodic_unit_helper.shift_and_merge_2D_periodic_mesh(m, finalMarkers, axis = 1)

    fusedVtx = get_fusedVtx_using_markers(len(m.vertices()), finalMarkers)
    ipu = inflation.InflatablePeriodicUnit(m, fusedVtx = fusedVtx, epsilon = 1e-9)

    experiment_log = experiment_helper.run_experiment(ipu, m, fusedVtx, stiffness_pressure, scale_factor_pressure, name, index, result_folder, useTFT = useTFT, use_low_pressure=False, allow_bending=False)

    with open('{}/experiment_result.json'.format(result_folder), 'w') as fp:
        json.dump(experiment_log, fp, indent=4)



if __name__ == '__main__':
    labels = np.arange(len(scales))

    num_thread = 1
    # For embarassingly parallel runs, using single thread should be more efficient
    parallelism.set_max_num_tbb_threads(num_thread)

    # Set num logical cores
    logical_core_count = int(multiprocessing.cpu_count() * 0.4)

    print("Using {} logical cores".format(logical_core_count))

    import time
    st = time.time()

    with multiprocessing.Pool(logical_core_count) as p:
        p.map(run_zigzag_experiment, iter(labels))
            
    # get the end time
    et = time.time()

    # get the execution time
    elapsed_time = et - st
    print('Execution time:', elapsed_time, 'seconds')

    experiment_logs = []
    for index in labels:
        result_folder = "{}/{}".format(base_folder, index)
        # Load the experiment log json and add it to the total experiment logs
        with open('{}/experiment_result.json'.format(result_folder), 'r') as fp:
            curr_json = json.load(fp)
            curr_json['name'] = index
            experiment_logs.append(curr_json)


    import json
    with open('{}/experiment_result.json'.format(base_folder), 'w') as fp:
        data = {}

        data['configuration'] = {
            'name': name,
            'base_folder': os.path.join(os.path.dirname(__file__), base_folder),
            'time_stamp': time_stamp,
            'useTFT': useTFT,
            'stiffness_pressure': stiffness_pressure,
            'scale_factor_pressure': scale_factor_pressure,
        }
        data['data'] = experiment_logs
        json.dump(data, fp, indent=4, cls = experiment_helper.NpEncoder)


    # visualize the experiment log json file with histogram plots on its fields
