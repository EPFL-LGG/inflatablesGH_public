import os
import sys; sys.path.append(os.path.join(os.path.dirname(__file__), '..'));

import experiment_helper
from periodic_simulation_setup import *

import parallelism, multiprocessing, itertools, setproctitle
import os, time, numpy as np

import json

name = "dash_line"
time_stamp = time.strftime("%Y_%m_%d_%H_%M")
# time_stamp = "2024_01_21_21_08"
base_folder = 'output/{}/{}'.format(name, time_stamp)
if not os.path.exists(base_folder):
    os.makedirs(base_folder)  

useTFT = True
stiffness_pressure = 0.4
scale_factor_pressure = 0.3
avg_len = 0.075

def get_label(r, angle):
    return "{:.2f}_{:.2f}".format(r, angle)

def run_mirror_cosine_dash(r, angle):
    result_folder = "{}/{}".format(base_folder, get_label(r, angle))
    if not os.path.exists(result_folder):
        os.makedirs(result_folder)  
    if (os.path.isfile('{}/experiment_result.json'.format(result_folder))):
        data = json.load(open('{}/experiment_result.json'.format(result_folder)))
        if data["Planar equilibrium"] == 1:
            # print("has {:.2f}_{:.2f}".format(r, angle))
            return
    print("running {:.2f}_{:.2f}".format(r, angle))

    h = 5
    dash_point = np.array([np.cos(angle / 180 * np.pi), np.sin(angle / 180 * np.pi)]) * r + np.array([0, 0])
    m, marker = pattern_generator_using_gmsh.get_single_dash(h, avg_len, avg_len, dash_point = dash_point)

    finalMarkers = np.where(np.array(marker) == 1)[0]
    m, finalMarkers = periodic_unit_helper.shift_and_merge_2D_periodic_mesh(m, finalMarkers, flip_orientation= 0)
    m, finalMarkers = periodic_unit_helper.shift_and_merge_2D_periodic_mesh(m, finalMarkers, axis = 1, flip_orientation= 1)

    fusedVtx = get_fusedVtx_using_markers(len(m.vertices()), finalMarkers)

    ipu = inflation.InflatablePeriodicUnit(m, fusedVtx = fusedVtx, epsilon = 1e-9)

    experiment_log = experiment_helper.run_experiment(ipu, m, fusedVtx, stiffness_pressure, scale_factor_pressure, name, "{:.2f}_{:.2f}".format(r, angle), result_folder, useTFT = useTFT, use_low_pressure=False)
    experiment_log['time_stamp'] = time.time()

    with open('{}/experiment_result.json'.format(result_folder), 'w') as fp:
        json.dump(experiment_log, fp, indent=4)


# import cProfile, pstats, io

# def profile_test():
#     run_mirror_cosine_dash(amplitudes[0], radius[0], angles[0])


# profiler = cProfile.Profile()
# profiler.enable()
# profile_test()
# profiler.disable()
# s = io.StringIO()
# stats = pstats.Stats(profiler, stream = s).sort_stats('tottime') # tottime, ncalls
# stats.print_stats()
# with open('cprofile_output.txt', 'w+') as f:
#     f.write(s.getvalue())

if __name__ == '__main__':
    radius = np.linspace(0.5, 2.5, 21)
    angles = np.linspace(45, 90, 16)

    num_thread = 1
    # For embarassingly parallel runs, using single thread should be more efficient
    parallelism.set_max_num_tbb_threads(num_thread)

    # Function arguments
    arg1 = radius
    arg2 = angles

    args = itertools.product(arg1, arg2)

    # Set num logical cores
    logical_core_count = int(multiprocessing.cpu_count() * 0.5)

    print("Using {} logical cores".format(logical_core_count))

    import time
    st = time.time()
    with open('{}/time_stamp.txt'.format(base_folder), 'w') as fp:
        fp.write(str(time.time()))
        
    with multiprocessing.Pool(logical_core_count) as p:
        p.starmap(run_mirror_cosine_dash, args)
            
    # get the end time
    et = time.time()

    # get the execution time
    elapsed_time = et - st
    print('Execution time:', elapsed_time, 'seconds')

    experiment_logs = []
    for r in radius:
        for angle in angles:
            result_folder = "{}/{}".format(base_folder, get_label(r, angle))
            # Load the experiment log json and add it to the total experiment logs
            with open('{}/experiment_result.json'.format(result_folder), 'r') as fp:
                curr_json = json.load(fp)
                curr_json['name'] = get_label(r, angle)
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
