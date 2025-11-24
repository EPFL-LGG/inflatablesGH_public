#!/usr/bin/env python
# coding: utf-8


import numpy as np
import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

sys.path.append("../../")
sys.path.append("../../Visualization/")
sys.path.append("../../../")
sys.path.append("..")
sys.path.append('../../periodic_patches/'); sys.path.append('../'); sys.path.append('../../../gmsh')
import inflation, sparse_matrices, mesh, numpy as np, importlib, pickle
import inflatables_parametrization as parametrization
from numpy.linalg import norm
from io_redirection import suppress_stdout
import visualization
import visualize_stiffness, importlib

sys.path.append('periodic_patches/')
sys.path.append('gmsh')


import experiment_pattern_helper

import pandas as pd
import json
import matplotlib.pyplot as plt


# ### Get scale function convex hull
import matplotlib.cm as cm
import matplotlib as mpl

from scipy.spatial import ConvexHull, convex_hull_plot_2d
import parametrization_helper

import parametrization_helper
import visualization
import time    
import os


from periodic_simulation_setup import *

from py_newton_optimizer import NewtonOptimizerOptions

import MeshFEM, parallelism, benchmark, utils
import utils, mesh_utilities

import parametrization_optimization_helper

import serialization_helper

def get_pattern_info(experiment_file, stiffness_path, pattern_name, param_index, data_path):
    with open(experiment_file, 'r') as fp:
        data = json.load(fp)

    df = pd.DataFrame(data['data'])
    valid_tags = np.array(df['name'][df['Planar equilibrium'] == 1])

    kappa_path = None

    bending_stiffness_data, stretching_stiffness_data, scale_factor_data, used_tags = visualize_stiffness.plot_all_data(kappa_path, stiffness_path, pattern_name, valid_tags, plot_data = False)

    parameters = []
    for index in param_index:
        parameters.append(data['pattern_parameters'][index]['values'])

    max_bending_stiffness = np.max(bending_stiffness_data, axis = 1)
    min_bending_stiffness = np.min(bending_stiffness_data, axis = 1)
    max_stretching_stiffness = np.max(stretching_stiffness_data, axis = 1)
    min_stretching_stiffness = np.min(stretching_stiffness_data, axis = 1)

    x_scale_factors, y_scale_factors = visualize_stiffness.get_axis_scale_factors(stiffness_path, pattern_name, valid_tags)

    min_scale_factors = np.min(np.concatenate((x_scale_factors.reshape(-1, 1), y_scale_factors.reshape(-1, 1)), axis = 1), axis = 1)
    max_scale_factors = np.max(np.concatenate((x_scale_factors.reshape(-1, 1), y_scale_factors.reshape(-1, 1)), axis = 1), axis = 1)


    points = np.concatenate((max_scale_factors.reshape((-1, 1)), min_scale_factors.reshape((-1, 1))), axis = 1)
    hull = ConvexHull(points)

    # ### Validate the max and min scale factors are aligned with the x and y axis
    eqns = hull.equations

    parametrization_helper.visualize_scale_factors(eqns, max_scale_factors, min_scale_factors, path = data_path + '/scale_factors.png')

    # ### Generate data without augmenting
    stiffness_coefficients = np.array(visualize_stiffness.get_stiffness_coefficients(stiffness_path, pattern_name, (used_tags)))
    # For patches with reflection symmetry:
    stiffness_coefficients[:, 1] *= 0
    stiffness_coefficients[:, 2] *= 0

    grid_shape = np.zeros(len(parameters) + 1, dtype = np.int64)
    grid_shape[0] = 9
    for i in range(len(parameters)):
        grid_shape[i + 1] = len(parameters[i])

    grid_data = np.zeros(tuple(grid_shape))

    stiffness_coefficients.reshape([*grid_shape[1:], 5]).shape

    grid_stiffness_coefficients = stiffness_coefficients.reshape([*grid_shape[1:], 5])
    grid_stiffness_coefficients = np.transpose(grid_stiffness_coefficients, (len(grid_stiffness_coefficients.shape)-1,) + tuple(range(len(grid_stiffness_coefficients.shape)-1)))

    grid_data[0] = max_scale_factors.reshape(grid_shape[1:])
    grid_data[1] = min_scale_factors.reshape(grid_shape[1:])
    grid_data[2] = x_scale_factors.reshape(grid_shape[1:])
    grid_data[3] = y_scale_factors.reshape(grid_shape[1:])
    grid_stiffness_coefficients = stiffness_coefficients.reshape([*grid_shape[1:], 5])
    grid_stiffness_coefficients = np.transpose(grid_stiffness_coefficients, (len(grid_stiffness_coefficients.shape)-1,) + tuple(range(len(grid_stiffness_coefficients.shape)-1)))
    for s in range(5):
        grid_data[4 + s] = grid_stiffness_coefficients[s]


    if len(parameters) == 1:
        splines = parametrization_helper.scipy_get_mat_params_over_one_pattern_params_grid_interpolation(parameters, grid_data)
    elif len(parameters) == 2:
        splines = parametrization_helper.scipy_get_mat_params_over_pattern_params_grid_interpolation(parameters[0], parameters[1], grid_data)
    else:
        splines = parametrization_helper.ndsplines_get_mat_params_over_pattern_params_grid_interpolation(grid_data, *parameters)

    return splines, hull, eqns



def run_experiment(shape_index, pattern_index, time_stamp = time.strftime("%Y_%m_%d_%H_%M"), rerun_experiment = False, use_knitro = True):
    # time_stamp = '2024_01_18_23_29'

    print(time_stamp)
    base_folder = '{}/../output/{}/'.format(os.path.dirname(os.path.abspath(__file__)), time_stamp)
    if not os.path.exists(base_folder):
        os.makedirs(base_folder)  

    experiment_log = {}

    start_time = time.time()

    experiment_file, stiffness_path, pattern_name, num_pattern_params, param_index, default_param, param_range, param_normalization_factor, fusing_curve_polyline, shape_name, shape_path, use_holes = experiment_pattern_helper.parse_input(shape_index, pattern_index)

    data_path = "{}/{}_{}/parametrization/".format(base_folder, shape_name, pattern_name)
    if not os.path.exists(data_path):
        os.makedirs(data_path)  

    if (not rerun_experiment) and os.path.exists('{}/experiment_result.json'.format(data_path)):
        print("Experiment already ran for {} {}".format(shape_index, pattern_index))
        return
    
    print("running {} {}".format(shape_index, pattern_index))

    splines, hull, eqns = get_pattern_info(experiment_file, stiffness_path, pattern_name, param_index, data_path)
    # ### Parametrization
    print(shape_path)
    target_surf = mesh.Mesh(shape_path)
    target_surf.setVertices(utils.prototypeScaleNormalization(target_surf.vertices(), placeAtopFloor=False))
    target_surf = mesh_utilities.subdivide_loop(target_surf, 1)
    target_surf.save(data_path + '/target_surf.obj')

    lines = np.array(eqns)
    # ### New local global with convex hull
    local_global_start_time = time.time()
    lg = parametrization.LocalGlobalGenericParametrizer(target_surf, parametrization.lscm(target_surf))
    lg.setLines(eqns)
    lg.alphaMin = hull.min_bound[0]
    lg.alphaMax = hull.max_bound[0]
    lg.betaMin = hull.min_bound[1]
    lg.betaMax = hull.max_bound[1]
    for i in range(1000): lg.runIteration()
    lg.runIteration()
    local_global_end_time = time.time()
    experiment_log['local_global_time'] = local_global_end_time - local_global_start_time
    visualization.visualize_both(lg, show_main = True, path = data_path + '/local_global_parametrization.png')
    parametrization_helper.visualize_scale_factors(eqns, lg.getAlphas(), lg.getBetas(), path = data_path + '/scale_factors_local_global.png')

    # ### Pattern parameters optimization
    default_pattern_params = []
    for i in range(num_pattern_params):
        default_pattern_params += [default_param[i]]  * len(lg.getAlphas())
    rparam = parametrization.RegularizedPatternParametrizer(lg, splines, default_pattern_params, num_pattern_params)
    rparam.patternParamBounds = np.array(param_range)
    rparam.diffRegW = 0.0
    visualization.visualize_both(rparam, height = 4)
    serialization_helper.save_parametrization_classes(target_surf, parametrization.lscm(target_surf), lg, splines, default_pattern_params, num_pattern_params, rparam, data_path + '/local_global_initialized_parametrization_classes.pkl.gz')
    
    parametrization_start_time = time.time()
    init_energy_values = parametrization_optimization_helper.initialize_pattern_parameters(rparam, lines, num_pattern_params, path = data_path + '/pattern_parameters_initialization.png', use_knitro = use_knitro)
    serialization_helper.save_parametrization_classes(target_surf, parametrization.lscm(target_surf), lg, splines, default_pattern_params, num_pattern_params, rparam, data_path + '/pattern_initialized_parametrization_classes.pkl.gz')

    phiRegW, after_phi_energy_values = parametrization_optimization_helper.add_phi_regularization(rparam, lines, num_pattern_params, path = data_path + '/phi_regularization.png', use_knitro = use_knitro)
    serialization_helper.save_parametrization_classes(target_surf, parametrization.lscm(target_surf), lg, splines, default_pattern_params, num_pattern_params, rparam, data_path + '/phi_regularized_parametrization_classes.pkl.gz')

    patternRegW, after_patternReg_energy_values = parametrization_optimization_helper.add_pattern_regularization(rparam, lines, num_pattern_params, phiRegW, path = data_path + '/pattern_regularization.png', use_knitro = use_knitro)
    serialization_helper.save_parametrization_classes(target_surf, parametrization.lscm(target_surf), lg, splines, default_pattern_params, num_pattern_params, rparam, data_path + '/pattern_regularized_parametrization_classes.pkl.gz')

    # ### Bending
    bendRegW, after_bending_energy_values = parametrization_optimization_helper.add_bending_energy(rparam, lines, num_pattern_params, phiRegW, patternRegW, path = data_path + '/bending', use_knitro = use_knitro)
    serialization_helper.save_parametrization_classes(target_surf, parametrization.lscm(target_surf), lg, splines, default_pattern_params, num_pattern_params, rparam, data_path + '/bending_parametrization_classes.pkl.gz')
    experiment_log['init_energy_values'] = init_energy_values
    experiment_log['after_phi_energy_values'] = after_phi_energy_values
    experiment_log['after_patternReg_energy_values'] = after_patternReg_energy_values
    experiment_log['after_bending_energy_values'] = after_bending_energy_values

    parametrization_end_time = time.time()
    experiment_log['parametrization_time'] = parametrization_end_time - parametrization_start_time

    experiment_log['time'] = time.time() - start_time
    with open('{}/experiment_result.json'.format(data_path), 'w') as fp:
        json.dump(experiment_log, fp, indent=4)

import parallelism, multiprocessing, itertools, setproctitle

if __name__ == '__main__':
    # shape_indices = np.arange(len(experiment_pattern_helper.Shape_data))
    # pattern_indices = np.arange(len(experiment_pattern_helper.Pattern_data))[:3]
    # # Function arguments
    # arg1 = shape_indices
    # arg2 = pattern_indices

    # args = itertools.product(arg1, arg2)

    # Crafted arguments
    # taller hill with holes, igloo with all, vest with cosine, neck brace with holes, cashew with cosine, squiward with dashline, two rings with dash lines
    # args = [(0, 2), (1, 0), (1, 1), (1, 2), (2, 1), (3, 2), (4, 1), (5, 0), (6, 0)]
    args = [(1, 1)]


    num_thread = 4

    parallelism.set_max_num_tbb_threads(num_thread)
    parallelism.set_gradient_assembly_num_threads(num_thread)
    parallelism.set_hessian_assembly_num_threads(num_thread)


    for arg in args:
        run_experiment(*arg)

    # # Set num logical cores
    # logical_core_count = int(multiprocessing.cpu_count() * 0.7)

    # print("Using {} logical cores".format(logical_core_count))

    # import time
    # st = time.time()

    # with multiprocessing.Pool(logical_core_count) as p:
    #     p.starmap(run_experiment, args)
            
    # # get the end time
    # et = time.time()

    # # get the execution time
    # elapsed_time = et - st
    # print('Execution time:', elapsed_time, 'seconds')

    # experiment_logs = []
    # # for shape_index in shape_indices:
    # #     for pattern_index in pattern_indices:
    # for shape_index, pattern_index in args:
    #         pattern = experiment_pattern_helper.Pattern_data[pattern_index]
    #         shape = experiment_pattern_helper.Shape_data[shape_index]
    #         shape_name = shape['name']
    #         pattern_name = pattern['name']

    #         data_path = "{}/{}_{}".format(base_folder, shape_name, pattern_name)
    #         # Load the experiment log json and add it to the total experiment logs
    #         with open('{}/experiment_result.json'.format(data_path), 'r') as fp:
    #             curr_json = json.load(fp)
    #             curr_json['name'] = "{}_{}".format(shape_name, pattern_name)
    #             experiment_logs.append(curr_json)


    # import json
    # with open('{}/experiment_result.json'.format(base_folder), 'w') as fp:
    #     data = {}
    #     data['data'] = experiment_logs
    #     data['elapsed_time'] = elapsed_time
    #     json.dump(data, fp, indent=4)

