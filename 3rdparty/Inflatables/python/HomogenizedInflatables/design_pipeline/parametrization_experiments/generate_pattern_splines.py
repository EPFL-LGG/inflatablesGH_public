#!/usr/bin/env python
# coding: utf-8


import numpy as np
import sys
sys.path.append("../../")
sys.path.append("../../Visualization/")
sys.path.append("../../../")
sys.path.append(".")
sys.path.append('../../periodic_patches/'); sys.path.append('../'); sys.path.append('../../../gmsh')
import inflation, sparse_matrices, mesh, numpy as np, importlib, pickle
import inflatables_parametrization as parametrization
from numpy.linalg import norm
from io_redirection import suppress_stdout
import visualization
import visualize_stiffness, importlib

sys.path.append('periodic_patches/')
sys.path.append('gmsh')


import parametrization_experiment_helper

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
parallelism.set_max_num_tbb_threads(32)
parallelism.set_gradient_assembly_num_threads(32)
parallelism.set_hessian_assembly_num_threads(32)
import utils, mesh_utilities

import parametrization_optimization_helper

time_stamp = time.strftime("%Y_%m_%d_%H_%M")
print(time_stamp)
import os

base_folder = 'output/{}/'.format(time_stamp)
if not os.path.exists(base_folder):
    os.makedirs(base_folder)  

import serialization_helper

def run_experiment(shape_index, pattern_index):
    print("running {} {}".format(shape_index, pattern_index))
    experiment_log = {}

    start_time = time.time()

    pattern = parametrization_experiment_helper.Pattern_data[pattern_index]

    experiment_file = pattern['experiment_file']
    stiffness_path = pattern['stiffness_path']
    pattern_name = pattern['name']
    num_pattern_params = pattern['num_pattern_params']
    param_index = pattern['param_index']
    default_param = pattern['default_param']
    param_range = pattern['param_range']
    param_normalization_factor = pattern['param_normalization_factor']
    fusing_curve_polyline = pattern['fusing_curve_polyline_function']
                                    
    shape = parametrization_experiment_helper.Shape_data[shape_index]

    shape_name = shape['name']
    shape_path = shape['path']

    print(shape_name, pattern_name)
    data_path = "{}/{}_{}".format(base_folder, shape_name, pattern_name)
    if not os.path.exists(data_path):
        os.makedirs(data_path)  

    target_surf = mesh.Mesh(shape_path)
    target_surf.setVertices(utils.prototypeScaleNormalization(target_surf.vertices(), placeAtopFloor=False))
    target_surf = mesh_utilities.subdivide_loop(target_surf, 1)
    target_surf.save(data_path + '/target_surf.obj')

    # ### Overview

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

    # ### Parametrization

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

    # lg.runIteration()

    # local_global_end_time = time.time()
    # experiment_log['local_global_time'] = local_global_end_time - local_global_start_time

    # visualization.visualize_both(lg, show_main = True, path = data_path + '/local_global_parametrization.png')

    # parametrization_helper.visualize_scale_factors(eqns, lg.getAlphas(), lg.getBetas(), path = data_path + '/scale_factors_local_global.png')

    # # ### Pattern parameters optimization
    # default_pattern_params = []
    # for i in range(len(parameters)):
    #     default_pattern_params += [default_param[i]]  * len(lg.getAlphas())
    # rparam = parametrization.RegularizedPatternParametrizer(lg, splines, default_pattern_params, len(grid_data.shape) - 1)
    # rparam.patternParamBounds = np.array(param_range)
    # rparam.diffRegW = 0.0
    # visualization.visualize_both(rparam, height = 4)

    # serialization_helper.save_parametrization_classes(target_surf, parametrization.lscm(target_surf), lg, splines, default_pattern_params, num_pattern_params, rparam, data_path + '/local_global_initialized_parametrization_classes.pkl.gz')

    # print("Initialize pattern parameters")
    # parametrization_start_time = time.time()
    # experiment_log['parametrization_time'] = parametrization_end_time - parametrization_start_time
    # init_energy_values = parametrization_optimization_helper.initialize_pattern_parameters(rparam, lines, num_pattern_params, path = data_path + '/pattern_parameters_initialization.png')
    # serialization_helper.save_parametrization_classes(target_surf, parametrization.lscm(target_surf), lg, splines, default_pattern_params, num_pattern_params, rparam, data_path + '/pattern_initialized_parametrization_classes.pkl.gz')

    # print("Add phi regularization")
    # phiRegW, after_phi_energy_values = parametrization_optimization_helper.add_phi_regularization(rparam, lines, num_pattern_params, path = data_path + '/phi_regularization.png')
    # serialization_helper.save_parametrization_classes(target_surf, parametrization.lscm(target_surf), lg, splines, default_pattern_params, num_pattern_params, rparam, data_path + '/phi_regularized_parametrization_classes.pkl.gz')

    # print("Add pattern regularization")
    # patternRegW, after_patternReg_energy_values = parametrization_optimization_helper.add_pattern_regularization(rparam, lines, num_pattern_params, phiRegW, path = data_path + '/pattern_regularization.png')
    # serialization_helper.save_parametrization_classes(target_surf, parametrization.lscm(target_surf), lg, splines, default_pattern_params, num_pattern_params, rparam, data_path + '/pattern_regularized_parametrization_classes.pkl.gz')

    # # ### Bending
    # print("Add bending energy")
    # bendRegW, after_bending_energy_values = parametrization_optimization_helper.add_bending_energy(rparam, lines, num_pattern_params, phiRegW, patternRegW, path = data_path + '/bending')
    # serialization_helper.save_parametrization_classes(target_surf, parametrization.lscm(target_surf), lg, splines, default_pattern_params, num_pattern_params, rparam, data_path + '/bending_parametrization_classes.pkl.gz')

    # parametrization_end_time = time.time()

    # experiment_log['init_energy_values'] = init_energy_values
    # experiment_log['after_phi_energy_values'] = after_phi_energy_values
    # experiment_log['after_patternReg_energy_values'] = after_patternReg_energy_values
    # experiment_log['after_bending_energy_values'] = after_bending_energy_values

    # # ## Upsampling and channel generation

    # pattern_generation_start_time = time.time()
    # fusing_lines = fusing_curve_polyline(default_param)[0].reshape(-1, 2)
    # fusing_edges = [[i, i + 1] for i in range(len(fusing_lines) - 1)]

    # boundary_vertices = [[0, 0], [np.pi, 0], [np.pi, np.pi], [0, np.pi]]

    # boundary_edges = np.array([[0, 1], [1, 2], [2, 3], [3, 0]]) + len(fusing_lines)

    # visualization.plot_line_segments(list(fusing_lines) + boundary_vertices, fusing_edges + list(boundary_edges))

    # default_frequency = 0.2
    # default_mesh_size = 4
    # default_edge_soup_threshold = 1e0

    # scale = 1

    # frequency = default_frequency * scale
    # # mesh_size = default_mesh_size / scale
    # edge_soup_threshold = default_edge_soup_threshold

    # print("Generate pattern")
    # sdfVertices, sdfTris, sdf, sheet_vxs, concatenated_polylines, sheet_edges_polylines,  boundaryVxs, boundaryEdges, upsampleMesh_vertices, upsampleMesh_triangles, upsampledAngles, upsampledPatternParams = parametrization_helper.get_polyline_from_pattern_parameters(rparam, fusing_curve_polyline, nsubdiv = 4, frequency=frequency, duplicates_removable_threshold=[1e-4, 1e-2, 1e-1, edge_soup_threshold, edge_soup_threshold * 2], path = data_path)

    # if len(boundaryEdges) > 1:
    #     concatenated_boundary_edges = []
    #     for polyline in boundaryEdges:
    #         concatenated_boundary_edges.extend(polyline)
    #     concatenated_boundary_edges = np.array(concatenated_boundary_edges)
    #     # Define a function to calculate the length of a sublist
    #     def sublist_length(sublist):
    #         return len(sublist)

    #     # Sort boundaryEdges in descending order of sublist length
    #     boundaryEdges = sorted(boundaryEdges, key=sublist_length, reverse=True)
    # else:
    #     concatenated_boundary_edges = np.array(boundaryEdges[0])

    # pattern_generation_end_time = time.time()
    # experiment_log['pattern_generation_time'] = pattern_generation_end_time - pattern_generation_start_time

    # visualization.scalarFieldPlotFast(sdfVertices, sdfTris, sdf, width = 5, height=5, path = data_path + '/sdf.png')
    # visualization.plot_line_segments(list(sheet_vxs) + list(boundaryVxs), list(concatenated_polylines) + list(concatenated_boundary_edges + len(sheet_vxs)), width = 5, height = 5, path = data_path + '/fusing_curves.png')
    # # plt.scatter(boundaryVxs[concatenated_boundary_edges[:,0], 0], boundaryVxs[concatenated_boundary_edges[:,0], 1], c = np.arange(len(concatenated_boundary_edges[:, 0])), cmap = mpl.colormaps['Greys'])

    # experiment_log['time'] = time.time() - start_time
    # with open('{}/experiment_result.json'.format(data_path), 'w') as fp:
    #     json.dump(experiment_log, fp, indent=4)

import parallelism, multiprocessing, itertools, setproctitle

if __name__ == '__main__':
    pattern_indices = np.arange(len(parametrization_experiment_helper.Pattern_data))[:3]