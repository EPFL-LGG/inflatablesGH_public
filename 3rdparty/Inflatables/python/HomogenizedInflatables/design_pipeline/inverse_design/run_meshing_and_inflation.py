#!/usr/bin/env python
# coding: utf-8
import sys
sys.path.append("../")
sys.path.append("../Visualization/")
sys.path.append("../../")

sys.path.append("../../design_pipeline/parametrization_experiments/")

sys.path.append('../../'); sys.path.append('../../../'); sys.path.append('../../../gmsh'); sys.path.append('../../Visualization/')


sys.path.append('periodic_patches/')
sys.path.append('gmsh')

import numpy as np
import inflation, sparse_matrices, mesh, numpy as np, pickle
import inflatables_parametrization as parametrization
from numpy.linalg import norm
from io_redirection import suppress_stdout
import visualization
import visualize_stiffness
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl
import parametrization_helper
import experiment_pattern_helper
import utils, mesh_utilities, benchmark
import pickle, mesh, wall_generation, visualization, numpy as np
import mesher_helper
import MeshFEM
import numpy as np
import copy
import shapely
from mesh_utilities import SurfaceSampler, tubeRemesh
import py_newton_optimizer
# from tri_mesh_viewer import TriMeshViewer
from visualization import TriMeshViewerWithSurface
import boundaries
from matplotlib import pyplot as plt
import sheet_optimizer, opt_config
import fabrication
import py_newton_optimizer
import boundaries
import time
import os
import parallelism, multiprocessing, itertools, setproctitle
import serialization_helper

import json


def generate_fusing_lines(fusing_curve_polyline, default_param, input_data_path, output_data_path, frequency, edge_soup_threshold):

    if os.path.exists(input_data_path + '/bending_parametrization_classes.pkl.gz'):
        rparam, lg, target_surf, splines, default_pattern_params, num_params = serialization_helper.load_parametrization_classes(input_data_path + '/bending_parametrization_classes.pkl.gz')

    # if os.path.exists(input_data_path + '/pattern_initialized_parametrization_classes.pkl.gz'):
        # rparam, lg, target_surf, splines, default_pattern_params, num_params = serialization_helper.load_parametrization_classes(input_data_path + '/pattern_initialized_parametrization_classes.pkl.gz')
    else:
        print("Parametrization data does not exist for " + input_data_path)
        return 
    
    
    pattern_generation_start_time = time.time()
    fusing_lines = fusing_curve_polyline(default_param)[0].reshape(-1, 2)
    fusing_edges = [[i, i + 1] for i in range(len(fusing_lines) - 1)]

    boundary_vertices = [[0, 0], [np.pi, 0], [np.pi, np.pi], [0, np.pi]]

    boundary_edges = np.array([[0, 1], [1, 2], [2, 3], [3, 0]]) + len(fusing_lines)

    visualization.plot_line_segments(list(fusing_lines) + boundary_vertices, fusing_edges + list(boundary_edges))

    sdfVertices, sdfTris, sdf, sheet_vxs, concatenated_polylines, sheet_edges_polylines,  boundaryVxs, boundaryEdges, upsampleMesh_vertices, upsampleMesh_triangles, upsampledAngles, upsampledPatternParams = parametrization_helper.get_polyline_from_pattern_parameters(rparam, fusing_curve_polyline, nsubdiv = 4, frequency=frequency, duplicates_removable_threshold=[1e-4, 1e-2, 1e-1, edge_soup_threshold, edge_soup_threshold + 1], path = output_data_path)

    if len(boundaryEdges) > 1:
        concatenated_boundary_edges = []
        for polyline in boundaryEdges:
            concatenated_boundary_edges.extend(polyline)
        concatenated_boundary_edges = np.array(concatenated_boundary_edges)
        # Define a function to calculate the length of a sublist
        def sublist_length(sublist):
            return len(sublist)

        # Sort boundaryEdges in descending order of sublist length
        boundaryEdges = sorted(boundaryEdges, key=sublist_length, reverse=True)
    else:
        concatenated_boundary_edges = np.array(boundaryEdges[0])

    visualization.scalarFieldPlotFast(sdfVertices, sdfTris, sdf, width = 5, height=5, path = output_data_path + '/sdf.png')
    visualization.plot_line_segments(list(sheet_vxs) + list(boundaryVxs), list(concatenated_polylines) + list(concatenated_boundary_edges + len(sheet_vxs)), width = 5, height = 5, path = output_data_path + '/fusing_curves.png')
    visualization.plot_line_segments(list(sheet_vxs) + list(boundaryVxs), list(concatenated_polylines) + list(concatenated_boundary_edges + len(sheet_vxs)), width = 5, height = 5, path = output_data_path + '/fusing_curves.svg')
    # Save to obj as well
    parametrization_helper.save_line_segments_as_obj(list(sheet_vxs) + list(boundaryVxs), list(concatenated_polylines) + list(concatenated_boundary_edges + len(sheet_vxs)), path = output_data_path + '/fusing_curves.obj')
    # plt.scatter(boundaryVxs[concatenated_boundary_edges[:,0], 0], boundaryVxs[concatenated_boundary_edges[:,0], 1], c = np.arange(len(concatenated_boundary_edges[:, 0])), cmap = mpl.colormaps['Greys'])
    return sdfVertices, sdfTris, sdf, sheet_vxs, concatenated_polylines, sheet_edges_polylines,  boundaryVxs, boundaryEdges, upsampleMesh_vertices, upsampleMesh_triangles, upsampledAngles, upsampledPatternParams

def run_experiment(shape_index, pattern_index, time_stamp = '2024_01_18_23_29', run_inflation = True, frequency = 0.2, rerun_experiment = False):
    # data_time_stamp = '2024_01_18_00_03'
    # data_time_stamp = '2024_01_18_23_29'
    output_time_stamp = time_stamp
    start_time = time.time()
    experiment_log = {}

    experiment_file, stiffness_path, pattern_name, num_pattern_params, param_index, default_param, param_range, param_normalization_factor, fusing_curve_polyline, shape_name, shape_path, use_holes = experiment_pattern_helper.parse_input(shape_index, pattern_index)

    # Parameters:
    # For igloo, neck brace, hill, squidward
    default_frequency = 0.2
    default_mesh_size = 2

    # # For lemonade, cashew
    # default_frequency = 0.3
    # default_mesh_size = 2

    default_edge_soup_threshold = 1.5e0 

    scale = (frequency / 0.2)
    mesh_size = default_mesh_size / scale
    edge_soup_threshold = default_edge_soup_threshold / scale

    subdiv_scale = int(np.log2(scale))
    print("subdiv_scale", subdiv_scale)
    ##########################################

    base_path = os.path.dirname(os.path.abspath(__file__))

    input_data_path = '{}/../output/{}/{}_{}/parametrization/'.format(base_path, time_stamp, shape_name, pattern_name)
    output_data_path = '{}/../output/{}/{}_{}/meshing/'.format(base_path, output_time_stamp, shape_name, pattern_name)
    if not os.path.exists(output_data_path):
        os.makedirs(output_data_path)


    if (not rerun_experiment) and os.path.exists(output_data_path + '/free_boundary_inflated_sheet_vars.npy'):
        print("Inflation data exists for shape {} and pattern {}".format(shape_name, pattern_name))
        return
    
    if (not rerun_experiment) and os.path.exists(output_data_path + '/sdf.npy'):
        sdfVertices, sdfTris, sdf, sheet_vxs, concatenated_polylines, sheet_edges_polylines,  boundaryVxs, boundaryEdges, upsampleMesh_vertices, upsampleMesh_triangles, upsampledAngles, upsampledPatternParams = parametrization_helper.get_polyline_from_pattern_parameters(None, fusing_curve_polyline, nsubdiv = 4 + subdiv_scale, frequency=frequency, duplicates_removable_threshold=[1e-4, 1e-2, 1e-1, edge_soup_threshold, edge_soup_threshold + 1], path = output_data_path, load_data = "final_results")
    else:
        try:
            sdfVertices, sdfTris, sdf, sheet_vxs, concatenated_polylines, sheet_edges_polylines,  boundaryVxs, boundaryEdges, upsampleMesh_vertices, upsampleMesh_triangles, upsampledAngles, upsampledPatternParams = generate_fusing_lines(fusing_curve_polyline, default_param, input_data_path, output_data_path, frequency, edge_soup_threshold)
        except:
            print("Fusing line generation failed for shape {} and pattern {}".format(shape_name, pattern_name))
            return

    print("size of concatenated polylines: ", len(concatenated_polylines))
    # if (shape_name == 'lilium' and pattern_index == 2):
    #     return

    target_surf = mesh.Mesh(shape_path)
    target_surf.setVertices(utils.prototypeScaleNormalization(target_surf.vertices(), placeAtopFloor=False))
    target_surf = mesh_utilities.subdivide_loop(target_surf, 1)


    if len(boundaryEdges) > 1:
        concatenated_boundary_edges = []
        for polyline in boundaryEdges:
            concatenated_boundary_edges.extend(polyline)
        concatenated_boundary_edges = np.array(concatenated_boundary_edges)
        # Define a function to calculate the length of a sublist
        def sublist_length(sublist):
            return len(sublist)

        # Sort boundaryEdges in descending order of sublist length
        boundaryEdges = sorted(boundaryEdges, key=sublist_length, reverse=True)
    else:
        concatenated_boundary_edges = np.array(boundaryEdges[0])

    fusing_curve_time = time.time()

    # ## Meshing and inflation simulation

    if (not rerun_experiment) and os.path.exists(output_data_path + '/parametrized_mesh.obj') and os.path.exists(output_data_path + '/fusing_data.npy'):
        print("Meshing data exists for shape {} and pattern {}".format(shape_name, pattern_name))
        m = MeshFEM.mesh.Mesh(output_data_path + '/parametrized_mesh.obj')
        vertices = m.vertices()
        m = MeshFEM.mesh.Mesh(np.concatenate((vertices, np.zeros((len(vertices), 1))), axis = 1), m.elements())
        fusing_data = np.load(output_data_path + '/fusing_data.npy')

    else:
        # ## Meshing and inflation simulation

        if use_holes:
            boundary_holes_vxs, non_boundary_holes_vxs, boundary_polygon, removed_holes_vxs, hole_marker = parametrization_helper.post_process_holes(boundaryVxs, boundaryEdges, sheet_vxs, sheet_edges_polylines, smoothing = 5.0, area_threshold=mesh_size **2, avg_len = mesh_size, distance_threshod=mesh_size)
        else:
            selected_elements = [np.array(sublist)[:, 0] for sublist in boundaryEdges[1:]]
            non_boundary_holes_vxs = list(boundaryVxs[selected_elements])

        boundary_curve = boundaryVxs[np.array(boundaryEdges[0])[:, 0]]

        # Check if the boundary curve is in counter-clockwise order
        if not parametrization_helper.is_counter_clockwise(boundary_curve):
            # If not, reverse it
            boundary_curve = boundary_curve[::-1]
        boundary_curve = parametrization_helper.smooth_polyline(boundary_curve, 0, default_mesh_size)[:-1]

        for i in range(len(non_boundary_holes_vxs)):
            if not parametrization_helper.is_counter_clockwise(non_boundary_holes_vxs[i]):
                # The interior holes should be in clockwise order
                non_boundary_holes_vxs[i] = non_boundary_holes_vxs[i][::-1]
            non_boundary_holes_vxs[i] = parametrization_helper.smooth_polyline(non_boundary_holes_vxs[i], 0, default_mesh_size)[:-1]
            if len(non_boundary_holes_vxs[i]) > 4:
                non_boundary_holes_vxs[i] = parametrization_helper.smooth_polyline(non_boundary_holes_vxs[i], 0, default_mesh_size)[:-1]

        if use_holes:
            fusing_lines, non_boundary_holes_vxs, boundary_polygon = parametrization_helper.ellipsify_holes(boundaryVxs, boundaryEdges, sheet_vxs, sheet_edges_polylines, smoothing = 5.0, area_threshold=mesh_size **2, avg_len = mesh_size, distance_threshod=mesh_size * 4)
            v, f, fusing_data = mesher_helper.generate_mesh_non_periodic(default_mesh_size, boundary_curve, [], non_boundary_holes_vxs, [], [], gui = False)
        else:
            v, f, fusing_data = mesher_helper.generate_mesh_non_periodic(default_mesh_size, boundary_curve, [], non_boundary_holes_vxs, sheet_vxs, concatenated_polylines, gui = False)


        # Use the function
        new_v, new_f, new_fusing_without_boundary = parametrization_helper.remove_dangling_vertices(v, f - 1, fusing_data)
        m = MeshFEM.mesh.Mesh(new_v, new_f)
        new_fusing = copy.copy(new_fusing_without_boundary)
        new_fusing[m.boundaryVertices()] = True

        m.save(output_data_path + '/parametrized_mesh.obj')
        np.save(output_data_path + '/fusing_data.npy', new_fusing)
        fusing_data = new_fusing

        visualization.plot_2d_mesh(m, pointList=np.where(np.array(new_fusing) == 1)[0], width=10, height=10, path=output_data_path + '/parametrized_mesh.png')

    meshing_time = time.time()

    if run_inflation:
        isheet = inflation.InflatableSheet(m, fusing_data)
        uv = np.load(output_data_path + '/rparam_uv.npy')

        paramSampler = SurfaceSampler(np.pad(uv, [(0, 0), (0, 1)], 'constant'), target_surf.triangles())
        liftedSheetPositions = paramSampler.sample(m.vertices(), target_surf.vertices())
        isheet.setUninflatedDeformation(liftedSheetPositions.transpose())

        opts = py_newton_optimizer.NewtonOptimizerOptions()
        opts.useIdentityMetric = True
        opts.beta = 1e-4
        opts.gradTol = 1e-10

        sys.path.append("../../")

        bdryVars = boundaries.getOuterBoundaryVars(isheet)
        isheet.setUseTensionFieldEnergy(True)
        isheet.setUseHessianProjectedEnergy(False)
        # First generate results with fixed boundary 
        viewer = TriMeshViewerWithSurface(isheet, target_surf, width=768, height=640)
        viewer.showWireframe(True)

        # viewer.setCameraParams(((1.613494603240345, -3.9332708615926393, 1.4922998234349831),
        # (-0.05948468564942635, 0.33267929672385665, 0.941162078339598),
        # (0.0, 0.0, 0.0)))

        viewer.update(scalarField=utils.getStrains(isheet)[:, 0])    

        framerate = 10
        def cb(it):
            if it % framerate == 0:
                viewer.update(scalarField=utils.getStrains(isheet)[:, 0])    


        # ### First solve with low pressure to get out of indefinite state
        isheet.pressure = 1e-7
        opts.niter = 20
        cr = inflation.inflation_newton(isheet, bdryVars, opts, hessianShift = 0, callback = cb)

        isheet.pressure = 1e-3
        opts.niter = 20
        cr = inflation.inflation_newton(isheet, bdryVars, opts, hessianShift = 0, callback = cb)

        isheet.pressure = 1e-2
        opts.niter = 20
        cr = inflation.inflation_newton(isheet, bdryVars, opts, hessianShift = 0, callback = cb)

        fixedVars_list = [bdryVars, []]
        tag_name = ['fixed_boundary', 'free_boundary']
        hessian_shifts = [0, 1e-6]


        for i in range(len(fixedVars_list)):
            fixedVars = fixedVars_list[i]
            tag = tag_name[i]
            hessian_shift = hessian_shifts[i]
            if os.path.exists(output_data_path + '/{}_inflated_sheet_vars.npy'.format(tag)):
                isheet.setVars(np.load(output_data_path + '/{}_inflated_sheet_vars.npy'.format(tag)))
                viewer.update(scalarField=utils.getStrains(isheet)[:, 0])    

            else:
                # ### Then inflate
                isheet.pressure = 0.025
                opts.niter = 100
                opts.gradTol = 1e-7

                benchmark.reset()
                cr = inflation.inflation_newton(isheet, fixedVars, opts, hessianShift = hessian_shift, callback = cb)
                benchmark.report()

                np.save(output_data_path + '/{}_inflated_sheet_vars.npy'.format(tag), isheet.getVars())

            orender = viewer.offscreenRenderer(width=1024,height=1024)
            orender.render()
            orender.save(output_data_path + '/{}_parametrized_mesh_inflated.png'.format(tag))

            def export_top_bottom_mesh(isheet, export_path, shape_name, pattern_name):
                mesh_3d = isheet.visualizationMesh(True)
                mesh_2d = isheet.mesh()
                vx_3d = mesh_3d.vertices()
                elements_3d = mesh_3d.elements()

                new_mesh_3d = MeshFEM.Mesh(vx_3d[:mesh_2d.numVertices()], elements_3d[:mesh_2d.numElements()])

                new_mesh_3d.save(export_path + '/{}_{}_{}_mesh_3d_top.obj'.format(tag, shape_name, pattern_name))

                new_mesh_3d = MeshFEM.Mesh(vx_3d[m.numVertices():], elements_3d[mesh_2d.numElements():] - mesh_2d.numVertices())
                new_mesh_3d.save(export_path + '/{}_{}_{}_mesh_3d_bottom.obj'.format(tag, shape_name, pattern_name))

            export_top_bottom_mesh(isheet, output_data_path, shape_name, pattern_name)
    inflation_time = time.time()

    experiment_log['fusing_curve_time'] = fusing_curve_time - start_time
    experiment_log['meshing_time'] = meshing_time - fusing_curve_time
    experiment_log['inflation_time'] = inflation_time - meshing_time
    with open('{}/experiment_result.json'.format(output_data_path), 'w') as fp:
        json.dump(experiment_log, fp, indent=4)

if __name__ == '__main__':
    num_thread = 10

    parallelism.set_max_num_tbb_threads(num_thread)
    parallelism.set_gradient_assembly_num_threads(num_thread)
    parallelism.set_hessian_assembly_num_threads(num_thread)

    args = [(1, 1)]

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




