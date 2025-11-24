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
import MeshFEM, mesh, sparse_matrices, benchmark, field_sampler, mesh_utilities
import inflatables_parametrization as parametrization, numpy as np, importlib, pickle, wall_generation
import utils
import py_newton_optimizer
from py_newton_optimizer import NewtonOptimizerOptions
from numpy.linalg import norm
from io_redirection import suppress_stdout
import visualization, wall_width_formulas as wwf
from tri_mesh_viewer import TriMeshViewer

import parallelism
parallelism.set_max_num_tbb_threads(4)
import experiment_pattern_helper
import os 
import inflation
from mesh_utilities import SurfaceSampler, tubeRemesh
import boundaries
import time 
import parametrization_helper
import json 
import sheet_optimizer, opt_config


num_thread = 32

parallelism.set_max_num_tbb_threads(num_thread)
parallelism.set_gradient_assembly_num_threads(num_thread)
parallelism.set_hessian_assembly_num_threads(num_thread)

def run_experiment(shape_index, pattern_index, time_stamp, num_iterations, run_free_boundary = True, fix_feet = False):

    experiment_file, stiffness_path, pattern_name, num_pattern_params, param_index, default_param, param_range, param_normalization_factor, fusing_curve_polyline, shape_name, shape_path, use_holes = experiment_pattern_helper.parse_input(shape_index, pattern_index)

    print("Running fine-tune optimization experiment for shape {} and pattern {}".format(shape_name, pattern_name))

    base_path = os.path.dirname(os.path.abspath(__file__))

    meshing_data_path = '{}/../output/{}/{}_{}/meshing/'.format(base_path, time_stamp, shape_name, pattern_name)
    print(meshing_data_path)
    optimization_data_path = '{}/../output/{}/{}_{}/fine_tune_opt/'.format(base_path, time_stamp, shape_name, pattern_name)
    if not os.path.exists(optimization_data_path):
        os.makedirs(optimization_data_path)

    target_surf = mesh.Mesh(shape_path)
    target_surf.setVertices(utils.prototypeScaleNormalization(target_surf.vertices(), placeAtopFloor=False))
    target_surf = mesh_utilities.subdivide_loop(target_surf, 1)

    m = MeshFEM.mesh.Mesh(meshing_data_path + '/parametrized_mesh.obj')
    vertices = m.vertices()
    m = MeshFEM.mesh.Mesh(np.concatenate((vertices, np.zeros((len(vertices), 1))), axis = 1), m.elements())
    fusing_data = np.load(meshing_data_path + '/fusing_data.npy')

    fusing_lines_info = parametrization_helper.get_fusing_lines_from_mesh(m, fusing_data)
    isheet = inflation.InflatableSheet(m, fusing_data, (fusing_lines_info[3], fusing_lines_info[2]))

    isheet.setUseTensionFieldEnergy(True)
    isheet.setUseHessianProjectedEnergy(False)
    isheet.pressure = 0.025


    uv = np.load(meshing_data_path + '/rparam_uv.npy')
    paramSampler = SurfaceSampler(np.pad(uv, [(0, 0), (0, 1)], 'constant'), target_surf.triangles())
    liftedSheetPositions = paramSampler.sample(m.vertices(), target_surf.vertices())

    opts = py_newton_optimizer.NewtonOptimizerOptions()
    opts.useIdentityMetric = True
    opts.beta = 1e-4
    opts.niter = 2000
    opts.gradTol = 1e-6

    isheet.setUninflatedDeformation(liftedSheetPositions.transpose(), prepareRigidMotionPinConstraints=False)
    # bdryVars = boundaries.getBoundaryVars(isheet)
    bdryVars = boundaries.getOuterBoundaryVars(isheet)

    if shape_name == 'lemonade_stand_modular':
        fix_feet = True

    if fix_feet:
        print("Fix feet!")
        m = isheet.mesh()
        V = m.vertices()
        BV = m.boundaryVertices()
        arclen = lambda l: np.linalg.norm(np.diff(V[BV[np.array(l)]], axis=0), axis=1).sum()
        outerLoopBdryVertices = max(m.boundaryLoops(), key=arclen)
        feet_fixed_vars = []
        for bvi in outerLoopBdryVertices:
            for c in range(3):
                feet_fixed_vars.append(isheet.varIdx(0, BV[bvi], c))
                
        feet_fixed_vars = np.array(feet_fixed_vars).reshape((-1, 3))

        sheet_vars = isheet.getVars()

        fixed_vars_values = sheet_vars[feet_fixed_vars[:, 1]]

        bottom_fixed_vars = []
        for i in range(len(feet_fixed_vars)):
            if np.abs(fixed_vars_values[i] - min(fixed_vars_values))< 1:
                bottom_fixed_vars.extend(feet_fixed_vars[i])
        bdryVars = np.array(bottom_fixed_vars)
        print(bdryVars)

    fixedVars_list = [bdryVars, []]
    tag_name = ['fixed_boundary', 'free_boundary']
    hessianShifts  = [0, 1e-7]

    viewer = TriMeshViewer(isheet, width=768, height=640)
    viewer.showWireframe(True)
    # viewer.setCameraParams(((1.613494603240345, -3.9332708615926393, 1.4922998234349831),
    # (-0.05948468564942635, 0.33267929672385665, 0.941162078339598),
    # (0.0, 0.0, 0.0)))
    viewer.update(scalarField=utils.getStrains(isheet)[:, 0])    

    framerate = 10
    def cb(it):
        if it % framerate == 0:
            viewer.update(scalarField=utils.getStrains(isheet)[:, 0])

    experiment_log = {}
    run_time_stamp = time.strftime("%Y_%m_%d_%H_%M")

    experiment_range = run_free_boundary + 1
    with open(optimization_data_path + 'stdout.txt', 'w') as f:
        original_stdout = sys.stdout
        sys.stdout = f
        for i in range(experiment_range):
        # for i in range(1):
            start_time = time.time()
            fixedVars = fixedVars_list[i]
            tag = tag_name[i]
            hessianShift = hessianShifts[i]
            print("Run tag: {}".format(tag))
            if os.path.exists(optimization_data_path + '{}_{}_{}.pkl.gz'.format(tag, shape_name, pattern_name)):
                print("Load existing optimization for tag {}".format(tag))
                sheet_opt = sheet_optimizer.load(optimization_data_path + '{}_{}_{}.pkl.gz'.format(tag, shape_name, pattern_name))
                targetAttractedSheet = sheet_opt.rso.targetAttractedInflation()
            else:
            
                # Reset the inflation and set up target-attraction forces
                isheet.setUninflatedDeformation(liftedSheetPositions.transpose(), prepareRigidMotionPinConstraints=False)
                if pattern_name == "square_with_ellipse_hole_angle_width_height":
                    wallVtxOnBoundary = boundaries.getWallVtxOnOuterBoundary(isheet)
                else:
                    # Let TargetAttractedInflation compute it.
                    wallVtxOnBoundary = []
                targetAttractedSheet = inflation.TargetAttractedInflation(isheet, target_surf, wallVtxOnBoundary)
                targetAttractedSheet.energy(targetAttractedSheet.EnergyType.Fitting)
                
                targetAttractedSheet.targetSurfaceFitter().holdClosestPointsFixed = True
                targetAttractedSheet.fittingWeight = 1e-5

                # Re-inflate, this time applying target-attraction forces.
                # First start with low pressure.
                benchmark.reset()
                opts.niter = 10
                isheet.pressure = 1e-7
                print("First run low pressure")
                cr = inflation.inflation_newton(targetAttractedSheet, fixedVars, opts, hessianShift = hessianShift, callback = cb)
                benchmark.report()
                benchmark.reset()
                opts.niter = 10
                isheet.pressure = 1e-5
                print("First run low pressure")
                cr = inflation.inflation_newton(targetAttractedSheet, fixedVars, opts, hessianShift = hessianShift, callback = cb)
                benchmark.report()
                benchmark.reset()
                opts.niter = 10
                isheet.pressure = 1e-2
                print("First run low pressure")
                cr = inflation.inflation_newton(targetAttractedSheet, fixedVars, opts, hessianShift = hessianShift, callback = cb)
                benchmark.report()

                benchmark.reset()
                opts.niter = 2000
                isheet.pressure = 0.1
                cr = inflation.inflation_newton(targetAttractedSheet, fixedVars, opts, hessianShift = hessianShift, callback = cb)
                benchmark.report()

                # Set up the sheet optimizer
                origDesignMesh = isheet.mesh().copy()

                def config(so):
                    fcs = so.rso.fusingCurveSmoothness()
                    fcs.interiorWeight = 1/10

                sheet_opt = sheet_optimizer.PySheetOptimizer(targetAttractedSheet, fixedVars, renderMode=sheet_optimizer.RenderMode.OFFSCREEN, screenshotPath=optimization_data_path + '{}_{}_{}_{}.mp4'.format(run_time_stamp, tag, shape_name, pattern_name),
                                                            detActivationThreshold=0.9, detActivationThresholdTubeTri=0.5,
                                                            originalDesignMesh=origDesignMesh, fusingCurveSmoothnessConfig=opt_config.FusingCurveSmoothnessParams(0.0, 0.0, 2.0, 2.0), customConfigCallback=config)
            
                # Configure some more weights
                sheet_opt.rso.compressionPenaltyWeight = 1e-6
                fcs = sheet_opt.rso.fusingCurveSmoothness()
                fcs.interiorWeight = 0.05

            # Tall hill
            sheet_opt.deploy_viewer.setCameraParams(((0.5651460114929754, -4.541740208827363, 1.4314586974945211),
 (-0.023135663753620576, 0.2979010983972417, 0.9543163399189993),
 (0.0, 0.0, 0.0)))
            # sheet_opt.flat_viewer.setCameraParams(((-0.005724654779466134, -0.15742800643532207, 3.663779058858404),
            #     (0.0, 1.0, 0.0),
            #     (-0.005724654779466134, -0.15742800643532207, 0.0)))
            if sheet_opt.flat_viewer is not None:
                sheet_opt.flat_viewer.showWireframe()
            sheet_opt.deploy_viewer.scalarFieldGetter = visualization.ISheetScalarField.TGT_DIST(targetAttractedSheet.sheet(), target_surf)

            sheet_opt.setSolver(sheet_optimizer.Solver.SCIPY, num_iterations)
            sheet_opt.rso.getEquilibriumSolver().options.niter = 20

            fcs = sheet_opt.rso.fusingCurveSmoothness()
            fcs.lengthScaleSmoothingWeight = 2
            fcs.curvatureWeight = 2

            sheet_opt.optimize()
            sheet_opt.save(optimization_data_path + '{}_{}_{}_{}.pkl.gz'.format(run_time_stamp, tag, shape_name, pattern_name))

            channelMargin = 0
            final_vertices, concatenated_polylines = parametrization_helper.get_fabrication_file_from_mesh(sheet_opt.rso.sheet(), sheet_opt.rso.sheet().mesh(), fusing_data, channelMargin, [1e-4, 1e-4, 1e-4], optimization_data_path + '{}_{}_{}_{}_design_optimized_sheet_pattern_margin_{}.obj'.format(run_time_stamp, tag, shape_name, pattern_name, channelMargin))
            final_vertices, concatenated_polylines = parametrization_helper.get_fabrication_file_from_mesh(sheet_opt.rso.sheet(), sheet_opt.rso.sheet().mesh(), fusing_data, channelMargin, [1e-4, 1e-4, 1e-4], optimization_data_path + '{}_{}_{}_{}_design_optimized_sheet_pattern_margin_{}.svg'.format(run_time_stamp, tag, shape_name, pattern_name, channelMargin), use_obj = False)
            parametrization_helper.export_top_bottom_mesh(isheet, optimization_data_path, shape_name, pattern_name)
            # visualization.plot_line_segments(final_vertices, concatenated_polylines, width = 20, height = 20, path = optimization_data_path + '{}_{}_{}_{}_design_optimized_sheet_pattern_margin_{}.png'.format(run_time_stamp, tag, shape_name, pattern_name, channelMargin))

            end_time = time.time()
            experiment_log[tag] = {'time': end_time - start_time, 'energy': sheet_opt.rso.energy()}
        
        sys.stdout = original_stdout

    with open('{}/experiment_result.json'.format(optimization_data_path), 'w') as fp:
        json.dump(experiment_log, fp, indent=4)

if __name__ == '__main__':
    # data_time_stamp = 'meshing_output_low_res_2024_01_18_23_29'
    # output_time_stamp = time.strftime("%Y_%m_%d_%H_%M")
    # output_time_stamp = 'low_res_2024_01_18_23_29'
    time_stamp = 'demo'
    num_iterations = 200

    # run_experiment(1, 1)
    # run_experiment(5, 0)
    # run_experiment(6, 0)
    # run_experiment(7, 1)
    # run_experiment(0, 2)
    run_experiment(1, 1, time_stamp, num_iterations)