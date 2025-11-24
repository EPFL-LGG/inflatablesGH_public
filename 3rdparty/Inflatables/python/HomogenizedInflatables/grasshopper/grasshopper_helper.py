import sys; sys.path.append('../../'); sys.path.append('../../periodic_patches/'); sys.path.append('../experiments/'); sys.path.append('../../gmsh')
import inflation, sparse_matrices, mesh, numpy as np, importlib, pickle
import inflatables_parametrization as parametrization
from numpy.linalg import norm
from io_redirection import suppress_stdout
import visualization
import json
import os
import os.path as osp

sys.path.append('periodic_patches/')
sys.path.append('gmsh')

from periodic_simulation_setup import *

from py_newton_optimizer import NewtonOptimizerOptions

import MeshFEM, parallelism, benchmark, utils
parallelism.set_max_num_tbb_threads(32)
parallelism.set_gradient_assembly_num_threads(32)
parallelism.set_hessian_assembly_num_threads(32)

import parametrization_helper
importlib.reload(parametrization_helper)

import scipy
from scipy.interpolate import RectBivariateSpline

import sheet_meshing, inflation
import py_newton_optimizer
import time

gh_dir = './'

def generate_target_mesh(model_name):

    ################################################
    ################################################
    # Init Target Surface
    filename = osp.join(gh_dir + 'inputs/{}.json'.format(model_name))

    # Import converted data
    with open(filename) as json_file:
        data = json.load(json_file)

    temp_vertices = data['TargetSurface']['Vertices']
    num_vertices = len(temp_vertices)
    vertices = np.ndarray(shape=(num_vertices, 3))
    for i in range(num_vertices):
        v = temp_vertices[i]
        vertices[i] = [float(v[0]), float(v[1]), float(v[2])]

    temp_faces = data['TargetSurface']['Trias']
    num_faces = len(temp_faces)
    faces = np.ndarray(shape=(num_faces, 3))
    for i in range(num_faces):
        f = temp_faces[i]
        faces[i] = [int(f[0]), int(f[1]), int(f[2])]
    
    input_mesh = mesh.Mesh(vertices, faces)
    num_iterations = int(data['NumIterations'])
    ################################################
    ################################################

    # TODO:
    lines = np.array([[-1.        , -1.        ,  2.415     ],
       [-2.30769231, -1.        ,  3.81923077],
       [ 0.77304965,  1.        , -2.27304965],
       [ 1.29357798,  1.        , -2.94036697],
       [-0.43333333, -1.        ,  1.655     ]])
    
    lg = parametrization.LocalGlobalGenericParametrizer(input_mesh, parametrization.lscm(input_mesh))

    for i in range(1000): lg.runIteration()

    # lg.alphaMin = 1.354641650129586
    # lg.alphaMax = 1.3936682623611498

    # lg.betaMin = 1.1282764273066557
    # lg.betaMax = 1.141976581731514

    # lg.alphaMin = 1.2349791815955107 #1.354641650129586
    # lg.alphaMax = 1.2827966225656238 #1.3936682623611498

    # lg.betaMin = 1.2349791815955107 #1.1282764273066557
    # lg.betaMax = 1.2827966225656238 #1.141976581731514


    lg.setLines(lines)
    lg.alphaMin = 1.0
    lg.alphaMax = 1.5

    lg.betaMin = 1.0
    lg.betaMax = 1.5


    (1.3936682623611498, 1.354641650129586, 1.141976581731514, 1.1282764273066557)

    #print(lg.energy())
    lg.runIteration()
    #print(lg.energy())

    for i in range(5000): lg.runIteration()
    #print(lg.energy())

    ################################################
    ################################################
    # Get splines
    grid_data = np.load("../Visualization/grid_data.npy")
    grid_pattern_1 = np.load("../Visualization/grid_pattern_1.npy")
    grid_pattern_2 = np.load("../Visualization/grid_pattern_2.npy")

    ################################################
    ################################################
    objective_splines = [RectBivariateSpline(grid_pattern_1, grid_pattern_2, grid_data[i], kx = 4, ky = 4) for i in range(7)]
    splines = parametrization_helper.get_mat_params_over_pattern_params_grid_interpolation(grid_pattern_1, grid_pattern_2, grid_data)
    default_pattern_params = [0.5 * (np.max(grid_pattern_1) + np.min(grid_pattern_1))] * len(lg.getAlphas()) + [0.5 * (np.max(grid_pattern_2) + np.min(grid_pattern_2))] * len(lg.getAlphas())
    rparam = parametrization.RegularizedPatternParametrizer(lg, splines, default_pattern_params)
    rparam.patternParamBounds = np.array([[0.5, 2.4], [0, 90]])

    PET = parametrization.RegularizedPatternParametrizer.PatternEnergyType
    list(map(rparam.energy, [PET.Full, PET.RGP, PET.Bending]))

    ################################################
    ################################################
    def optimize_rparam(param, patternRegW, phiRegW, bendRegW = 0.0, update_uv = True, niter = 100):
        param.patternRegW = patternRegW
        param.phiRegW = phiRegW
        param.bendRegW = bendRegW
        param.diffRegW = 0.0
        opts = NewtonOptimizerOptions()
        opts.useIdentityMetric = True
        opts.beta = 1e-4
        opts.niter = niter
        opts.gradTol = 1e-9
        opts.factorizer = opts.factorizer.CatamariNesdis
        benchmark.reset()
        
        if update_uv:
            fixedvars = [param.uOffset(), param.vOffset(), param.phiOffset()]
        else:
            fixedvars = range(param.stretchOffset())

        cr = parametrization.pattern_parametrization_knitro(param, opts.niter, fixedvars)
        benchmark.report()
        return cr
    
    ################################################
    ################################################

    rparam.bendRegW = 1e-4

    ################################################
    ################################################
    PET = parametrization.RegularizedPatternParametrizer.PatternEnergyType
    list(map(rparam.energy, [PET.Full, PET.PatternRegularization, PET.Bending, PET.DEBUG_Fitting, PET.DEBUG_StretchRegularization, PET.DEBUG_PhiRegularization]))

    ################################################
    ################################################
    benchmark.reset()
    with suppress_stdout(): report = optimize_rparam(rparam, 0, 0, bendRegW = 0, update_uv = False, niter = 200)
    benchmark.report()

    ################################################
    ################################################
    PET = parametrization.RegularizedPatternParametrizer.PatternEnergyType
    list(map(rparam.energy, [PET.Full, PET.PatternRegularization, PET.Bending, PET.DEBUG_Fitting, PET.DEBUG_StretchRegularization, PET.DEBUG_PhiRegularization]))

    ################################################
    ################################################
    benchmark.reset()
    with suppress_stdout(): report = optimize_rparam(rparam, patternRegW = 0, phiRegW = 0, bendRegW = 1e-2, update_uv = True, niter = 100)
    benchmark.report()

    ################################################
    ################################################
    PET = parametrization.RegularizedPatternParametrizer.PatternEnergyType
    list(map(rparam.energy, [PET.Full, PET.PatternRegularization, PET.Bending, PET.DEBUG_Fitting, PET.DEBUG_StretchRegularization, PET.DEBUG_PhiRegularization]))

    ################################################
    ################################################
    # Upsampling and channel generation
    nsubdiv=8
    upsampledMesh, upsampledAngles, upsampledPatternParams = rparam.upsampledVertexLeftStretchAnglesAndPatternParameters(nsubdiv)

    ################################################
    ################################################
    radius_data =  upsampledPatternParams[0]
    radius_data = radius_data / 2.5 * (np.pi / 2)

    angles_data =  upsampledPatternParams[1]
    angles_data = angles_data / 180 * np.pi

    (sdfVertices, sdfTris, sdf) = wall_generation.evaluate_cross_field(upsampledMesh.vertices(), upsampledMesh.triangles(), upsampledAngles, radius_data, angles_data, frequency=60, margin = 0.07)


    ################################################
    ################################################
    # plot 
    pts, edges = wall_generation.extract_contours(sdfVertices, sdfTris, sdf,
                                              targetEdgeSpacing=0.002,
                                              minContourLen=0.02)

    ################################################
    ################################################
    # Meshing and inflation simulation

    m, iwv, iwbv = sheet_meshing.newMeshingAlgorithm(sdfVertices, sdfTris, sdf, pts, edges, triArea=0.0001)
    #visualization.plot_2d_mesh(m, pointList=np.where(np.array(iwv) == 1)[0], width=10, height=10)

    isheet = inflation.InflatableSheet(m, np.array(iwv) != 0)

    from mesh_utilities import SurfaceSampler, tubeRemesh


    paramSampler = SurfaceSampler(np.pad(rparam.uv(), [(0, 0), (0, 1)], 'constant'), input_mesh.triangles())
    liftedSheetPositions = paramSampler.sample(m.vertices(), input_mesh.vertices())

    isheet.setUninflatedDeformation(liftedSheetPositions.transpose())
    isheet.getVars()

    ################################################
    ################################################
    iterations_per_output = 10
    opts = py_newton_optimizer.NewtonOptimizerOptions()
    opts.useIdentityMetric = True
    opts.beta = 1e-4
    opts.gradTol = 1e-10
    opts.niter = iterations_per_output

    isheet.setUseTensionFieldEnergy(True)
    isheet.setUseHessianProjectedEnergy(False)

    fixedVars, hessianShift = [], 1e-6

    isheet.pressure = 1e-5  
    opts.niter = 5

    cr = inflation.inflation_newton(isheet, fixedVars, opts, hessianShift = hessianShift)

    # Then inflate
    isheet.pressure = 1e-1

    opts.niter = num_iterations

    benchmark.reset()
    cr = inflation.inflation_newton(isheet, fixedVars, opts, hessianShift = hessianShift)
    benchmark.report()

    out_filename = osp.join(gh_dir + 'outputs/{}_optimized.json'.format(model_name))
    write_optimization_info_json(isheet, out_filename)

    print('Optimization completed!')

def write_optimization_info_json(isheet : inflation.InflatableSheet, filename):
    
    data = isheet.visualizationGeometry()
    v_list = [ [float(vertices[0]), float(vertices[1]), float(vertices[2])] for vertices in data[0] ]
    f_list = [ [int(face[0]), int(face[1]), int(face[2])] for face in data[1] ]

    model_info = {'vertices' : v_list, 'faces' : f_list }
    
    with open(filename, 'w') as f:
        json.dump(model_info, f, indent=4)
