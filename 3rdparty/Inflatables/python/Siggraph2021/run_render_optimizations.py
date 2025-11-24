import sys; sys.path.append('..')
import MeshFEM, mesh, sparse_matrices, benchmark, field_sampler, mesh_utilities
import inflatables_parametrization as parametrization, numpy as np, importlib, pickle, wall_generation
import inflation, utils
import py_newton_optimizer
from py_newton_optimizer import NewtonOptimizerOptions
from numpy.linalg import norm
from io_redirection import suppress_stdout
import visualization, wall_width_formulas as wwf
from tri_mesh_viewer import TriMeshViewer
import copy

import os
os.makedirs('optimization_renders', exist_ok=True)

name, uvPath, optPath = sys.argv[1:]
outPrefix = f'optimization_renders/{name}'

opt = utils.load(optPath)
targetAttractedSheet = copy.deepcopy(opt.rso.targetAttractedInflation())
isheet = targetAttractedSheet.sheet()
oldOptISheet = copy.deepcopy(targetAttractedSheet.sheet())
target_surf = utils.getTargetSurf(targetAttractedSheet)
origMesh = opt.rso.originalMesh().copy()
optDoF = opt.rso.numVars(),
del opt

import json
model_settings = json.load(open('model_settings.json', 'r')).get(name, {})
print('model_settings: ', model_settings)

################################################################################
# Initial Equilibrium Solve (Fixed Boundary)
################################################################################
import copy
uv = utils.load(uvPath)

isheet = targetAttractedSheet.sheet()
m = isheet.mesh()
isheet.mesh().setVertices(origMesh.vertices()) # Revert to the original design sheet stored in the optimizer

paramSampler = field_sampler.FieldSampler(np.pad(uv, [(0, 0), (0, 1)], 'constant'), target_surf.triangles())
liftedSheetPositions = paramSampler.sample(m.vertices(), target_surf.vertices())

isheet.setUninflatedDeformation(liftedSheetPositions.transpose(), prepareRigidMotionPinConstraints=False)
targetAttractedSheet = inflation.TargetAttractedInflation(isheet, target_surf)

def getModelSetting(setting, default):
    return model_settings.get(name, {}).get(setting, default)

fittingWeight = getModelSetting('fittingWeight', 0.0)
targetAttractedSheet.fittingWeight = fittingWeight
print('Using fitting weight: ', fittingWeight)

targetAttractedSheet.targetSurfaceFitter().holdClosestPointsFixed = fittingWeight != 0.0

import py_newton_optimizer
opts = py_newton_optimizer.NewtonOptimizerOptions()
opts.useIdentityMetric = True
opts.beta = 1e-4
opts.gradTol = 1e-7
opts.niter = 5000

import boundaries
fixedVars = boundaries.getOuterBoundaryVars(isheet)

benchmark.start_timer_section('equilibrium')

cr = inflation.inflation_newton(targetAttractedSheet, fixedVars, opts)

################################################################################
# Sheet Optimization with Default Settings
################################################################################
import boundaries
fixedVars = boundaries.getOuterBoundaryVars(isheet)

import sheet_optimizer, opt_config
origDesignMesh = isheet.mesh().copy()

def config(so):
    fcs = so.rso.fusingCurveSmoothness()
    fcs.interiorWeight             = model_settings.get('interiorWeight',             0.1)
    fcs.boundaryWeight             = model_settings.get('boundaryWeight',             1.0)
    fcs.lengthScaleSmoothingWeight = model_settings.get('lengthScaleSmoothingWeight', 1.0)
    fcs.curvatureWeight            = model_settings.get('curvatureWeight',            1.0)
    print(f'Smoothness regularization settings: {fcs.interiorWeight} {fcs.lengthScaleSmoothingWeight} {fcs.curvatureWeight} ') 

print('running sheet optimization')
new_sheet_opt = sheet_optimizer.PySheetOptimizer(targetAttractedSheet, fixedVars, renderMode=sheet_optimizer.RenderMode.OFFSCREEN, screenshotPath=f'{outPrefix}.mp4',
                                             detActivationThreshold=0.9, detActivationThresholdTubeTri=0.5,
                                             originalDesignMesh=origDesignMesh, fusingCurveSmoothnessConfig=opt_config.FusingCurveSmoothnessParams(0.0, 0.0, 1.0, 1.0),
                                             customConfigCallback=config)

new_sheet_opt.flat_viewer.showWireframe(True)
if 'cameraParams' in model_settings:
    new_sheet_opt.deploy_viewer.setCameraParams(model_settings['cameraParams'])
if 'flatCam' in model_settings:
    new_sheet_opt.flat_viewer.setCameraParams(model_settings['flatCam'])

new_sheet_opt.setSolver(sheet_optimizer.Solver.SCIPY, maxIters=200)
new_sheet_opt.optimize()

# Save out results of benchmarked optimization
new_sheet_opt.save(f'{outPrefix}.opt.pkl.gz')
