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
os.makedirs('benchmark_results', exist_ok=True)

name, uvPath, optPath = sys.argv[1:]

opt = utils.load(optPath)
targetAttractedSheet = copy.deepcopy(opt.rso.targetAttractedInflation())
oldOptISheet = copy.deepcopy(targetAttractedSheet.sheet())
target_surf = utils.getTargetSurf(targetAttractedSheet)
origMesh = opt.rso.originalMesh().copy()
optDoF = opt.rso.numVars(),
del opt

import json
model_settings = json.load(open('model_settings.json', 'r'))

################################################################################
# Run Default Parametrization Settings
################################################################################

# Choose reasonable stretching bounds
alphaMin = wwf.stretchFactorForCanonicalWallWidth(wwf.canonicalWallWidthForGeometry(2, 10))
alphaMax = wwf.stretchFactorForCanonicalWallWidth(wwf.canonicalWallWidthForGeometry(1, 10))
print(alphaMin, alphaMax)

lg = parametrization.LocalGlobalParametrizer(target_surf, parametrization.lscm(target_surf))

benchmark.start_timer_section('param')

for i in range(1000): lg.runIteration()
lg.alphaMin = alphaMin
lg.alphaMax = alphaMax
for i in range(5000): lg.runIteration()

rparam = parametrization.RegularizedParametrizerSVD(target_surf, lg.uv())
rparam.alphaMin = alphaMin
rparam.alphaMax = alphaMax

def optimize_rparam(param, alphaRegW, phiRegW, bendRegW):
    param.alphaRegW = alphaRegW
    param.phiRegW = phiRegW
    param.bendRegW = bendRegW
    opts = NewtonOptimizerOptions()
    opts.niter = 1000
    opts.hessianProjectionController = py_newton_optimizer.HessianProjectionAdaptive()
    #opts.hessianProjectionController = py_newton_optimizer.HessianProjectionNever()
    cr = parametrization.regularized_parametrization_newton(param, param.rigidMotionPinVars, opts)

# Default optmiization weight sequence
with suppress_stdout(): optimize_rparam(rparam, 1e2, 1e1, 500.0)
with suppress_stdout(): optimize_rparam(rparam, 1e1, 1e0, 250.0)
with suppress_stdout(): optimize_rparam(rparam, 1e0, 1e-1, 125.0)
with suppress_stdout(): optimize_rparam(rparam, 1e-1, 1e-2, 62.5)
with suppress_stdout(): optimize_rparam(rparam, 1e-2, 1e-3, 31.25)

benchmark.stop_timer_section('param')

################################################################################
# Stripe Pattern Generation and Meshing with Default Settings
################################################################################
benchmark.start_timer_section('stripe')
nsubdiv=3
upsampledMesh, upsampledAngles, upsampledStretches = rparam.upsampledVertexLeftStretchAnglesAndMagnitudes(nsubdiv)
upsampledStretches = np.clip(upsampledStretches, alphaMin, alphaMax)
(sdfVertices, sdfTris, sdf) = wall_generation.evaluate_stripe_field(upsampledMesh.vertices(), upsampledMesh.triangles(), upsampledAngles,
                                                                    wwf.canonicalWallWidthForStretchFactor(upsampledStretches), frequency=0.2)
benchmark.stop_timer_section('stripe')

import sheet_meshing
triArea = 6.0
benchmark.start_timer_section('meshing')
pts, edges = wall_generation.extract_contours(sdfVertices, sdfTris, sdf,
                                              targetEdgeSpacing=4.0,
                                              minContourLen=10)
m, iwv, iwbv = sheet_meshing.newMeshingAlgorithm(sdfVertices, sdfTris, sdf, pts, edges, triArea=triArea)
benchmark.stop_timer_section('meshing')

################################################################################
# Initial Equilibrium Solve (Fixed Boundary)
################################################################################
import copy
uv = utils.load(uvPath)

isheet = targetAttractedSheet.sheet()
m = isheet.mesh()
#isheet.mesh().save('opt_design.msh')
isheet.mesh().setVertices(origMesh.vertices()) # Revert to the original design sheet stored in the optimizer
#isheet.mesh().save('orig_design.msh')

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

print('running equilibrium')

benchmark.start_timer_section('equilibrium')

cr = inflation.inflation_newton(targetAttractedSheet, fixedVars, opts)

benchmark.stop_timer_section('equilibrium')

print('done')

################################################################################
# Sheet Optimization with Default Settings
################################################################################
import boundaries
fixedVars = boundaries.getOuterBoundaryVars(isheet)

import sheet_optimizer, opt_config
origDesignMesh = isheet.mesh().copy()

def config(so):
    fcs = so.rso.fusingCurveSmoothness()
    fcs.interiorWeight             = getModelSetting('interiorWeight',             0.1)
    fcs.lengthScaleSmoothingWeight = getModelSetting('lengthScaleSmoothingWeight', 1.0)
    fcs.curvatureWeight            = getModelSetting('curvatureWeight',            1.0)
    print(f'Smoothness regularization settings: {fcs.interiorWeight} {fcs.lengthScaleSmoothingWeight} {fcs.curvatureWeight} ') 

print('running sheet optimization')
new_sheet_opt = sheet_optimizer.PySheetOptimizer(targetAttractedSheet, fixedVars, renderMode=sheet_optimizer.RenderMode.NONE,
                                             detActivationThreshold=0.9, detActivationThresholdTubeTri=0.5,
                                             originalDesignMesh=origDesignMesh, fusingCurveSmoothnessConfig=opt_config.FusingCurveSmoothnessParams(0.0, 0.0, 1.0, 1.0),
                                             customConfigCallback=config)
#new_sheet_opt.setSolver(sheet_optimizer.Solver.SCIPY, maxIters=2) # early termination for debugging this script
benchmark.start_timer_section('optimization')
new_sheet_opt.optimize(doBenchmark=False)
benchmark.stop_timer_section('optimization')

# Save out results of benchmarked optimization
new_sheet_opt.save(f'benchmark_results/opt_{name}.pkl.gz')

################################################################################
# Model Statistics
################################################################################
import measurements

# Measure distances using the original optimization result
dists = measurements.getWallDistances(oldOptISheet, target_surf, relative=True)

report = {
    'sim_dof':    targetAttractedSheet.numVars(),
    'opt_dof':    optDoF,
    'fit_max':    dists.max(),
    'fit_98_pct': np.percentile(dists, 98),
    'fit_dists':  dists,
    'benchmarks': benchmark.to_dict()
}

utils.save(report, f'benchmark_results/{name}.pkl.gz')
