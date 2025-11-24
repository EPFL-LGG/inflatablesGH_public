#!/usr/bin/env python
# coding: utf-8


import numpy as np
import sys
sys.path.append("../../")
sys.path.append("../../Visualization/")



sys.path.append("../../../")



import sys; sys.path.append('..')
import MeshFEM, mesh, sparse_matrices, benchmark, field_sampler, mesh_utilities
import inflatables_parametrization as parametrization, numpy as np, importlib, pickle, wall_generation
import utils
import py_newton_optimizer
from py_newton_optimizer import NewtonOptimizerOptions
from numpy.linalg import norm
from io_redirection import suppress_stdout
import visualization, wall_width_formulas as wwf

import time
time_stamp = time.strftime("%Y_%m_%d_%H_%M")

import os

# name = "cashew_planar_coarse"
name = "igloo"

optimization_data_path = 'output/parallel_tube_optimization/{}/{}/'.format(time_stamp, name)
if not os.path.exists(optimization_data_path):
    os.makedirs(optimization_data_path)

print("running optimization for {} at {}".format(name, optimization_data_path))

target_surf = mesh.Mesh("../../../../examples/{}.obj".format(name))
target_surf.setVertices(utils.prototypeScaleNormalization(target_surf.vertices(), placeAtopFloor=True))
target_surf = mesh_utilities.subdivide_loop(target_surf, 1)



# Choose reasonable stretching bounds in terms of the relative fusing curve widths.
alphaMin = wwf.stretchFactorForCanonicalWallWidth(wwf.canonicalWallWidthForGeometry(2, 10))
alphaMax = wwf.stretchFactorForCanonicalWallWidth(wwf.canonicalWallWidthForGeometry(1, 10))
print(alphaMin, alphaMax)



# Run some iterations of the local-global algorithm to ensure a good separation between singular values.
# This step can also be used as a prediction of the feasiblity of a design surface:
# if it is unable to nearly satisfy the singular value constraints,
# the surface is probably infeasible.
lg = parametrization.LocalGlobalParametrizer(target_surf, parametrization.lscm(target_surf))

lg.alphaMin = 1.4
lg.alphaMax = np.pi / 2
print(lg.energy())
for i in range(1000): lg.runIteration()

print(lg.energy())
lg.runIteration()
print(lg.energy())



visualization.visualize(lg)



rparam = parametrization.RegularizedParametrizerSVD(target_surf, lg.uv())
rparam.alphaMin = alphaMin
rparam.alphaMax = alphaMax



def optimize_rparam(param, alphaRegW, phiRegW, bendRegW):
    param.alphaRegW = alphaRegW
    param.phiRegW = phiRegW
    param.bendRegW = bendRegW
    opts = NewtonOptimizerOptions()
    opts.niter = 200
    opts.hessianProjectionController = py_newton_optimizer.HessianProjectionAdaptive()
    #opts.hessianProjectionController = py_newton_optimizer.HessianProjectionNever()
    cr = parametrization.regularized_parametrization_newton(param, param.rigidMotionPinVars, opts)



# Rerunning this cell a couple times can improve the results
benchmark.reset()
with suppress_stdout(): optimize_rparam(rparam, 100.0, 10.0, 500.0)
with suppress_stdout(): optimize_rparam(rparam, 10.0, 1.0, 250.0)
with suppress_stdout(): optimize_rparam(rparam, 1.0, 0.1, 125.0)
with suppress_stdout(): optimize_rparam(rparam, 0.1, 0.01, 62.5)
with suppress_stdout(): optimize_rparam(rparam, 0.1, 0.01, 31.25)
benchmark.report()



importlib.reload(utils)



# Report the values and gradients of each objective term
print(f'Energies: {utils.allEnergies(rparam)}')
print(f'Gradient Norms: {utils.allGradientNorms(rparam)}')



# Visualize the flattening
visualization.visualize(rparam)



importlib.reload(visualization)



visualization.visualizeChannelOrientation(rparam, quiver=visualization.QuiverVisualization.PER_VTX, orientationHue=False, width = 10, height = 10)


# ## Upsampling and channel generation


import time



nsubdiv=4
start_time = time.time()
upsampledMesh, upsampledAngles, upsampledStretches = rparam.upsampledVertexLeftStretchAnglesAndMagnitudes(nsubdiv)
upsampledStretches = np.clip(upsampledStretches, alphaMin, alphaMax)
(sdfVertices, sdfTris, sdf) = wall_generation.evaluate_stripe_field(upsampledMesh.vertices(), upsampledMesh.triangles(), upsampledAngles,
                                                                    wwf.canonicalWallWidthForStretchFactor(upsampledStretches), frequency=0.7)
print("computing sdf takes: ", time.time() - start_time)



import pickle, mesh, wall_generation, visualization, numpy as np



visualization.scalarFieldPlotFast(sdfVertices, sdfTris, sdf, height=12)



pts, edges = wall_generation.extract_contours(sdfVertices, sdfTris, sdf,
                                              targetEdgeSpacing=4.0,
                                              minContourLen=10)



visualization.plot_line_segments(pts, edges, width=10, height=16)


# ## Meshing and inflation simulation


import sheet_meshing, inflation



import importlib
importlib.reload(sheet_meshing)



m, iwv, iwbv = sheet_meshing.newMeshingAlgorithm(sdfVertices, sdfTris, sdf, pts, edges, triArea=2)



isheet = inflation.InflatableSheet(m, iwv)
isheet.setRelaxedStiffnessEpsilon(1e-6)
uv = rparam.uv()



# Are the flat region causing a problem? They might not actually control the metric...
# Try replacing them with single wall...
# Analyze the actual stretching factor (much easier to do with skeleton walls)

from tri_mesh_viewer import TriMeshViewer
viewer = TriMeshViewer(isheet, width=768, height=640)
viewer.showWireframe(True)

viewer.show()



# Manually stretch the sheet onto the target surface by applying the inverse of the parametrization
paramSampler = field_sampler.FieldSampler(np.pad(uv, [(0, 0), (0, 1)], 'constant'), target_surf.triangles())
liftedSheetPositions = paramSampler.sample(m.vertices(), target_surf.vertices())

isheet.setUninflatedDeformation(liftedSheetPositions.transpose(), prepareRigidMotionPinConstraints=False)



import py_newton_optimizer
iterations_per_output = 10
opts = py_newton_optimizer.NewtonOptimizerOptions()
opts.useIdentityMetric = True
opts.beta = 1e-4
opts.gradTol = 1e-4
opts.niter = 1000



# Fix the boundary positions
import boundaries
bdryVars = boundaries.getBoundaryVars(isheet)
fixedVars = bdryVars


viewer.update(scalarField=utils.getStrains(isheet)[:, 0])    



isheet.setUseTensionFieldEnergy(True)

isheet.setUseHessianProjectedEnergy(False)

fixedVars, hessianShift = bdryVars, 1e-6

framerate = 20
def cb(it):
    if it % framerate == 0:
        viewer.update(scalarField=utils.getStrains(isheet)[:, 0])    


# ### First solve with low pressure to get out of indefinite state


isheet.pressure = 1e-5



opts.niter = 5

import time
cr = inflation.inflation_newton(isheet, fixedVars, opts, hessianShift = hessianShift, callback = cb)


# ### Then inflate


import time
benchmark.reset()
cr = inflation.inflation_newton(isheet, fixedVars, opts, hessianShift = hessianShift, callback = cb)
benchmark.report()

isheet.tensionStateHistogram()



# Plot maximum tensile strains in the sheet to verify the pressure is reasonable
from matplotlib import pyplot as plt
plt.hist(utils.getStrains(isheet)[:, 0], bins=1000);
plt.xlim(-0.04, 0.1);


# ## Shape Optimization


# Reset the inflation and set up target-attraction forces
isheet.setUninflatedDeformation(liftedSheetPositions.transpose(), prepareRigidMotionPinConstraints=False)
targetAttractedSheet = inflation.TargetAttractedInflation(isheet, target_surf)
targetAttractedSheet.energy(targetAttractedSheet.EnergyType.Fitting)



targetAttractedSheet.targetSurfaceFitter().holdClosestPointsFixed = True
targetAttractedSheet.fittingWeight = 1e-5



# Re-inflate, this time applying target-attraction forces.
import time
isheet.pressure = 0.025


benchmark.reset()
cr = inflation.inflation_newton(targetAttractedSheet, fixedVars, opts, hessianShift = hessianShift, callback = cb)
benchmark.report()



# Set up the sheet optimizer
import sheet_optimizer, opt_config
origDesignMesh = isheet.mesh().copy()

def config(so):
    fcs = so.rso.fusingCurveSmoothness()
    fcs.interiorWeight = 1/10


sheet_opt = sheet_optimizer.PySheetOptimizer(targetAttractedSheet, fixedVars, renderMode=sheet_optimizer.RenderMode.OFFSCREEN, screenshotPath=optimization_data_path + '{}.mp4'.format(name),
                                            detActivationThreshold=0.9, detActivationThresholdTubeTri=0.5,
                                            originalDesignMesh=origDesignMesh, fusingCurveSmoothnessConfig=opt_config.FusingCurveSmoothnessParams(0.0, 0.0, 2.0, 2.0), customConfigCallback=config)

# Configure some more weights
sheet_opt.rso.compressionPenaltyWeight = 1e-6
fcs = sheet_opt.rso.fusingCurveSmoothness()
fcs.interiorWeight = 0.05


sheet_opt.flat_viewer.showWireframe()

# cashew
sheet_opt.deploy_viewer.setCameraParams(((1.6079927490763757, -1.3197729912161424, 4.546708542514496),
 (-0.2864801503602479, 0.8882263079150637, 0.3591422411484639),
 (0.0, 0.0, 0.0)))

# # igloo
# sheet_opt.deploy_viewer.setCameraParams(((1.8419607977129837, -3.919084387618971, 2.4995915631178147),
#  (-0.2492987224778718, 0.434753560302595, 0.8653551229264879),
#  (0.0, 0.0, 0.0)))

# sheet_opt.flat_viewer.setCameraParams(((-0.005724654779466134, -0.15742800643532207, 3.663779058858404),
#     (0.0, 1.0, 0.0),
#     (-0.005724654779466134, -0.15742800643532207, 0.0)))
sheet_opt.flat_viewer.showWireframe()
sheet_opt.deploy_viewer.scalarFieldGetter = visualization.ISheetScalarField.TGT_DIST(targetAttractedSheet.sheet(), target_surf)

sheet_opt.setSolver(sheet_optimizer.Solver.SCIPY, 500)
sheet_opt.rso.getEquilibriumSolver().options.niter = 20


sheet_opt.optimize()


sheet_opt.save('{}/optimized_sheet_opt.pkl.gz'.format(optimization_data_path))



scaleFactor = 1 # Factor for fine-tuning size to fit the machine's build area
channelMargin = 0
tabMargin = 0



isheet = sheet_opt.rso.sheet()
optMesh = sheet_opt.rso.mesh().copy()
origMesh = sheet_opt.rso.originalMesh().copy()
import inflation
tas = sheet_opt.rso.targetAttractedInflation()
tsf = tas.targetSurfaceFitter()
targetSurf = mesh.Mesh(tsf.targetSurfaceV, tsf.targetSurfaceF)
iwv = [isheet.isWallVtx(i) for i in range(isheet.mesh().numVertices())]


import fabrication
importlib.reload(fabrication)
fabrication.writeFabricationData('{}/fabrication_data/{}/parallel_tube'.format(optimization_data_path, name), origMesh, optMesh, iwv, targetSurf, uv,
                                 scale=scaleFactor, numTabs=20, inletOffset=0, tabOffset=0.6 / 20,
                                 channelMargin=channelMargin, tabMargin=tabMargin, tabWidth=5, tabHeight=8, fuseSeamWidth=0.01, inletScale=0,
                                 overlap=0.0, smartOuterChannel=True)









