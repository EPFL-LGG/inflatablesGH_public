import sys; sys.path.append('..')
import os
import MeshFEM, mesh, sparse_matrices, benchmark, field_sampler, mesh_utilities, parallelism
import inflatables_parametrization as parametrization, numpy as np, importlib, pickle, wall_generation, visualization
import utils
import py_newton_optimizer
from py_newton_optimizer import NewtonOptimizerOptions
from io_redirection import suppress_stdout
from py_newton_optimizer import NewtonOptimizerOptions
import wall_width_formulas as wwf
import numpy as np
from matplotlib import pyplot as plt
from tri_mesh_viewer import TriMeshViewer

parallelism.set_max_num_tbb_threads(6)

################################################################################

threshold = float(sys.argv[1])
out_dir = f'cap_{threshold}'
os.makedirs(out_dir, exist_ok=True)

################################################################################

sphere = mesh.Mesh('../../examples/full_sphere.msh')
sphere = mesh_utilities.subdivide_loop(sphere, 1)
################################################################################
bb = utils.bbox(sphere.vertices())
zThreshold = bb[1][2] * threshold + (1 - threshold) * bb[0][2]

# Delete the vertices below the Z threshold and then construct the triangles
# induced by this vertex subset.
keepVtx = sphere.vertices()[:, 2] >= zThreshold
V = sphere.vertices()[keepVtx, :]

vtxRenumber = -1 * np.ones(sphere.numVertices())
vtxRenumber[np.arange(sphere.numVertices())[keepVtx]] = np.arange(V.shape[0])

F = vtxRenumber[sphere.triangles()]
F = F[np.min(F, axis=1) >= 0]
target_surf = mesh.Mesh(V, F)
target_surf.setVertices(utils.prototypeScaleNormalization(target_surf.vertices(), placeAtopFloor=True))
################################################################################

benchmark.reset()

# Choose reasonable stretching bounds
alphaMin = wwf.stretchFactorForCanonicalWallWidth(wwf.canonicalWallWidthForGeometry(2, 10))
alphaMax = wwf.stretchFactorForCanonicalWallWidth(wwf.canonicalWallWidthForGeometry(1, 10))
print(alphaMin, alphaMax)

lg = parametrization.LocalGlobalParametrizer(target_surf, parametrization.lscm(target_surf))

for i in range(1000): lg.runIteration()
print(lg.energy())
lg.alphaMin = alphaMin
lg.alphaMax = alphaMax

print(lg.energy())
for i in range(1000): lg.runIteration()
print(lg.energy())


print(lg.energy())
for i in range(8000): lg.runIteration()
print(lg.energy())

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

for _ in range(3):
    with suppress_stdout(): optimize_rparam(rparam, 1e2, 1e0, 50.0)
    with suppress_stdout(): optimize_rparam(rparam, 1e1, 1e-1, 50.0)

    with suppress_stdout(): optimize_rparam(rparam, 1e0, 1e-2, 250.0)
    with suppress_stdout(): optimize_rparam(rparam, 1e-1, 5e-3, 100.0)
    with suppress_stdout(): optimize_rparam(rparam, 1e-1, 5e-3, 25.0)
    with suppress_stdout(): optimize_rparam(rparam, 1e-2, 2.5e-3, 12.5)

benchmark.report()

utils.save(rparam.uv(), f'{out_dir}/uv.pkl.gz')

visualization.singularValueHistogram(rparam)
plt.savefig(f'{out_dir}/sv_hist.pdf')
plt.close()

visualization.visualize(rparam)
plt.savefig(f'{out_dir}/rparam.png')
plt.close()

# ## Upsampling and channel generation

benchmark.reset()
nsubdiv=2
upsampledMesh, upsampledAngles, upsampledStretches = rparam.upsampledVertexLeftStretchAnglesAndMagnitudes(nsubdiv)
upsampledStretches = np.clip(upsampledStretches, alphaMin, alphaMax)
(sdfVertices, sdfTris, sdf) = wall_generation.evaluate_stripe_field(upsampledMesh.vertices(), upsampledMesh.triangles(), upsampledAngles,
                                                                    wwf.canonicalWallWidthForStretchFactor(upsampledStretches), frequency=0.2, nsubdiv=3)

visualization.scalarFieldPlotFast(sdfVertices, sdfTris, sdf, height=12)

pts, edges = wall_generation.extract_contours(sdfVertices, sdfTris, sdf,
                                              targetEdgeSpacing=4.0,
                                              minContourLen=10)

visualization.plot_line_segments(pts, edges, width=20, height=16)
plt.savefig(f'{out_dir}/segments.png')
plt.close()

triArea = 12.0
import sheet_meshing
m, iwv, iwbv = sheet_meshing.newMeshingAlgorithm(sdfVertices, sdfTris, sdf, pts, edges, triArea=triArea)
visualization.plot_2d_mesh(m, pointList=iwv, width=20, height=18)
plt.savefig(f'{out_dir}/mesh.png')
plt.close()

benchmark.report()

## Meshing and inflation simulation

import sparse_matrices, mesh, numpy as np, importlib, pickle, wall_generation
import inflation
import vis, matplotlib
from py_newton_optimizer import NewtonOptimizerOptions
from numpy.linalg import norm
from io_redirection import suppress_stdout
import visualization, wall_width_formulas as wwf

benchmark.reset()

isheet = inflation.InflatableSheet(m, iwv)
uv = rparam.uv()

isheet.setIdentityDeformation()
paramSampler = field_sampler.FieldSampler(np.pad(uv, [(0, 0), (0, 1)], 'constant'), target_surf.triangles())
liftedSheetPositions = paramSampler.sample(m.vertices(), target_surf.vertices())

import py_newton_optimizer
niter = 2000
iterations_per_output = 10
opts = py_newton_optimizer.NewtonOptimizerOptions()
opts.useIdentityMetric = True
opts.beta = 1e-4
opts.gradTol = 1e-7
opts.niter = iterations_per_output

isheet.setUseTensionFieldEnergy(True)


utils.allEnergies(isheet)


import opt_config
fixedVars = opt_config.FixedVarsBoundary.get(isheet)
freeVars = ~utils.maskForIndexList(fixedVars, isheet.numVars())

isheet.setRelaxedStiffnessEpsilon(1e-6)


isheet.setUninflatedDeformation(liftedSheetPositions.transpose(), prepareRigidMotionPinConstraints=False)
targetAttractedSheet = inflation.TargetAttractedInflation(isheet, target_surf)
targetAttractedSheet.energy(targetAttractedSheet.EnergyType.Fitting)

targetAttractedSheet.targetSurfaceFitter().holdClosestPointsFixed = True
targetAttractedSheet.fittingWeight = 1e-5
targetAttractedSheet.fittingWeight = 0.0

import time
isheet.pressure = 0.025
for step in range(int(niter / iterations_per_output)):
    cr = inflation.inflation_newton(targetAttractedSheet, fixedVars, opts)
    if cr.numIters() < iterations_per_output: break
    # viewer.update(False, isheet.visualizationMesh())
    time.sleep(0.05) # Allow some mesh synchronization time for pythreejs
benchmark.report()

# # Plot maximum tensile strains in the sheet
# plt.hist(utils.getStrains(isheet)[:, 0], bins=1000);
# plt.xlim(-0.04, 0.06);
# plt.savefig(f'{out_dir}/max_strains_init.pdf')
# plt.close()
# 
# 
# # Plot maximum tensile strains in the sheet
# plt.hist(utils.getStrains(isheet)[:, 1], bins=1000);
# plt.xlim(-0.5, 0.06);
# plt.savefig(f'{out_dir}/min_strains_init.pdf')
# plt.close()

utils.allEnergies(targetAttractedSheet)
utils.allGradientNorms(targetAttractedSheet, freeVariables=freeVars)


import sheet_optimizer, opt_config
origDesignMesh = isheet.mesh().copy()

benchmark.reset()
sheet_opt = sheet_optimizer.PySheetOptimizer(targetAttractedSheet, fixedVars, renderMode=sheet_optimizer.RenderMode.NONE,
                                             detActivationThreshold=0.9, detActivationThresholdTubeTri=0.5,
                                             originalDesignMesh=origDesignMesh, fusingCurveSmoothnessConfig=opt_config.FusingCurveSmoothnessParams(0.0, 0.0, 1.0, 1.0))

sheet_opt.rso.fusingCurveSmoothness().interiorWeight = 1 / 10

# sheet_opt.deploy_viewer.setCameraParams(((-2.8127683958467866, -3.179300313721665, 1.274256905127606),
#                                          (0.24888447315371076, 0.1624925314484115, 0.9548050566720383),
#                                          (0.0, 0.0, 0.0)))

sheet_opt.optimize()

benchmark.report()
sheet_opt.save(f'{out_dir}/opt.pkl.gz')
