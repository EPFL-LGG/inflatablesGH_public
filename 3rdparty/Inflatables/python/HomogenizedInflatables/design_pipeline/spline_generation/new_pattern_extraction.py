#!/usr/bin/env python
# coding: utf-8


import numpy as np
import sys
sys.path.append("../../")
sys.path.append("../../Visualization/")



experiment_file = '../../experiments/parallelized_experiments/output/cosine_curve_amplitude_full_period/2023_12_17_15_52//experiment_result.json'

stiffness_path = '../../experiments/parallelized_experiments/output/cosine_curve_amplitude_full_period/2023_12_17_15_52/'


# ### Overview


import pandas as pd
import json
import matplotlib.pyplot as plt
import numpy as np



with open(experiment_file, 'r') as fp:
    data = json.load(fp)



df = pd.DataFrame(data['data'])
fig, axes = plt.subplots(nrows = 1, ncols = 3, figsize = (20, 5))
a = (df.hist('Ipu simulation succeed', ax = axes[0]), df.hist('Planar equilibrium', ax = axes[1]), df.hist('Simulation Kappa value', ax = axes[2]))



import visualize_stiffness
import importlib
importlib.reload(visualize_stiffness)



valid_tags = np.array(df['name'][df['Planar equilibrium'] == 1])



kappa_path = None



name = 'cosine_curve_amplitude_full_period'



bending_stiffness_data, stretching_stiffness_data, scale_factor_data, used_tags = visualize_stiffness.plot_all_data(kappa_path, stiffness_path, name, valid_tags, plot_data = False)



parameters = (np.array(data['pattern_parameters'][0]['values']))



max_bending_stiffness = np.max(bending_stiffness_data, axis = 1)
min_bending_stiffness = np.min(bending_stiffness_data, axis = 1)
max_stretching_stiffness = np.max(stretching_stiffness_data, axis = 1)
min_stretching_stiffness = np.min(stretching_stiffness_data, axis = 1)



parameters[np.argmax(min_bending_stiffness)]



x_scale_factors, y_scale_factors = visualize_stiffness.get_axis_scale_factors(stiffness_path, name, valid_tags)



min_scale_factors = np.min(np.concatenate((x_scale_factors.reshape(-1, 1), y_scale_factors.reshape(-1, 1)), axis = 1), axis = 1)
max_scale_factors = np.max(np.concatenate((x_scale_factors.reshape(-1, 1), y_scale_factors.reshape(-1, 1)), axis = 1), axis = 1)



angle_offsets = visualize_stiffness.get_max_flattening_factor_offset(stiffness_path, name, valid_tags)


# ### Get scale function convex hull


import matplotlib.cm as cm
import matplotlib as mpl



from scipy.spatial import ConvexHull, convex_hull_plot_2d
import numpy as np
rng = np.random.default_rng()
points = rng.random((30, 2))   # 30 random points in 2-D
# points = np.concatenate((min_scale_factor.reshape((-1, 1)), max_scale_factor.reshape((-1, 1))), axis = 1)

points = np.concatenate((max_scale_factors.reshape((-1, 1)), min_scale_factors.reshape((-1, 1))), axis = 1)
hull = ConvexHull(points)



hull


# ### Validate the max and min scale factors are aligned with the x and y axis


import visualize_stiffness
importlib.reload(visualize_stiffness)



eqns = hull.equations



hull.max_bound, hull.min_bound



import parametrization_helper, importlib
importlib.reload(parametrization_helper)



parametrization_helper.visualize_scale_factors(eqns, max_scale_factors, min_scale_factors)


# ### Generate data without augmenting


import visualize_stiffness
importlib.reload(visualize_stiffness)



stiffness_coefficients = np.array(visualize_stiffness.get_stiffness_coefficients(stiffness_path, name, (used_tags)))
# For patches with reflection symmetry:
stiffness_coefficients[:, 1] *= 0
stiffness_coefficients[:, 2] *= 0



# for i in range(5):
#     for j in range(30):
#         stiffness_coefficients[:, i] = parametrization_helper.savitzky_golay(stiffness_coefficients[:, i], 11, 3) # window size 51, polynomial order 3



np.set_printoptions(suppress=True, precision=4)



np.argmax(stiffness_coefficients[:, 1]), np.argmax(stiffness_coefficients[:, 2])



def get_stiffness_polynomial(s, theta):
    return s[0] * np.cos(theta)**2 * np.sin(theta)**2 + s[1] * np.cos(theta)**3 * np.sin(theta) + s[2] * np.cos(theta) * np.sin(theta)**3 + s[3] * np.cos(theta)**4 + s[4] * np.sin(theta)**4



grid_data = np.zeros((9, len(parameters)))



for i in range(len(parameters)):
    grid_data[0][i] = max_scale_factors[i]
    grid_data[1][i] = min_scale_factors[i]
    grid_data[2][i] = x_scale_factors[i]
    grid_data[3][i] = y_scale_factors[i]
    for s in range(5):
        grid_data[4 + s][i] = stiffness_coefficients[i][s]



# np.save("grid_pattern_1.npy", grid_pattern_1)
# np.save("grid_pattern_2.npy", grid_pattern_2)
# np.save("grid_data.npy", grid_data)



importlib.reload(parametrization_helper)



splines = parametrization_helper.ndsplines_get_mat_params_over_pattern_params_grid_interpolation(grid_data, (parameters))



grid_data.shape



scale_factors_grid_data = np.zeros((2, len(parameters)))
for i in range(len(parameters)):
    scale_factors_grid_data[0][i] = x_scale_factors[i]
    scale_factors_grid_data[1][i] = y_scale_factors[i]
scale_factors_splines = parametrization_helper.ndsplines_get_mat_params_over_pattern_params_grid_interpolation(scale_factors_grid_data, (parameters))



test_parameters = np.linspace(0, 0.9, 100)



fig, axes = plt.subplots(1, 9, figsize=(45, 8))
titles = ['max scale factors', 'min scale factors', 'x scale factors', 'y scale factors', 's1', 's2', 's3', 's4', 's5']
# titles = ['max scale factors', 'min scale factors', 's1', 's2', 's3', 's4', 's5']

for i in range(9):
    axes[i].plot(test_parameters, splines[i * 3 + 1](test_parameters))
    axes[i].set_title(titles[i], fontsize=21)



stiffness_coefficients = np.array(stiffness_coefficients)



stiffness_coefficients.shape



fig, axes = plt.subplots(1, 9, figsize=(45, 8))
titles = ['max scale factors', 'min scale factors', 'x scale factors', 'y scale factors', 's1', 's2', 's3', 's4', 's5']
data = [max_scale_factors, min_scale_factors, x_scale_factors, y_scale_factors, stiffness_coefficients[:, 0], stiffness_coefficients[:, 1], stiffness_coefficients[:, 2], stiffness_coefficients[:, 3], stiffness_coefficients[:, 4]]

for i in range(9):
    axes[i].plot(parameters, data[i])
    axes[i].set_title(titles[i], fontsize=21)


# ### End data generating

# ### Parametrization


import sys; sys.path.append('../../../'); sys.path.append('../../../periodic_patches/'); sys.path.append('../../experiments/'); sys.path.append('../../../gmsh')
import inflation, sparse_matrices, mesh, numpy as np, importlib, pickle
import inflatables_parametrization as parametrization
from numpy.linalg import norm
from io_redirection import suppress_stdout
import visualization



sys.path.append('periodic_patches/')
sys.path.append('gmsh')



from periodic_simulation_setup import *



from py_newton_optimizer import NewtonOptimizerOptions



import MeshFEM, parallelism, benchmark, utils
parallelism.set_max_num_tbb_threads(32)
parallelism.set_gradient_assembly_num_threads(32)
parallelism.set_hessian_assembly_num_threads(32)



import utils, mesh_utilities
importlib.reload(utils)



target_surf = mesh.Mesh("../../../../examples/igloo.obj")
target_surf.setVertices(utils.prototypeScaleNormalization(target_surf.vertices(), placeAtopFloor=False))
# target_surf = mesh_utilities.subdivide_loop(target_surf, 1)



lines = np.array(eqns)


# ### New local global with convex hull


lg = parametrization.LocalGlobalGenericParametrizer(target_surf, parametrization.lscm(target_surf))

lg.setLines(eqns)

lg.alphaMin = hull.min_bound[0]
lg.alphaMax = hull.max_bound[0]

lg.betaMin = hull.min_bound[1]
lg.betaMax = hull.max_bound[1]

print(lg.energy())
for i in range(1): lg.runIteration()

print(lg.energy())
lg.runIteration()
print(lg.energy())


default_pattern_params = [0.2]  * len(lg.getAlphas())



mat_info = np.array(default_pattern_params).reshape((1, len(lg.getAlphas())))



rparam = parametrization.RegularizedPatternParametrizer(lg, splines, default_pattern_params, len(grid_data.shape) - 1)
rparam.patternParamBounds = np.array([[0.03, 0.4]])
rparam.diffRegW = 0.0



PET = parametrization.RegularizedPatternParametrizer.PatternEnergyType



rparam.bendRegW = 1




PET = parametrization.RegularizedPatternParametrizer.PatternEnergyType
list(map(rparam.energy, [PET.Full, PET.RGP, PET.Bending]))



def optimize_rparam(param, patternRegW, phiRegW, bendRegW = 0.0, update_uv = True, niter = 100):
    param.patternRegW = patternRegW
    param.phiRegW = phiRegW
    param.bendRegW = bendRegW
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



PET = parametrization.RegularizedPatternParametrizer.PatternEnergyType
list(map(rparam.energy, [PET.Full, PET.PatternRegularization, PET.Bending, PET.DEBUG_Fitting, PET.DEBUG_StretchRegularization, PET.DEBUG_PhiRegularization]))



benchmark.reset()
with suppress_stdout(): report = optimize_rparam(rparam, 0, 0, bendRegW = 0, update_uv = False, niter = 2)
benchmark.report()



PET = parametrization.RegularizedPatternParametrizer.PatternEnergyType
list(map(rparam.energy, [PET.Full, PET.PatternRegularization, PET.Bending, PET.DEBUG_Fitting, PET.DEBUG_StretchRegularization, PET.DEBUG_PhiRegularization]))

importlib.reload(visualization)
importlib.reload(parametrization_helper)



rparam.get_stretch_angle_offset_from_pattern_params(rparam.getMatInfoArgs())



max(parametrization_helper.get_stretch_angle_offset_from_pattern_params(scale_factors_splines, rparam.getPatternParams().reshape(1, -1)))



rparam.get_stretch_angle_offset_from_pattern_params(rparam.getPatternParams().reshape(1, -1))


# ## Upsampling and channel generation


nsubdiv=4
upsampledMesh, upsampledAngles, upsampledPatternParams = rparam.upsampledVertexLeftStretchAnglesAndPatternParameters(nsubdiv)
upsampleMesh_vertices = upsampledMesh.vertices()
upsampleMesh_triangles = upsampledMesh.triangles()
angle_offset = parametrization_helper.get_stretch_angle_offset_from_pattern_params(scale_factors_splines, np.array(upsampledPatternParams).reshape((1, -1)))
upsampledAngles += angle_offset



# min(angle_offset), max(angle_offset)



np.save("upsampleMesh_vertices.npy", upsampledMesh.vertices())
np.save("upsampleMesh_triangles.npy", upsampledMesh.triangles())
np.save("upsampleAngles.npy", upsampledAngles)
np.save("upsampledPatternParams.npy", upsampledPatternParams)



upsampleMesh_vertices = np.load("upsampleMesh_vertices.npy")
upsampleMesh_triangles = np.load("upsampleMesh_triangles.npy")
upsampleAngles = np.load("upsampleAngles.npy")
upsampledPatternParams = np.load("upsampledPatternParams.npy")



np.set_printoptions(suppress=True)



max(upsampledPatternParams[0])



amp_data =  upsampledPatternParams[0]



import igl



from parametrization_helper import get_distance_to_line_segments



def fusing_curve_polyline(patternParams):
#     Draw cosine curves.
    amp = patternParams[0]
    def get_y_from_x(x):
        return amp * np.cos(x) * 0.5 * np.pi + np.pi / 2
    
    x_coords = np.linspace(-np.pi, np.pi, 10)
    y_coords = get_y_from_x(x_coords)
    x_coords += np.pi
    x_coords /= 2
    polyline = np.concatenate(((y_coords).reshape(-1, 1), (x_coords).reshape(-1, 1)), axis = 1)
    return polyline



def pattern_function(theta, gamma, patternParams, margin, draw_boundary = False):
    if draw_boundary:
#         This is for debugging only and shouldn't be used for generating the inflatable mesh.
        if (theta < 0.1):
            return - margin
        if (gamma < 0.1):
            return - margin
    # Gamma is y, theta is x
    # Gamma theta are between 0 and pi
    polyline = fusing_curve_polyline(patternParams)    
    polyline_dist = get_distance_to_line_segments(np.array([theta, gamma]), polyline)
    return polyline_dist - margin



def pattern_polyline_function(theta_gamma, patternParams):
    return [([0, 0, 1], [0, 1, 0])]



pattern_function(np.pi / 4, np.pi / 10, [0.], 0.)



# upsampledPatternParams = np.ones_like(upsampledPatternParams) * 0.4



import time
start_time = time.time()
(sdfVertices, sdfTris, sdf, edge_soup) = wall_generation.evaluate_cross_field_custom_pattern(upsampleMesh_vertices, upsampleMesh_triangles, upsampleAngles, upsampledPatternParams, pattern_function, pattern_polyline_function, frequency=0.1, margin = 0.07)
print(time.time() - start_time)

# pickle.dump((sdfVertices, sdfTris, sdf), open('stripe_sdf_ns4_f100.pkl', 'wb'))

# import pickle, mesh, wall_generation, visualization, numpy as np
# (sdfVertices, sdfTris, sdf) = pickle.load(open('stripe_sdf_ns4_f100.pkl', 'rb'))



importlib.reload(visualization)

import matplotlib as mpl



visualization.scalarFieldPlotFast(sdfVertices, sdfTris, sdf, width = 10, height=10)



importlib.reload(visualization)
visualization.scalarFieldPlotZeroContourFast(sdfVertices, sdfTris, sdf, width = 15, height=10, cmap = mpl.colormaps["PiYG"])



pts, edges = wall_generation.extract_contours(sdfVertices, sdfTris, sdf,
                                              targetEdgeSpacing=0.4,
                                              minContourLen=1)

visualization.plot_line_segments(pts, edges, width=15, height=15)


# ## Meshing and inflation simulation


import sheet_meshing, inflation



m, iwv, iwbv = sheet_meshing.newMeshingAlgorithm(sdfVertices, sdfTris, sdf, pts, edges, triArea=1e0)



visualization.plot_2d_mesh(m, pointList=np.where(np.array(iwv) == 1)[0], width=10, height=10)



import inflation
isheet = inflation.InflatableSheet(m, np.array(iwv) != 0)



from mesh_utilities import SurfaceSampler, tubeRemesh


paramSampler = SurfaceSampler(np.pad(rparam.uv(), [(0, 0), (0, 1)], 'constant'), target_surf.triangles())
liftedSheetPositions = paramSampler.sample(m.vertices(), target_surf.vertices())

isheet.setUninflatedDeformation(liftedSheetPositions.transpose())

isheet.getVars()



import py_newton_optimizer
niter = 2000
iterations_per_output = 10
opts = py_newton_optimizer.NewtonOptimizerOptions()
opts.useIdentityMetric = True
opts.beta = 1e-4
opts.gradTol = 1e-10
opts.niter = iterations_per_output



# Are the flat region causing a problem? They might not actually control the metric...
# Try replacing them with single wall...
# Analyze the actual stretching factor (much easier to do with skeleton walls)

from tri_mesh_viewer import TriMeshViewer
viewer = TriMeshViewer(isheet, width=768, height=640)
viewer.showWireframe(True)

viewer.show()



viewer.update(scalarField=utils.getStrains(isheet)[:, 0])    



isheet.setUseTensionFieldEnergy(True)

isheet.setUseHessianProjectedEnergy(False)

fixedVars, hessianShift = [], 1e-6

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


isheet.pressure = 1e-2



opts.niter = 2000
opts.gradTol = 1e-7

import time
benchmark.reset()
cr = inflation.inflation_newton(isheet, fixedVars, opts, hessianShift = hessianShift, callback = cb)
benchmark.report()

isheet.tensionStateHistogram()



# Plot maximum tensile strains in the sheet to verify the pressure is reasonable
from matplotlib import pyplot as plt
plt.hist(utils.getStrains(isheet)[:, 0], bins=1000);
plt.xlim(-0.04, 0.1);



import gzip



pickle.dump(isheet,  gzip.open("igloo_pattern_optimized_2023_12_17_low_frequency_with_bending_high_resolution.pkl.gz", 'wb'))


# ### Generate Fabrication Files


scaleFactor = 1 # Factor for fine-tuning size to fit the machine's build area
channelMargin = 0 / scaleFactor # 8mm channel margin
tabMargin = 2 / scaleFactor # 2mm tab margin



import inflation

targetSurf = target_surf
iwv = [isheet.isWallVtx(i) for i in range(isheet.mesh().numVertices())]



isheet.mesh()



mesh = isheet.mesh()



mesh.save("igloo_2D.obj")



np.save("igloo_is_wall.npy", iwv)







uv = rparam.uv()



# !pip install shapely==1.7.0



import shapely



import fabrication



importlib.reload(utils)
importlib.reload(fabrication)
importlib.reload(shapely)
import shapely.geometry as shp



shapely.__version__



# new_fabrication.writeFabricationData('fabrication_data/igloo/free_bdry', isheet.mesh(), iwv, targetSurf, uv,
#                                  scale=scaleFactor,
#                                  channelMargin=channelMargin, fuseSeamWidth=None,
#                                  overlap=0.0, smartOuterChannel=True)



import fabrication
fabrication.writeFabricationData('fabrication_data/igloo_2023_12_17/free_bdry_low_frequency', isheet.mesh(), isheet.mesh(), iwv, targetSurf, uv,
                                 scale=scaleFactor, numTabs=0, inletOffset=0, tabOffset=0.60 / 80,
                                 channelMargin=channelMargin, tabMargin=tabMargin, tabWidth=5, tabHeight=8, fuseSeamWidth=None, inletScale=None,
                                 overlap=0.0, smartOuterChannel=False)



# import fabrication
# fabrication.writeFabricationData('fabrication_data/igloo/free_bdry', isheet.mesh(), isheet.mesh(), iwv, targetSurf, uv,
#                                  scale=scaleFactor, numTabs=0, inletOffset=0.742, tabOffset=0.60 / 80,
#                                  channelMargin=channelMargin, tabMargin=tabMargin, tabWidth=5, tabHeight=8, fuseSeamWidth=1.0, inletScale=12 / channelMargin / scaleFactor,
#                                  overlap=0.0, smartOuterChannel=True)



from IPython.display import SVG, display



display(SVG(filename = "fabrication_data/igloo/free_bdry_low_frequency/orig.wall_boundaries.svg"))









