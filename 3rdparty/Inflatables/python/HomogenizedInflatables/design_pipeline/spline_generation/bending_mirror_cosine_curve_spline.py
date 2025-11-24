#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import sys
sys.path.append("../../")
sys.path.append("../../Visualization/")


# In[2]:


experiment_file = '../../experiments/parallelized_experiments/output/cosine_curve_amplitude_full_period/2023_12_08_15_48/experiment_result.json'

stiffness_path = '../../experiments/parallelized_experiments/output/cosine_curve_amplitude_full_period/2023_12_08_15_48/'


# ### Overview

# In[3]:


import pandas as pd
import json
import matplotlib.pyplot as plt
import numpy as np


# In[4]:


with open(experiment_file, 'r') as fp:
    data = json.load(fp)


# In[5]:


df = pd.DataFrame(data['data'])
fig, axes = plt.subplots(nrows = 1, ncols = 3, figsize = (20, 5))
a = (df.hist('Ipu simulation succeed', ax = axes[0]), df.hist('Planar equilibrium', ax = axes[1]), df.hist('Simulation Kappa value', ax = axes[2]))


# In[6]:


import visualize_stiffness
import importlib
importlib.reload(visualize_stiffness)


# In[7]:


valid_tags = np.array(df['name'][df['Planar equilibrium'] == 1])


# In[8]:


kappa_path = None


# In[9]:


name = 'cosine_curve_amplitude_full_period'


# In[10]:


bending_stiffness_data, stretching_stiffness_data, scale_factor_data, used_tags = visualize_stiffness.plot_all_data(kappa_path, stiffness_path, name, valid_tags, plot_data = False)


# In[11]:


parameters = (np.array(data['pattern_parameters'][0]['values']))


# In[12]:


max_bending_stiffness = np.max(bending_stiffness_data, axis = 1)
min_bending_stiffness = np.min(bending_stiffness_data, axis = 1)
max_stretching_stiffness = np.max(stretching_stiffness_data, axis = 1)
min_stretching_stiffness = np.min(stretching_stiffness_data, axis = 1)


# In[13]:


x_scale_factors, y_scale_factors = visualize_stiffness.get_axis_scale_factors(stiffness_path, name, valid_tags)


# In[14]:


min_scale_factors = np.min(np.concatenate((x_scale_factors.reshape(-1, 1), y_scale_factors.reshape(-1, 1)), axis = 1), axis = 1)
max_scale_factors = np.max(np.concatenate((x_scale_factors.reshape(-1, 1), y_scale_factors.reshape(-1, 1)), axis = 1), axis = 1)


# In[15]:


angle_offsets = visualize_stiffness.get_max_flattening_factor_offset(stiffness_path, name, valid_tags)


# In[16]:


parameters


# In[17]:


parameters[(np.where(min_bending_stiffness > 0.1))]


# In[18]:


parameters[(np.where(min_stretching_stiffness > 4))]


# ### Get scale function convex hull

# In[19]:


import matplotlib.cm as cm
import matplotlib as mpl


# In[20]:


from scipy.spatial import ConvexHull, convex_hull_plot_2d
import numpy as np
rng = np.random.default_rng()
points = rng.random((30, 2))   # 30 random points in 2-D
# points = np.concatenate((min_scale_factor.reshape((-1, 1)), max_scale_factor.reshape((-1, 1))), axis = 1)

points = np.concatenate((max_scale_factors.reshape((-1, 1)), min_scale_factors.reshape((-1, 1))), axis = 1)
hull = ConvexHull(points)


# In[21]:


hull


# In[22]:



# In[23]:


import numpy.linalg as la

# Need to plot the patches over the min and max scale factors, so we can get the polygon that constrain the singular values
# The scale factors we are considering during the parametrization are from the flattening, so it's the change from the inflated state to the fabricated state, hence we need to take one over the factors we have from the average deformation gradient from homogenization.
# max_scale_factor = 1 / np.array(scale_factor_data)[:, 0]
# min_scale_factor = 1 / np.array(scale_factor_data)[:, 1]

fig, ax = plt.subplots(figsize = (10, 10))

# plt.scatter(x_scale_factor, y_scale_factor, label = data_info[i][1], s = 50, alpha = 0.3)
# plt.scatter(y_scale_factor, x_scale_factor, label = data_info[i][1], s = 50, alpha = 0.3)

# plt.scatter(y_scale_factor, x_scale_factor, label = data_info[i][1], s = 50, alpha = 0.8, c = min_stiffness)


points = np.concatenate((max_scale_factors.reshape((-1, 1)), min_scale_factors.reshape((-1, 1))), axis = 1)
hull = ConvexHull(points)

for simplex in hull.simplices:
    plt.plot(points[simplex, 0], points[simplex, 1], 'k-')

ax.title.set_text("Scale factors")
plt.xlabel("x scale factors")
plt.ylabel("y scale factors")

plt.scatter(max_scale_factors, min_scale_factors, label = 'min_stiffness', s = 200, alpha = 1, c = min_bending_stiffness)
# plt.scatter(x_scale_factors, y_scale_factors, label = 'max_stiffness', s = 200, alpha = 1, c = max_bending_stiffness)

# Plot x = y line
lims = [
np.min([ax.get_xlim(), ax.get_ylim()]),  # min of both axes
np.max([ax.get_xlim(), ax.get_ylim()]),  # max of both axes
]

# now plot both limits against eachother
ax.plot(lims, lims, 'k-', alpha=0.75, zorder=0)
ax.set_aspect('equal')
ax.set_xlim(lims)
ax.set_ylim(lims)
ax.legend()
fig.tight_layout()
plt.savefig('scale_factor_values_{}.png'.format(name), dpi = 300)


# ### Validate the max and min scale factors are aligned with the x and y axis

# In[24]:


import visualize_stiffness
importlib.reload(visualize_stiffness)


# In[25]:


eqns = hull.equations


# In[26]:


hull.max_bound, hull.min_bound


# In[27]:


import parametrization_helper, importlib
importlib.reload(parametrization_helper)


# In[28]:


parametrization_helper.visualize_scale_factors(eqns, max_scale_factors, min_scale_factors)


# ### Generate data without augmenting

# In[29]:


import visualize_stiffness
importlib.reload(visualize_stiffness)


# In[30]:


stiffness_coefficients = visualize_stiffness.get_stiffness_coefficients(stiffness_path, name, (used_tags))


# In[31]:


grid_data = np.zeros((9, len(parameters)))


# In[32]:


for i in range(len(parameters)):
    grid_data[0][i] = max_scale_factors[i]
    grid_data[1][i] = min_scale_factors[i]
    grid_data[2][i] = x_scale_factors[i]
    grid_data[3][i] = y_scale_factors[i]
    for s in range(5):
        grid_data[4 + s][i] = stiffness_coefficients[i][s]


# In[33]:


# np.save("grid_pattern_1.npy", grid_pattern_1)
# np.save("grid_pattern_2.npy", grid_pattern_2)
# np.save("grid_data.npy", grid_data)


# In[34]:


importlib.reload(parametrization_helper)


# In[35]:


splines = parametrization_helper.ndsplines_get_mat_params_over_pattern_params_grid_interpolation(grid_data, (parameters))


# In[36]:


grid_data.shape


# In[37]:


scale_factors_grid_data = np.zeros((2, len(parameters)))
for i in range(len(parameters)):
    scale_factors_grid_data[0][i] = x_scale_factors[i]
    scale_factors_grid_data[1][i] = y_scale_factors[i]
scale_factors_splines = parametrization_helper.ndsplines_get_mat_params_over_pattern_params_grid_interpolation(scale_factors_grid_data, (parameters))


# In[38]:


test_parameters = np.linspace(0, 0.9, 100)


# In[39]:


fig, axes = plt.subplots(1, 9, figsize=(45, 8))
titles = ['max scale factors', 'min scale factors', 'x scale factors', 'y scale factors', 's1', 's2', 's3', 's4', 's5']
# titles = ['max scale factors', 'min scale factors', 's1', 's2', 's3', 's4', 's5']

for i in range(9):
    axes[i].plot(test_parameters, splines[i * 3](test_parameters))
    axes[i].set_title(titles[i], fontsize=21)


# ### End data generating

# ### Parametrization

# In[40]:


import sys; sys.path.append('../../../'); sys.path.append('../../../periodic_patches/'); sys.path.append('../../experiments/'); sys.path.append('../../../gmsh')
import inflation, sparse_matrices, mesh, numpy as np, importlib, pickle
import inflatables_parametrization as parametrization
from numpy.linalg import norm
from io_redirection import suppress_stdout
import visualization


# In[41]:


sys.path.append('periodic_patches/')
sys.path.append('gmsh')


# In[42]:


from periodic_simulation_setup import *


# In[43]:


from py_newton_optimizer import NewtonOptimizerOptions


# In[44]:


import MeshFEM, parallelism, benchmark, utils
parallelism.set_max_num_tbb_threads(32)
parallelism.set_gradient_assembly_num_threads(32)
parallelism.set_hessian_assembly_num_threads(32)


# In[45]:


import utils, mesh_utilities
importlib.reload(utils)


# In[46]:


target_surf = mesh.Mesh("../../../../examples/igloo.obj")
target_surf.setVertices(utils.prototypeScaleNormalization(target_surf.vertices(), placeAtopFloor=False))
# target_surf = mesh_utilities.subdivide_loop(target_surf, 1)


# In[47]:


lines = np.array(eqns)


# ### SIGGRAPH 21 Local global

# In[48]:




# ### New local global with convex hull

# In[54]:


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


# In[55]:


lg.alphaMin, lg.alphaMax, lg.betaMin, lg.betaMax


# In[56]:


visualization.visualize_both(lg)


# In[57]:


parametrization_helper.visualize_scale_factors(eqns, lg.getAlphas(), lg.getBetas())


# ### Pattern parameters optimization

# In[58]:


default_pattern_params = [0.2]  * len(lg.getAlphas())


# In[59]:


mat_info = np.array(default_pattern_params).reshape((1, len(lg.getAlphas())))


# In[60]:


rparam = parametrization.RegularizedPatternParametrizer(lg, splines, default_pattern_params, len(grid_data.shape) - 1)
rparam.patternParamBounds = np.array([[0.03, 0.4]])
rparam.diffRegW = 0.0


# In[61]:


visualization.visualize_both(rparam, height = 4, showBarriers=True)


# In[62]:


PET = parametrization.RegularizedPatternParametrizer.PatternEnergyType


# In[63]:


rparam.bendRegW = 1


# In[64]:


rparam.energy(PET.RGP)


# In[65]:


rparam.energy(PET.Bending)


# In[66]:


PET = parametrization.RegularizedPatternParametrizer.PatternEnergyType
list(map(rparam.energy, [PET.Full, PET.RGP, PET.Bending]))


# In[67]:


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


# In[68]:


PET = parametrization.RegularizedPatternParametrizer.PatternEnergyType
list(map(rparam.energy, [PET.Full, PET.PatternRegularization, PET.Bending, PET.DEBUG_Fitting, PET.DEBUG_StretchRegularization, PET.DEBUG_PhiRegularization]))


# In[69]:


benchmark.reset()
with suppress_stdout(): report = optimize_rparam(rparam, 0, 0, bendRegW = 0, update_uv = False, niter = 200)
benchmark.report()


# In[70]:


PET = parametrization.RegularizedPatternParametrizer.PatternEnergyType
list(map(rparam.energy, [PET.Full, PET.PatternRegularization, PET.Bending, PET.DEBUG_Fitting, PET.DEBUG_StretchRegularization, PET.DEBUG_PhiRegularization]))


# In[71]:


parametrization_helper.visualize_scale_factors(eqns, rparam.getAlphas(), rparam.getBetas())


# In[72]:


importlib.reload(visualization)
visualization.visualize_both(rparam, height = 4, showBarriers=True)
visualization.visualize_pattern(rparam, height = 4, showBarriers=False, num_pattern_vars=1)


# In[73]:


benchmark.reset()
with suppress_stdout(): report = optimize_rparam(rparam, patternRegW = 0, phiRegW = 0, bendRegW = 0, update_uv = True, niter = 200)
benchmark.report()


# In[74]:


PET = parametrization.RegularizedPatternParametrizer.PatternEnergyType
list(map(rparam.energy, [PET.Full, PET.PatternRegularization, PET.Bending, PET.DEBUG_Fitting, PET.DEBUG_StretchRegularization, PET.DEBUG_PhiRegularization]))


# In[75]:


importlib.reload(visualization)
visualization.visualize_both(rparam, height = 4, showBarriers=True)
visualization.visualize_pattern(rparam, height = 4, showBarriers=False, num_pattern_vars=1)


# In[76]:


importlib.reload(parametrization_helper)
parametrization_helper.visualize_scale_factors(lines, rparam.getAlphas(), rparam.getBetas())


# In[77]:


# no_regularization_var = rparam.getVars()


# In[78]:


# two_separate_optimization_vars = rparam.getVars()


# In[79]:


benchmark.reset()
with suppress_stdout(): report = optimize_rparam(rparam, patternRegW = 0, phiRegW = 1e-5, bendRegW = 0, update_uv = True, niter = 100)
benchmark.report()

benchmark.reset()
with suppress_stdout(): report = optimize_rparam(rparam, patternRegW = 0, phiRegW = 1e-5, bendRegW = 0, update_uv = True, niter = 100)
benchmark.report()

benchmark.reset()
with suppress_stdout(): report = optimize_rparam(rparam, patternRegW = 0, phiRegW = 1e-5, bendRegW = 0, update_uv = True, niter = 100)
benchmark.report()


# In[80]:


PET = parametrization.RegularizedPatternParametrizer.PatternEnergyType
list(map(rparam.energy, [PET.Full, PET.PatternRegularization, PET.Bending, PET.DEBUG_Fitting, PET.DEBUG_StretchRegularization, PET.DEBUG_PhiRegularization]))


# In[81]:


importlib.reload(visualization)
visualization.visualize_both(rparam, height = 4, showBarriers=True)
visualization.visualize_pattern(rparam, height = 4, showBarriers=False, num_pattern_vars=1)


# In[82]:


# rparam.bendRegW = 1e-1
# PET = parametrization.RegularizedPatternParametrizer.PatternEnergyType
# list(map(rparam.energy, [PET.Full, PET.PatternRegularization, PET.Bending, PET.DEBUG_Fitting, PET.DEBUG_StretchRegularization, PET.DEBUG_PhiRegularization]))


# In[83]:


# benchmark.reset()
# with suppress_stdout(): report = optimize_rparam(rparam, patternRegW = 0, phiRegW = 1e-5, bendRegW = -1e-1, update_uv = True, niter = 100)
# benchmark.report()


# In[84]:


# rparam.bendRegW = 1e-1
# PET = parametrization.RegularizedPatternParametrizer.PatternEnergyType
# list(map(rparam.energy, [PET.Full, PET.PatternRegularization, PET.Bending, PET.DEBUG_Fitting, PET.DEBUG_StretchRegularization, PET.DEBUG_PhiRegularization]))


# In[85]:


# importlib.reload(visualization)
# visualization.visualize_both(rparam, height = 4, showBarriers=True)
# visualization.visualize_pattern(rparam, height = 4, showBarriers=False, num_pattern_vars=1)


# In[86]:


importlib.reload(parametrization_helper)
parametrization_helper.visualize_scale_factors(lines, rparam.getAlphas(), rparam.getBetas())


# In[88]:


visualization.visualizeChannelOrientation(rparam, quiver=visualization.QuiverVisualization.PER_VTX, orientationHue=False, width = 10, height = 10)


# In[91]:


visualization.visualizeChannelOrientationWithIsotropicPoints(rparam, quiver=visualization.QuiverVisualization.PER_TRI, orientationHue=False, width = 10, height = 10, use_x_axis=True)


# In[104]:


importlib.reload(visualization)
importlib.reload(parametrization_helper)


# In[105]:


rparam.get_stretch_angle_offset_from_pattern_params(rparam.getMatInfoArgs())


# In[106]:


scale_factors_splines[0](rparam.getMatInfoArgs()), scale_factors_splines[1](rparam.getMatInfoArgs())


# In[107]:


splines[3 * 2](rparam.getMatInfoArgs()), splines[3 * 3](rparam.getMatInfoArgs())


# In[111]:


max(parametrization_helper.get_stretch_angle_offset_from_pattern_params(scale_factors_splines, rparam.getPatternParams().reshape(1, -1)))


# In[112]:


rparam.getPatternParams().shape


# In[113]:


max(rparam.getPatternParams())


# In[114]:


rparam.get_stretch_angle_offset_from_pattern_params(rparam.getPatternParams().reshape(1, -1))


# ## Upsampling and channel generation

# In[120]:


nsubdiv=2
upsampledMesh, upsampledAngles, upsampledPatternParams = rparam.upsampledVertexLeftStretchAnglesAndPatternParameters(nsubdiv)
upsampleMesh_vertices = upsampledMesh.vertices()
upsampleMesh_triangles = upsampledMesh.triangles()
angle_offset = parametrization_helper.get_stretch_angle_offset_from_pattern_params(scale_factors_splines, np.array(upsampledPatternParams).reshape((1, -1)))
upsampledAngles += angle_offset


# In[121]:


# min(angle_offset), max(angle_offset)


# In[122]:


np.save("upsampleMesh_vertices.npy", upsampledMesh.vertices())
np.save("upsampleMesh_triangles.npy", upsampledMesh.triangles())
np.save("upsampleAngles.npy", upsampledAngles)
np.save("upsampledPatternParams.npy", upsampledPatternParams)


# In[123]:


upsampleMesh_vertices = np.load("upsampleMesh_vertices.npy")
upsampleMesh_triangles = np.load("upsampleMesh_triangles.npy")
upsampleAngles = np.load("upsampleAngles.npy")
upsampledPatternParams = np.load("upsampledPatternParams.npy")


# In[124]:


np.set_printoptions(suppress=True)


# In[125]:


max(upsampledPatternParams[0])


# In[126]:


amp_data =  upsampledPatternParams[0]


# In[127]:


import igl


# In[128]:


from parametrization_helper import get_distance_to_line_segments


# In[129]:


def fusing_curve_polyline(patternParams):
#     Draw cosine curves.
    amp = patternParams[0]
    def get_y_from_x(x):
        return amp * np.cos(x) * 0.5 * np.pi + np.pi / 2
    
    x_coords = np.linspace(-np.pi, np.pi, 10)
    y_coords = get_y_from_x(x_coords)
    x_coords += np.pi
    x_coords /= 2
    polyline = np.concatenate(((x_coords).reshape(-1, 1), (y_coords).reshape(-1, 1)), axis = 1)
    return polyline


# In[130]:


def pattern_function(theta, gamma, patternParams, margin, draw_boundary = True):
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


# In[131]:


pattern_function(np.pi / 4, np.pi / 10, [0.], 0.)


# In[132]:


# upsampledPatternParams = np.ones_like(upsampledPatternParams) * 0.4


# In[133]:


import time
start_time = time.time()
(sdfVertices, sdfTris, sdf) = wall_generation.evaluate_cross_field_custom_pattern(upsampleMesh_vertices, upsampleMesh_triangles, upsampleAngles, upsampledPatternParams, pattern_function, frequency=0.1, margin = 0.07)
print(time.time() - start_time)

# pickle.dump((sdfVertices, sdfTris, sdf), open('stripe_sdf_ns4_f100.pkl', 'wb'))

# import pickle, mesh, wall_generation, visualization, numpy as np
# (sdfVertices, sdfTris, sdf) = pickle.load(open('stripe_sdf_ns4_f100.pkl', 'rb'))


# In[134]:


importlib.reload(visualization)

import matplotlib as mpl


# In[135]:


visualization.scalarFieldPlotFast(sdfVertices, sdfTris, sdf, width = 10, height=10)


# In[136]:


importlib.reload(visualization)
visualization.scalarFieldPlotZeroContourFast(sdfVertices, sdfTris, sdf, width = 15, height=10, cmap = mpl.colormaps["PiYG"])


# In[219]:


pts, edges = wall_generation.extract_contours(sdfVertices, sdfTris, sdf,
                                              targetEdgeSpacing=0.4,
                                              minContourLen=1)

visualization.plot_line_segments(pts, edges, width=15, height=15)


# ## Meshing and inflation simulation

# In[134]:


import sheet_meshing, inflation


# In[135]:


m, iwv, iwbv = sheet_meshing.newMeshingAlgorithm(sdfVertices, sdfTris, sdf, pts, edges, triArea=1e0)


# In[136]:


visualization.plot_2d_mesh(m, pointList=np.where(np.array(iwv) == 1)[0], width=100, height=100)


# In[137]:


import inflation
isheet = inflation.InflatableSheet(m, np.array(iwv) != 0)


# In[138]:


from mesh_utilities import SurfaceSampler, tubeRemesh


paramSampler = SurfaceSampler(np.pad(rparam.uv(), [(0, 0), (0, 1)], 'constant'), target_surf.triangles())
liftedSheetPositions = paramSampler.sample(m.vertices(), target_surf.vertices())

isheet.setUninflatedDeformation(liftedSheetPositions.transpose())

isheet.getVars()


# In[139]:


import py_newton_optimizer
niter = 2000
iterations_per_output = 10
opts = py_newton_optimizer.NewtonOptimizerOptions()
opts.useIdentityMetric = True
opts.beta = 1e-4
opts.gradTol = 1e-10
opts.niter = iterations_per_output


# In[140]:


# Are the flat region causing a problem? They might not actually control the metric...
# Try replacing them with single wall...
# Analyze the actual stretching factor (much easier to do with skeleton walls)

from tri_mesh_viewer import TriMeshViewer
viewer = TriMeshViewer(isheet, width=768, height=640)
viewer.showWireframe(True)

viewer.show()


# In[141]:


viewer.update(scalarField=utils.getStrains(isheet)[:, 0])    


# In[142]:


isheet.setUseTensionFieldEnergy(True)

isheet.setUseHessianProjectedEnergy(False)

fixedVars, hessianShift = [], 1e-6

framerate = 20
def cb(it):
    if it % framerate == 0:
        viewer.update(scalarField=utils.getStrains(isheet)[:, 0])    


# ### First solve with low pressure to get out of indefinite state

# In[ ]:


isheet.pressure = 1e-5


# In[ ]:


opts.niter = 5

import time
cr = inflation.inflation_newton(isheet, fixedVars, opts, hessianShift = hessianShift, callback = cb)


# ### Then inflate

# In[ ]:


isheet.pressure = 1e-2


# In[ ]:


opts.niter = 2000
opts.gradTol = 1e-7

import time
benchmark.reset()
cr = inflation.inflation_newton(isheet, fixedVars, opts, hessianShift = hessianShift, callback = cb)
benchmark.report()

isheet.tensionStateHistogram()


# In[ ]:


# Plot maximum tensile strains in the sheet to verify the pressure is reasonable
from matplotlib import pyplot as plt
plt.hist(utils.getStrains(isheet)[:, 0], bins=1000);
plt.xlim(-0.04, 0.1);


# In[ ]:


import gzip


# In[1014]:


pickle.dump(isheet,  gzip.open("igloo_pattern_optimized_2023_12_14_high_frequency.pkl.gz", 'wb'))


# ### Generate Fabrication Files

# In[1112]:


scaleFactor = 1 # Factor for fine-tuning size to fit the machine's build area
channelMargin = 0 / scaleFactor # 8mm channel margin
tabMargin = 2 / scaleFactor # 2mm tab margin


# In[1113]:


import inflation

targetSurf = target_surf
iwv = [isheet.isWallVtx(i) for i in range(isheet.mesh().numVertices())]


# In[1114]:


isheet.mesh()


# In[1115]:


mesh = isheet.mesh()


# In[1116]:


mesh.save("igloo_2D.obj")


# In[1117]:


np.save("igloo_is_wall.npy", iwv)


# In[ ]:





# In[1118]:


uv = rparam.uv()


# In[1119]:


# !pip install shapely==1.7.0


# In[1120]:


import shapely


# In[1121]:


import fabrication


# In[1122]:


importlib.reload(utils)
importlib.reload(fabrication)
importlib.reload(shapely)
import shapely.geometry as shp


# In[1123]:


shapely.__version__


# In[1124]:


# new_fabrication.writeFabricationData('fabrication_data/igloo/free_bdry', isheet.mesh(), iwv, targetSurf, uv,
#                                  scale=scaleFactor,
#                                  channelMargin=channelMargin, fuseSeamWidth=None,
#                                  overlap=0.0, smartOuterChannel=True)


# In[1125]:


import fabrication
fabrication.writeFabricationData('fabrication_data/igloo/free_bdry_low_frequency', isheet.mesh(), isheet.mesh(), iwv, targetSurf, uv,
                                 scale=scaleFactor, numTabs=0, inletOffset=0, tabOffset=0.60 / 80,
                                 channelMargin=channelMargin, tabMargin=tabMargin, tabWidth=5, tabHeight=8, fuseSeamWidth=None, inletScale=None,
                                 overlap=0.0, smartOuterChannel=False)


# In[1078]:


# import fabrication
# fabrication.writeFabricationData('fabrication_data/igloo/free_bdry', isheet.mesh(), isheet.mesh(), iwv, targetSurf, uv,
#                                  scale=scaleFactor, numTabs=0, inletOffset=0.742, tabOffset=0.60 / 80,
#                                  channelMargin=channelMargin, tabMargin=tabMargin, tabWidth=5, tabHeight=8, fuseSeamWidth=1.0, inletScale=12 / channelMargin / scaleFactor,
#                                  overlap=0.0, smartOuterChannel=True)


# In[1079]:


from IPython.display import SVG, display


# In[1126]:


display(SVG(filename = "fabrication_data/igloo/free_bdry_low_frequency/orig.wall_boundaries.svg"))


# In[ ]:





# In[ ]:




