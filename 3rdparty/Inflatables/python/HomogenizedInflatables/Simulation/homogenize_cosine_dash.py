#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import sys; sys.path.append(os.path.join(os.path.abspath(''), '../experiments'));

import experiment_helper
import igl
from periodic_simulation_setup import *
import json

import parallelism, multiprocessing, itertools, setproctitle
import os, time, numpy as np

import json


# In[2]:


allowBending = False
useTFT = True
disableFusedRegionTFT = False
stiffness_pressure = 0.4
scale_factor_pressure = 0.01
avg_len = 1


# In[3]:


tag = '0.30_1.50_60.00'
amp = float(tag.split('_')[0])
r = float(tag.split('_')[1])
angle = float(tag.split('_')[2])


# In[4]:


h = 5
dash_point = np.array([np.cos(angle / 180 * np.pi), np.sin(angle / 180 * np.pi)]) * r + np.array([0, 0])
ipu, m, marker = pattern_generator_using_gmsh.get_cosine_dash(h, avg_len, avg_len, amplitude=amp, dash_point = dash_point)

finalMarkers = np.where(np.array(marker) == 1)[0]
m, finalMarkers = periodic_unit_helper.shift_and_merge_2D_periodic_mesh(m, finalMarkers, flip_orientation= 0)
m, finalMarkers = periodic_unit_helper.shift_and_merge_2D_periodic_mesh(m, finalMarkers, axis = 1, flip_orientation= 1)

fusedVtx = get_fusedVtx_using_markers(len(m.vertices()), finalMarkers)


# In[5]:


ipu = inflation.InflatablePeriodicUnit(m, fusedVtx = fusedVtx, epsilon = 1e-9)


# In[6]:


max(m.edgeLengths()), min(m.edgeLengths())


# In[7]:


visualization.plot_2d_mesh(m, pointList=finalMarkers, width=10, height=10)


# In[8]:


viewer = TriMeshViewer(ipu, width=768, height=640)
viewer.showWireframe(True)
viewer.show()


# In[9]:


configure_solver_parallelism()


# In[10]:


framerate = 5 # Update every 5 iterations
def cb(it):
    if it % framerate == 0:
        viewer.update(scalarField=utils.getStrains(ipu.sheet)[:, 0])


# In[11]:


hessianShiftForRigidMotion = 1e-10
hessianShiftForAlphainPlanar = 1e-12


# In[12]:


ipu.sheet.setUseTensionFieldEnergy(useTFT)
ipu.sheet.setUseHessianProjectedEnergy(False)
if (disableFusedRegionTFT):
    ipu.sheet.disableFusedRegionTensionFieldTheory(False)
ipu.sheet.pressure = stiffness_pressure


# In[13]:


fixedVars, hessianShift = list(periodic_unit_helper.get_center_fixedVars(ipu)), hessianShiftForAlphainPlanar


# In[14]:


opts.niter = 500
opts.gradTol = 1e-10

cr = inflation.inflation_newton(ipu, fixedVars, opts, callback=cb, hessianShift = hessianShift)


# In[15]:


opts.niter = 500
opts.gradTol = 1e-10

cr = inflation.inflation_newton(ipu, fixedVars, opts, callback=cb, hessianShift = hessianShift)


# In[16]:


experiment_log = {}
experiment_log["Ipu simulation succeed"] = int(cr.success)
experiment_log["Simulation Kappa value"] = (ipu.getVars()[-2])
if np.abs(ipu.getVars()[-2]) > 1e-6:
    experiment_log["Planar equilibrium"] = 0
    print("Warning: Can not compute stiffness due to non-planar equilibrium!")
else:
    experiment_log["Planar equilibrium"] = 1


# In[17]:


ipu.getVars()[:3]


# In[18]:


az_ipu = get_az_ipu_from_ipu(ipu, m, fusedVtx, useTFT, disableFusedRegionTFT)


# In[19]:


az_viewer = TriMeshViewer(az_ipu, width=768, height=640)
az_viewer.showWireframe(True)
az_viewer.show()


# In[20]:


framerate = 5 # Update every 5 iterations
def az_cb(it):
    if it % framerate == 0:
        az_viewer.update(scalarField=utils.getStrains(az_ipu.ipu.sheet)[:, 0])


# In[21]:


az_ipu.getVars()[az_ipu.get_average_z_idx()]


# In[22]:


az_ipu.getVars()[-2:]


# In[23]:


periodic_unit_helper.getNumpyArrayFromCSC(az_ipu.hessian())[-2:, -2:]


# In[24]:


vars = az_ipu.getVars()
vars[-2] = 0
vars[-1] += 0.1
az_ipu.setVars(vars)
print(az_ipu.gradient()[-2:])
periodic_unit_helper.getNumpyArrayFromCSC(az_ipu.hessian())[-2:, -2:]


# In[25]:


import fd_validation


# In[26]:


F_indices = [np.arange(0, 3)]
F_indices = np.array(F_indices).flatten()

u_indices = [np.arange(3, az_ipu.numVars() - 2)]
u_indices = np.array(u_indices).flatten()


R_star_indices = [np.arange(3 + az_ipu.ipu.numFluctuationDisplacementVars(), az_ipu.numVars())]
R_star_indices = np.array(R_star_indices).flatten()




# In[27]:


var_types = ['F', 'u', 'R']
var_indices = {'F': F_indices,
               'u': u_indices, 
               'R': R_star_indices}


# In[28]:


az_ipu.ipu.sheet.usingTensionFieldEnergy(False)


# In[29]:


# fd_validation.hessConvergencePlot(az_ipu.ipu.sheet)


# In[30]:


# fd_validation.hessian_convergence_block_plot(az_ipu, var_types, var_indices, customArgs={'energyType': ipu.EnergyType.Pressure})


# In[31]:


# import mode_viewer

# import compute_vibrational_modes
# class ModalAnalysisWrapper:
#     def __init__(self, sheet):
#         self.sheet = sheet
#     def hessian(self):
#         return self.sheet.hessian(inflation.InflatablePeriodicUnit.EnergyType.Elastic)


# In[32]:


# lambdas, modes = compute_vibrational_modes.compute_vibrational_modes(ModalAnalysisWrapper(az_ipu), mtype=compute_vibrational_modes.MassMatrixType.FULL, n=10, sigma=-1e-10, fixedVars = [az_ipu.get_average_z_idx(), az_ipu.numVars() - 1])


# In[33]:


# import mode_viewer, importlib
# mview = mode_viewer.ModeViewer(az_ipu, modes, lambdas, amplitude=200)
# mview.show()


# In[34]:


# if not allowBending:
#     fixedVars, hessianShift = [az_ipu.numVars() - 2, az_ipu.numVars() - 1], 1e-6
# else:
fixedVars, hessianShift = [az_ipu.get_average_z_idx(), periodic_unit_helper.get_center_fixedVars(ipu)[0], periodic_unit_helper.get_center_fixedVars(ipu)[1]], hessianShiftForAlphainPlanar
opts.niter = 1000
opts.gradTol = 1e-10
cr = inflation.inflation_newton(az_ipu, fixedVars, opts, callback=az_cb, hessianShift = hessianShift)


# In[35]:


def get_label(amp, r, angle):
    return "{:.2f}_{:.2f}_{:.2f}".format(amp, r, angle)


# In[36]:


name = "mirror_cosine_dash"
time_stamp = time.strftime("%Y_%m_%d_%H_%M")
# time_stamp = "2023_10_25_14_54"
base_folder = 'output/{}/{}'.format(name, time_stamp)
    
result_folder = "{}/{}".format(base_folder, get_label(amp, r, angle))
if not os.path.exists(result_folder):
    os.makedirs(result_folder)  
    
render_images = True
variable = get_label(amp, r, angle)


# In[37]:


stiffness_fixedVars = [az_ipu.get_average_z_idx(), periodic_unit_helper.get_center_fixedVars(ipu)[0], periodic_unit_helper.get_center_fixedVars(ipu)[1], az_ipu.numVars() - 1, az_ipu.numVars() - 2]


# In[38]:


az_optimizer = inflation.get_inflation_optimizer(az_ipu, stiffness_fixedVars, opts, callback=az_cb, hessianShift = 0)


# In[39]:


stiffness_values_14, sampled_alphas, stiffness_coefficient = visualize_sampled_bending_stiffness(az_ipu, 1000, az_optimizer, hessianShift = 1e-14, fixedVars = [az_ipu.get_average_z_idx()], filename = "{}/stiffness_{}_{}.png".format(result_folder, name, variable), generate_images = render_images)
stiffness_values_10, sampled_alphas, stiffness_coefficient = visualize_sampled_bending_stiffness(az_ipu, 1000, az_optimizer, hessianShift = 1e-10, fixedVars = [az_ipu.get_average_z_idx()], filename = "{}/stiffness_{}_{}.png".format(result_folder, name, variable), generate_images = render_images)
stiffness_values_6, sampled_alphas, stiffness_coefficient = visualize_sampled_bending_stiffness(az_ipu, 1000, az_optimizer, hessianShift = 1e-6, fixedVars = [az_ipu.get_average_z_idx()], filename = "{}/stiffness_{}_{}.png".format(result_folder, name, variable), generate_images = render_images)


# In[40]:


stiffness_values, sampled_alphas, stiffness_coefficient = visualize_sampled_bending_stiffness(az_ipu, 1000, az_optimizer, hessianShift = 0, fixedVars = stiffness_fixedVars, filename = "{}/stiffness_{}_{}.png".format(result_folder, name, variable), generate_images = render_images)
np.save("{}/stiffness_values_{}_{}.npy".format(result_folder, name, variable), stiffness_values)
np.save("{}/sampled_alphas_{}_{}.npy".format(result_folder, name, variable), sampled_alphas)
np.save("{}/stiffness_coefficient_{}_{}.npy".format(result_folder, name, variable), stiffness_coefficient)


# In[ ]:


la.norm(stiffness_values - stiffness_values_6), la.norm(stiffness_values - stiffness_values_10), la.norm(stiffness_values - stiffness_values_14)


# In[ ]:


from IPython.display import Image
Image(filename="{}/stiffness_{}_{}.png".format(result_folder, name, variable)) 


# In[35]:


np.save("{}/scale_factors_{}_{}.npy".format(result_folder, name, variable), get_deformation_scale_factors(az_ipu.ipu))
np.save("{}/average_deformation_gradient_matrix_{}_{}.npy".format(result_folder, name, variable), get_deformation_matrix(az_ipu.ipu))


# In[36]:


min(stiffness_values), max(stiffness_values)


# In[37]:


import periodic_simulation_setup
import importlib
importlib.reload(periodic_simulation_setup)


# In[38]:


bs_obj = periodic_simulation_setup.bending_stiffness_class(az_ipu, az_ipu.ipu.sheet, az_optimizer, az_viewer, fixedVars = [])

bs_obj.setVars(bs_obj.getVars())

fd_validation.secondDerivativeConvergencePlot(bs_obj, epsilons = np.logspace(-6, 1, 50))


# In[ ]:




