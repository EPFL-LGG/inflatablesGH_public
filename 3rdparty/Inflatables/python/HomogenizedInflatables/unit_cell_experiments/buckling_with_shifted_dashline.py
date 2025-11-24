#!/usr/bin/env python
# coding: utf-8

# In[1]:


import sys; sys.path.append('../..'); sys.path.append('../../..')


# In[2]:


from periodic_simulation_setup import *


# In[3]:


h = 3
w = 3
res = 100
triArea = h * w / res
avg_len = 0.4

shift = np.array([0., 0.])

ipu, points, segment_edges, m, marker= periodic_unit_helper.get_shifted_dashline(h, w, avg_len, shift, False)

# visualization.plot_line_segments(points, segment_edges)

finalMarkers = np.where(np.array(marker) == 1)[0]

m, finalMarkers = periodic_unit_helper.shift_and_merge_2D_periodic_mesh(m, finalMarkers)
m, finalMarkers = periodic_unit_helper.shift_and_merge_2D_periodic_mesh(m, finalMarkers, axis = 1)

visualization.plot_2d_mesh(m, pointList=finalMarkers, width=5, height=5)


# In[4]:


fusedVtx = get_fusedVtx_using_markers(len(m.vertices()), finalMarkers)


# In[5]:


fuse_boundary = False


# In[6]:


# Fuse boundary
if fuse_boundary:
    bbox = igl.bounding_box(m.vertices())

    max_x = max(bbox[0][:, 0])
    min_x = min(bbox[0][:, 0])
    max_y = max(bbox[0][:, 1])
    min_y = min(bbox[0][:, 1])

    vxs = m.vertices()
    for i, vx in enumerate(m.vertices()):
        if np.abs(vx[0] - max_x) < 1e-6:
            fusedVtx[i] = True
        if np.abs(vx[0] - min_x) < 1e-6:
            fusedVtx[i] = True    
        if np.abs(vx[1] - max_y) < 1e-6:
            fusedVtx[i] = True
        if np.abs(vx[1] - min_y) < 1e-6:
            fusedVtx[i] = True    


# In[7]:


ipu = inflation.InflatablePeriodicUnit(m, fusedVtx = fusedVtx, epsilon = 1e-9)


# In[8]:


viewer = TriMeshViewer(ipu, width=768, height=640)
viewer.showWireframe(True)
viewer.show()


# In[9]:


# Choose strategy for constraining rigid motion
fixedVars, hessianShift = periodic_unit_helper.get_center_fixedVars(ipu), 0
fixedVars, hessianShift = [], 1e-6


# In[10]:


ipu.sheet.setUseTensionFieldEnergy(True)
ipu.sheet.setUseHessianProjectedEnergy(False)
ipu.sheet.disableFusedRegionTensionFieldTheory(False)

ipu.sheet.pressure = 3


# In[11]:


np.set_printoptions(precision=4, suppress=True)


# In[12]:


benchmark.reset()

opts.niter = 100
framerate = 5 # Update every 5 iterations
def cb(it):
    if it % framerate == 0:
        viewer.update(scalarField=utils.getStrains(ipu.sheet)[:, 0])
cr = inflation.inflation_newton(ipu, fixedVars, opts, callback=cb, hessianShift = hessianShift)
benchmark.report()


# In[ ]:


ipu.sheet.energy(inflation.InflatableSheet.EnergyType.Elastic)


# In[ ]:


ipu.energy()


# In[ ]:


ipu.energy(inflation.InflatableSheet.EnergyType.Elastic)


# ### Vibrational Mode analysis
# 

# In[ ]:


lambdas, modes = compute_vibrational_modes.compute_vibrational_modes(ModalAnalysisWrapper(ipu), mtype=compute_vibrational_modes.MassMatrixType.FULL, n=16, sigma=-1e-10, fixedVars = [])

import mode_viewer, importlib
mview = mode_viewer.ModeViewer(ipu, modes, lambdas, amplitude=100)
mview.show()

