# import sys; sys.path.append("../experiments")
# import experiment_helper
# import igl
# from periodic_simulation_setup import *
# import json

# allowBending = False
# useTFT = True
# useMirror = True

# name = "three_star"
# time_stamp = time.strftime("%Y_%m_%d_%H_%M")
# result_folder = 'output/{}/{}'.format(name, time_stamp)
# if not os.path.exists(result_folder):
#     os.makedirs(result_folder)  

# # pressure = 0.8
# stiffness_pressure = 0.4
# scale_factor_pressure = 0.01

# # radius = np.linspace(1, 2.4, 15)
# # angles = np.linspace(0, 15, 16)

# radius = np.linspace(1, 2.4, 8)
# angles = np.linspace(0, 15, 6)


# for r in radius:
#     for angle in angles:
#         h = 5
#         avg_len = 0.1
#         dash_point = np.array([np.cos(angle / 180 * np.pi), np.sin(angle / 180 * np.pi)]) * r + np.array([0, 0])
#         ipu, m, marker = pattern_generator_using_gmsh.get_three_star(h, avg_len, avg_len, dash_point = dash_point)        

#         finalMarkers = np.where(np.array(marker) == 1)[0]
#         m, finalMarkers = periodic_unit_helper.shift_and_merge_2D_periodic_mesh(m, finalMarkers, flip_orientation= 0)
#         m, finalMarkers = periodic_unit_helper.shift_and_merge_2D_periodic_mesh(m, finalMarkers, axis = 1, flip_orientation= 1)

#         fusedVtx = get_fusedVtx_using_markers(len(m.vertices()), finalMarkers)

#         ipu = inflation.InflatablePeriodicUnit(m, fusedVtx = fusedVtx, epsilon = 1e-9)

#         az_ipu = get_az_ipu_from_ipu(ipu, m, fusedVtx, useTFT, False)
        
#         dof = np.load("../experiments/output/three_star/2023_04_26_19_43/high_pressure_dofs_three_star_{:.1f}_{:.2f}.npy".format(r, angle))
#         az_ipu.setVars(dof)

#         viewer = TriMeshViewer(az_ipu, width=768, height=640)
#         viewer.showWireframe(True)

#         viewer.update(scalarField=utils.getStrains(az_ipu.ipu.sheet)[:, 0])


#         render = viewer.offscreenRenderer(1000, 1000)
#         render.render()
#         render.save("{}/render_{}_{}.png".format(result_folder, name, "{:.1f}_{:.2f}".format(r, angle)))


# import sys; sys.path.append("../experiments")
# import experiment_helper
# import igl
# from periodic_simulation_setup import *
# import json

# import parallelism, multiprocessing, itertools, setproctitle
# import os, time, numpy as np

# allowBending = False
# useTFT = True
# useMirror = True

# name = "grid_dots"
# time_stamp = time.strftime("%Y_%m_%d_%H_%M")
# result_folder = 'output/{}/{}'.format(name, time_stamp)
# if not os.path.exists(result_folder):
#     os.makedirs(result_folder)  

# # pressure = 0.8
# stiffness_pressure = 0.4
# scale_factor_pressure = 0.01

# # radius = np.linspace(1, 2.4, 15)
# # angles = np.linspace(0, 15, 16)

# radius = np.linspace(1, 2.4, 8)
# angles = np.linspace(0, 15, 6)

# def render_image(i):
#     m = MeshFEM.Mesh('../experiments/output/grid_dots/2023_05_02_10_19/mesh_grid_dots_{}.obj'.format(i))
#     vx = m.vertices()
#     new_vx = np.zeros((vx.shape[0], 3))
#     new_vx[:, :2] = vx
#     m = MeshFEM.Mesh(new_vx, m.elements())
#     fusedVtx = np.load('../experiments/output/grid_dots/2023_05_02_10_19/fusedVtx_grid_dots_{}.npy'.format(i))
#     ipu = inflation.InflatablePeriodicUnit(m, fusedVtx = fusedVtx, epsilon = 1e-9)

#     az_ipu = get_az_ipu_from_ipu(ipu, m, fusedVtx, useTFT, False)
        
#     dof = np.load("../experiments/output/grid_dots/2023_05_02_10_19/high_pressure_dofs_grid_dots_{}.npy".format(i))
#     az_ipu.setVars(dof)

#     viewer = TriMeshViewer(az_ipu, width=768, height=640)
#     viewer.showWireframe(True)

#     viewer.update(scalarField=utils.getStrains(az_ipu.ipu.sheet)[:, 0])


#     render = viewer.offscreenRenderer(1000, 1000)
#     render.render()
#     render.save("{}/render_{}_{}.png".format(result_folder, name, i))


# def mute():
#     "Suppress stdout output"
#     sys.stdout = open(os.devnull, 'w')
    

# num_thread = 1
# # For embarassingly parallel runs, using single thread should be more efficient
# parallelism.set_max_num_tbb_threads(num_thread)

# # Set num logical cores
# logical_core_count = 10

# print("Using {} logical cores".format(logical_core_count))


# import time
# st = time.time()

# with multiprocessing.Pool(logical_core_count, initializer=mute) as p:
#     p.map(render_image, iter(range(300)))
        
# # get the end time
# et = time.time()

# # get the execution time
# elapsed_time = et - st
# print('Execution time:', elapsed_time, 'seconds')



import sys; sys.path.append("../experiments")
import experiment_helper
import igl
from periodic_simulation_setup import *
import json

import parallelism, multiprocessing, itertools, setproctitle
import os, time, numpy as np

allowBending = False
useTFT = True
useMirror = True

name = "mirror_cosine_dash"
time_stamp = time.strftime("%Y_%m_%d_%H_%M")
result_folder = 'output/{}/{}'.format(name, time_stamp)
if not os.path.exists(result_folder):
    os.makedirs(result_folder)  

# pressure = 0.8
stiffness_pressure = 0.4
scale_factor_pressure = 0.01


amplitudes = np.linspace(0, 1, 11)
amplitudes = [amplitudes[0]]
radius = np.linspace(0.5, 2.4, 16)
angles = np.linspace(0, 45, 31)
angles = [angles[-1]]
avg_len = 0.2

def render_image(amp, r, angle):
    m = MeshFEM.Mesh('../experiments/output/mirror_cosine_dash/2023_05_04_18_51/mesh_mirror_cosine_dash_{}.obj'.format("{:.1f}_{:.1f}_{:.2f}".format(amp, r, angle)))
    vx = m.vertices()
    new_vx = np.zeros((vx.shape[0], 3))
    new_vx[:, :2] = vx
    m = MeshFEM.Mesh(new_vx, m.elements())
    fusedVtx = np.load('../experiments/output/mirror_cosine_dash/2023_05_04_18_51/fusedVtx_mirror_cosine_dash_{}.npy'.format("{:.1f}_{:.1f}_{:.2f}".format(amp, r, angle)))
    ipu = inflation.InflatablePeriodicUnit(m, fusedVtx = fusedVtx, epsilon = 1e-9)

    az_ipu = get_az_ipu_from_ipu(ipu, m, fusedVtx, useTFT, False)
        
    dof = np.load("../experiments/output/mirror_cosine_dash/2023_05_04_18_51/high_pressure_dofs_mirror_cosine_dash_{}.npy".format("{:.1f}_{:.1f}_{:.2f}".format(amp, r, angle)))
    az_ipu.setVars(dof)
    viewer = TriMeshViewer(az_ipu, width=768, height=640)
    viewer.showWireframe(True)
    viewer.update(scalarField=utils.getStrains(az_ipu.ipu.sheet)[:, 0])
    render = viewer.offscreenRenderer(1000, 1000)
    render.render()
    render.save("{}/render_{}_{}.png".format(result_folder, name, "{:.1f}_{:.1f}_{:.2f}".format(amp, r, angle)))


def mute():
    "Suppress stdout output"
    sys.stdout = open(os.devnull, 'w')
    

num_thread = 1
# For embarassingly parallel runs, using single thread should be more efficient
parallelism.set_max_num_tbb_threads(num_thread)

# Set num logical cores
logical_core_count = 10

print("Using {} logical cores".format(logical_core_count))
arg1 = amplitudes
arg2 = radius
arg3 = angles

args = itertools.product(arg1, arg2, arg3)


import time
st = time.time()

with multiprocessing.Pool(logical_core_count, initializer=mute) as p:
    p.starmap(render_image, args)
        
# get the end time
et = time.time()

# get the execution time
elapsed_time = et - st
print('Execution time:', elapsed_time, 'seconds')