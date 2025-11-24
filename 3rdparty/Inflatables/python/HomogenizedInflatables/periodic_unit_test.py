import sys; sys.path.append('..')
import inflation, numpy as np, importlib, fd_validation, visualization, parametric_pillows, wall_generation
from numpy.linalg import norm
import MeshFEM, parallelism, benchmark, utils

import matplotlib.pyplot as plt

import periodic_unit_helper

m, fuseMarkers, brdyWallMarkers = periodic_unit_helper.get_mesh_input()

isheet = inflation.InflatableSheet(m, np.array(fuseMarkers) != 0)

ipu = inflation.InflatablePeriodicUnit(m, np.array(fuseMarkers) != 0)


ipu.getVars()
ipu.setVars(np.ones(ipu.numVars()))