import sys; sys.path.append('..')
import inflation, numpy as np, importlib, fd_validation, visualization, parametric_pillows, wall_generation
from numpy.linalg import norm
import MeshFEM, parallelism, benchmark, utils
import periodic_unit_helper
import numpy.linalg as la

Pressure = inflation.InflatableSheet.EnergyType.Pressure
Elastic  = inflation.InflatableSheet.EnergyType.Elastic
Full     = inflation.InflatableSheet.EnergyType.Full


n_vx = [[0, 0], [0, 1], [0, 2],
        [1, 0], [1, 1], [1, 2],
        [2, 0], [2, 1], [2, 2]]
n_edge = [(0, 1), (1, 2), 
          (3, 4), (4, 5),
          (6, 7), (7, 8),
          (0, 3), (3, 6),
          (1, 4), (4, 7),
          (2, 5), (5, 8)]
triArea = 0.5


m, fuseMarkers, fuseSegments = wall_generation.triangulate_channel_walls(n_vx, n_edge, triArea, flags="Y")

fuseMarkers = [0] * 9

fuseMarkers[4] = 1

ipu = inflation.InflatablePeriodicUnit(m, np.array(fuseMarkers) == 1, epsilon = 1e-5)

ipu.setVars(ipu.getVars())

ipu.energy()
ipu.hessianSparsityPattern()

ipu.hessian(energyType = Elastic)

# H = ipu.bent_sheet_hessian()
