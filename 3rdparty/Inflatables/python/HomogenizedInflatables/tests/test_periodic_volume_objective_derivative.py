import sys
sys.path.append('../..')
import inflation, numpy as np, importlib, fd_validation, visualization, parametric_pillows, wall_generation
from numpy.linalg import norm
import MeshFEM, parallelism, benchmark, utils
import periodic_unit_helper
import numpy.linalg as la

from test_constructors import *
import fd_validation
import pytest

# The square pillow with 5 vertices can not be used for this test because it has zero boundary volume.
constructors = [construct_two_disconnected_square, construct_two_squares_connected_at_center]

@pytest.mark.parametrize("constructor", constructors)
def test_periodic_boundary_volume(constructor):
    ipu = constructor()
    ipu.setVars(np.random.uniform(-1e-3, 1e-3, ipu.numVars()) + ipu.getVars())
    sheet = ipu.sheet
    mesh = ipu.sheet.mesh()
    def get_triangles_of_boundary_edges(edge):
        tri1 = [sheet.getDeformedVtxPosition(edge[0], 0), sheet.getDeformedVtxPosition(edge[1], 0), sheet.getDeformedVtxPosition(edge[1], 1)]
        tri2 = [sheet.getDeformedVtxPosition(edge[0], 0), sheet.getDeformedVtxPosition(edge[1], 1), sheet.getDeformedVtxPosition(edge[0], 1)]
        return [tri1, tri2]
    tris = []
    for edge in mesh.boundaryElements():
        tris += get_triangles_of_boundary_edges(edge)

    volume = 0
    for tri in tris:
        volume += la.det(tri)
    assert(np.abs(volume / 6 - ipu.periodicVolume()) < 1e-8)

class fd_wrapper:
    def __init__(self, ipu):
        self.ipu = ipu

    def setVars(self, v):
        self.ipu.sheet.setVars(v)
    def numVars(self):
        return self.ipu.sheet.numVars()

    def getVars(self):
        return self.ipu.sheet.getVars()

    def energy(self):   return self.ipu.energyPeriodicPressurePotential()
    def gradient(self): return self.ipu.gradientPeriodicPressurePotential()    
    def hessian(self): return self.ipu.hessianPeriodicPressurePotential()

# Old tests for using the boundary integral for computing the periodic volumes. 

# @pytest.mark.parametrize("constructor", constructors)
# def test_periodic_boundary_volume_gradient(constructor):
#     ipu = constructor()
#     ipu.setVars(np.random.uniform(-1e-3, 1e-3, ipu.numVars()) + ipu.getVars())
#     ipu.sheet.pressure = 1
#     (_, errors, _) = fd_validation.gradConvergence(fd_wrapper(ipu))
#     assert(np.min(errors) < 1e-6)

# @pytest.mark.parametrize("constructor", constructors)
# def test_periodic_boundary_volume_hessian(constructor):
#     ipu = constructor()
#     ipu.setVars(np.random.uniform(-1e-3, 1e-3, ipu.numVars()) + ipu.getVars())
#     ipu.sheet.pressure = 1
#     (_, errors, _) = fd_validation.hessConvergence(fd_wrapper(ipu))
#     assert(np.min(errors) < 1e-6)