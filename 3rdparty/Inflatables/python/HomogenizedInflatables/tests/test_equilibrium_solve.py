import sys
sys.path.append('../..')
import inflation, numpy as np, importlib, fd_validation, visualization, parametric_pillows, wall_generation
from numpy.linalg import norm
import MeshFEM, parallelism, benchmark, utils
import periodic_unit_helper
import numpy.linalg as la

from test_constructors import *
import fd_validation

def test_equilibrium_solve():
    import py_newton_optimizer
    opts = py_newton_optimizer.NewtonOptimizerOptions()
    opts.useIdentityMetric = True
    opts.beta = 1e-4
    opts.gradTol = 1e-10

    ipu = construct_two_squares_connected_at_center()
    ipu.sheet.pressure = 1
    ipu.sheet.setUseTensionFieldEnergy(True)
    ipu.sheet.setUseHessianProjectedEnergy(False)

    fixedVars, hessianShift = [], 1e-6
    opts.niter = 1000
    cr = inflation.inflation_newton(ipu, [], opts, callback=None, hessianShift = hessianShift)
    assert(cr.success)