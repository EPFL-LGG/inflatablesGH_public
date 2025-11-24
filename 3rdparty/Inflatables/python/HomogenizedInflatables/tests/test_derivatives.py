import sys
sys.path.append('../..')
sys.path.append('..')
import inflation, numpy as np, importlib, fd_validation, visualization, parametric_pillows, wall_generation
from numpy.linalg import norm
import MeshFEM, parallelism, benchmark, utils
import periodic_unit_helper
import periodic_simulation_setup
import numpy.linalg as la

from test_constructors import *
import fd_validation

epsilon = 1e-5

import pytest

Pressure = inflation.InflatableSheet.EnergyType.Pressure
Elastic  = inflation.InflatableSheet.EnergyType.Elastic
Full     = inflation.InflatableSheet.EnergyType.Full

constructors = [construct_square_pillow_five_vertices, construct_two_disconnected_square, construct_two_squares_connected_at_center, construct_rectangle_connected_at_one_side, construct_dashline]

energyTypes = [Pressure, Elastic, Full]

test_data = [(constructor, energyType) for constructor in constructors for energyType in energyTypes]

@pytest.mark.parametrize("constructor, energyType", test_data)
def test_gradient(constructor, energyType):
    ipu = constructor()
    ipu.sheet.setUseTensionFieldEnergy(False)
    ipu.sheet.setUseHessianProjectedEnergy(False)
    ipu.sheet.pressure = 3
    ipu.setVars(ipu.getVars() + np.random.uniform(-1e-2, 1e-2, ipu.numVars()))
    customArgs = {"energyType": energyType}
    (_, errors, _) = fd_validation.gradConvergence(ipu, customArgs = customArgs)
    assert np.min(errors) < epsilon

@pytest.mark.parametrize("constructor, energyType", test_data)
def test_hessian(constructor, energyType):
    ipu = constructor()
    ipu.sheet.setUseTensionFieldEnergy(False)
    ipu.sheet.setUseHessianProjectedEnergy(False)
    ipu.sheet.pressure = 3
    ipu.setVars(ipu.getVars() + np.random.uniform(-1e-2, 1e-2, ipu.numVars()))
    customArgs = {"energyType": energyType}
    (_, errors, _) = fd_validation.hessConvergence(ipu, customArgs = customArgs)
    assert np.min(errors) < epsilon

@pytest.mark.parametrize("constructor, energyType", test_data)
def test_bending_gradient(constructor, energyType):
    ipu = constructor()
    ipu.sheet.setUseTensionFieldEnergy(False)
    ipu.sheet.setUseHessianProjectedEnergy(False)
    ipu.sheet.pressure = 3
    ipu.setVars(ipu.getVars() + np.random.uniform(-1e-2, 1e-2, ipu.numVars()))
    customArgs = {"energyType": energyType}
    (_, errors, _) = fd_validation.gradConvergence(periodic_simulation_setup.bent_sheet_wrapper(ipu), customArgs = customArgs)
    assert np.min(errors) < epsilon

@pytest.mark.parametrize("constructor, energyType", test_data)
def test_bending_hessian(constructor, energyType):
    ipu = constructor()
    ipu.sheet.setUseTensionFieldEnergy(False)
    ipu.sheet.setUseHessianProjectedEnergy(False)
    ipu.sheet.pressure = 3
    ipu.setVars(ipu.getVars() + np.random.uniform(-1e-2, 1e-2, ipu.numVars()))
    customArgs = {"energyType": energyType}
    (_, errors, _) = fd_validation.hessConvergence(periodic_simulation_setup.bent_sheet_wrapper(ipu), customArgs = customArgs)
    assert np.min(errors) < epsilon