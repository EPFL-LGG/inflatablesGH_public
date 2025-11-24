import sys; sys.path.append('../'); sys.path.append('.'); sys.path.append('experiments/'); sys.path.append('../gmsh')
import inflation, sparse_matrices, mesh, numpy as np, importlib, pickle
import inflatables_parametrization as parametrization
from numpy.linalg import norm
from io_redirection import suppress_stdout
import visualization

import MeshFEM, parallelism, benchmark, utils
parallelism.set_max_num_tbb_threads(32)
parallelism.set_gradient_assembly_num_threads(32)
parallelism.set_hessian_assembly_num_threads(32)

from py_newton_optimizer import NewtonOptimizerOptions
import parametrization_helper

PET = parametrization.RegularizedPatternParametrizer.PatternEnergyType
names = ["Full", "PReg", "Bend", "Fitt", "SReg", "PhiR"]
import os
import numpy.linalg as la

def optimize_rparam(param, patternRegW, phiRegW, bendRegW = 0.0, update_uv = True, niter = 100, use_knitro = True):
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
        fixedvars = range(param.patternOffset())

    if use_knitro:
        param.useBarrier = False
        cr = parametrization.pattern_parametrization_knitro(param, opts.niter, fixedvars)
    else:
        cr = parametrization.pattern_parametrization_newton(param, fixedvars, opts)
    benchmark.report()
    return cr

def initialize_pattern_parameters(rparam, lines, num_pattern_params, path = None, use_knitro = True):
    rparam.patternRegW = 0
    rparam.phiRegW = 0
    rparam.bendRegW = 0
    energy_values = list(map(rparam.energy, [PET.Full, PET.PatternRegularization, PET.Bending, PET.DEBUG_Fitting, PET.DEBUG_StretchRegularization, PET.DEBUG_PhiRegularization]))
    print("Before: ", list(zip(names, energy_values)))

    benchmark.reset()
    with suppress_stdout(): report = optimize_rparam(rparam, 0, 0, bendRegW = 0, update_uv = False, niter = 200, use_knitro = use_knitro)
    benchmark.report()

    benchmark.reset()
    with suppress_stdout(): report = optimize_rparam(rparam, 0, 0, bendRegW = 0, update_uv = False, niter = 200, use_knitro = use_knitro)
    benchmark.report()

    benchmark.reset()
    with suppress_stdout(): report = optimize_rparam(rparam, 0, 0, bendRegW = 0, update_uv = False, niter = 200, use_knitro = use_knitro)
    benchmark.report()

    energy_values = list(map(rparam.energy, [PET.Full, PET.PatternRegularization, PET.Bending, PET.DEBUG_Fitting, PET.DEBUG_StretchRegularization, PET.DEBUG_PhiRegularization]))
    print("After: ", list(zip(names, energy_values)))

    parametrization_helper.visualize_scale_factors(lines, rparam.getAlphas(), rparam.getBetas())

    importlib.reload(visualization)
    visualization.visualize_both(rparam, height = 4)
    visualization.visualize_pattern(rparam, height = 4, num_pattern_vars=num_pattern_params, path = path)
    return energy_values

def add_phi_regularization(rparam, lines, num_pattern_params, path = None, use_knitro = True):
    rparam.phiRegW = 1
    phiRegW = rparam.phiRegW / la.norm(rparam.gradient(PET.DEBUG_PhiRegularization)) * la.norm(rparam.gradient(PET.DEBUG_Fitting)) * 4

    rparam.phiRegW = phiRegW
    energy_values = list(map(rparam.energy, [PET.Full, PET.PatternRegularization, PET.Bending, PET.DEBUG_Fitting, PET.DEBUG_StretchRegularization, PET.DEBUG_PhiRegularization]))
    print("Before: ", list(zip(names, energy_values)))

    benchmark.reset()
    with suppress_stdout(): report = optimize_rparam(rparam, patternRegW = 0, phiRegW = phiRegW, bendRegW = 0, update_uv = True, niter = 200, use_knitro = use_knitro)
    benchmark.report()

    energy_values = list(map(rparam.energy, [PET.Full, PET.PatternRegularization, PET.Bending, PET.DEBUG_Fitting, PET.DEBUG_StretchRegularization, PET.DEBUG_PhiRegularization]))
    print("After: ", list(zip(names, energy_values)))

    visualization.visualize_both(rparam, height = 4)
    visualization.visualize_pattern(rparam, height = 4, num_pattern_vars=num_pattern_params, width = 30, path = path)

    parametrization_helper.visualize_scale_factors(lines, rparam.getAlphas(), rparam.getBetas())
    return phiRegW, energy_values

def add_pattern_regularization(rparam, lines, num_pattern_params, phiRegW, path = None, use_knitro = True):
    rparam.patternRegW = 1
    patternRegW = rparam.patternRegW / la.norm(rparam.gradient(PET.PatternRegularization)) * la.norm(rparam.gradient(PET.DEBUG_Fitting)) * 4
    rparam.patternRegW = patternRegW

    energy_values = list(map(rparam.energy, [PET.Full, PET.PatternRegularization, PET.Bending, PET.DEBUG_Fitting, PET.DEBUG_StretchRegularization, PET.DEBUG_PhiRegularization]))
    print("Before: ", list(zip(names, energy_values)))

    benchmark.reset()
    with suppress_stdout(): report = optimize_rparam(rparam, patternRegW = patternRegW, phiRegW = phiRegW, bendRegW = 0, update_uv = True, niter = 100, use_knitro = use_knitro)
    benchmark.report()

    benchmark.reset()
    with suppress_stdout(): report = optimize_rparam(rparam, patternRegW = patternRegW, phiRegW = phiRegW / 2, bendRegW = 0, update_uv = True, niter = 100, use_knitro = use_knitro)
    benchmark.report()

    energy_values = list(map(rparam.energy, [PET.Full, PET.PatternRegularization, PET.Bending, PET.DEBUG_Fitting, PET.DEBUG_StretchRegularization, PET.DEBUG_PhiRegularization]))
    print("After: ", list(zip(names, energy_values)))

    visualization.visualize_both(rparam, height = 4)
    visualization.visualize_pattern(rparam, height = 4, num_pattern_vars=num_pattern_params, width = 30, path = path)
    return patternRegW, energy_values


def add_bending_energy(rparam, lines, num_pattern_params, phiRegW, patternRegW, path = None, use_knitro = True):
    rparam.bendRegW = 1
    bendRegW = rparam.bendRegW / la.norm(rparam.gradient(PET.Bending)) * la.norm(rparam.energy(PET.DEBUG_Fitting)) * 2
    rparam.bendRegW = bendRegW

    energy_values = list(map(rparam.energy, [PET.Full, PET.PatternRegularization, PET.Bending, PET.DEBUG_Fitting, PET.DEBUG_StretchRegularization, PET.DEBUG_PhiRegularization]))
    print("Before: ", list(zip(names, energy_values)))

    benchmark.reset()
    with suppress_stdout(): report = optimize_rparam(rparam, patternRegW = patternRegW / 2, phiRegW = phiRegW / 4, bendRegW = bendRegW, update_uv = False, niter = 400, use_knitro = use_knitro)
    benchmark.report()

    energy_values = list(map(rparam.energy, [PET.Full, PET.PatternRegularization, PET.Bending, PET.DEBUG_Fitting, PET.DEBUG_StretchRegularization, PET.DEBUG_PhiRegularization]))
    print("After fixed uv: ", list(zip(names, energy_values)))

    visualization.visualize_both(rparam, height = 4)
    visualization.visualize_pattern(rparam, height = 4, num_pattern_vars=num_pattern_params, width = 30, path = (path + '_fixed_uv') if path is not None else None)

    benchmark.reset()
    with suppress_stdout(): report = optimize_rparam(rparam, patternRegW = patternRegW / 2, phiRegW = phiRegW / 4, bendRegW = bendRegW, update_uv = True, niter = 200, use_knitro = use_knitro)
    benchmark.report()

    benchmark.reset()
    with suppress_stdout(): report = optimize_rparam(rparam, patternRegW = patternRegW / 2, phiRegW = phiRegW / 8, bendRegW = bendRegW / 2, update_uv = True, niter = 200, use_knitro = use_knitro)
    benchmark.report()

    # benchmark.reset()
    # with suppress_stdout(): report = optimize_rparam(rparam, patternRegW = patternRegW, phiRegW = phiRegW / 8, bendRegW = bendRegW / 4, update_uv = True, niter = 200)
    # benchmark.report()

    energy_values = list(map(rparam.energy, [PET.Full, PET.PatternRegularization, PET.Bending, PET.DEBUG_Fitting, PET.DEBUG_StretchRegularization, PET.DEBUG_PhiRegularization]))
    print("After freed uv: ", list(zip(names, energy_values)))

    visualization.visualize_both(rparam, height = 4)
    visualization.visualize_pattern(rparam, height = 4, num_pattern_vars=num_pattern_params, width = 30, path = (path + '_freed_uv') if path is not None else None)

    parametrization_helper.visualize_scale_factors(lines, rparam.getAlphas(), rparam.getBetas(), path = (path + '_freed_uv_scale_factors') if path is not None else None)

    visualization.visualizeChannelOrientation(rparam, quiver=visualization.QuiverVisualization.PER_VTX, orientationHue=False, width = 5, height = 5)
    visualization.visualizeChannelOrientation(rparam, quiver=visualization.QuiverVisualization.PER_TRI, orientationHue=False, width = 10, height = 10)
    return bendRegW, energy_values