import sys
sys.path.append('..')
import mesh, inflation, numpy as np, importlib, fd_validation, visualization, parametric_pillows, wall_generation, py_newton_optimizer
from numpy.linalg import norm

m, fuseMarkers = wall_generation.triangulate_channel_walls(*parametric_pillows.concentricCircles(8, 50), 0.001)
isheet = inflation.InflatableSheet(m, np.array(fuseMarkers) != 0)
isheet.pressure = 30

targetSurf = mesh.Mesh('data/pringle.obj')
targetAttractedSheet = inflation.TargetAttractedInflation(isheet, targetSurf)

targetAttractedSheet.fittingWeight = 0.1

opts = py_newton_optimizer.NewtonOptimizerOptions()
opts.useIdentityMetric = True
opts.beta = 1e-4
opts.gradTol = 1e-7
opts.niter = 50

rso = inflation.ReducedSheetOptimizer(targetAttractedSheet, opts)
xorig = rso.getVars()
rso.useFirstOrderPrediction = True
rso.setVars(xorig + 1e-3 * np.random.uniform(low=-1, high=1, size=xorig.shape))
rso.useFirstOrderPrediction = False
rso.setVars(xorig + 1e-3 * np.random.uniform(low=-1, high=1, size=xorig.shape))
rso.useFirstOrderPrediction = True
rso.setVars(xorig + 1e-3 * np.random.uniform(low=-1, high=1, size=xorig.shape))
rso.useFirstOrderPrediction = False
rso.setVars(xorig + 1e-3 * np.random.uniform(low=-1, high=1, size=xorig.shape))
