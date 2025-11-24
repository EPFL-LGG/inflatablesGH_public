import sys; sys.path.append('../..')
import MeshFEM, mesh, sparse_matrices, benchmark, field_sampler, mesh_utilities
import inflatables_parametrization as parametrization, numpy as np, importlib, pickle, wall_generation
import utils
import py_newton_optimizer
from py_newton_optimizer import NewtonOptimizerOptions
from numpy.linalg import norm
from io_redirection import suppress_stdout
import visualization, wall_width_formulas as wwf
from tri_mesh_viewer import TriMeshViewer

import parallelism
parallelism.set_max_num_tbb_threads(4)

variants = ['large', 'medium', 'small']

if len(sys.argv) != 2 or sys.argv[1] not in variants:
    print("usage: python run_opt.py variant")
    print("where variant is one of " + ",".join(variants))
    sys.exit(-1)

variant = sys.argv[1]

target_surf = mesh.Mesh('../../SiggraphExamples/Meshes/20200118_bike_helmet_v1_R00.obj')
target_surf.setVertices(utils.prototypeScaleNormalization(target_surf.vertices(), placeAtopFloor=True))
target_surf = mesh_utilities.subdivide_loop(target_surf, 1)

targetAttractedSheet = utils.load(f'data/helmet_{variant}_tas_init.pkl.gz')
isheet = targetAttractedSheet.sheet()

import boundaries
fixedVars = boundaries.getOuterBoundaryVars(isheet)

import sheet_optimizer, opt_config
origDesignMesh = isheet.mesh().copy()

def config(so):
    fcs = so.rso.fusingCurveSmoothness()
    fcs.interiorWeight = 1/10

sheet_opt = sheet_optimizer.PySheetOptimizer(targetAttractedSheet, fixedVars, renderMode=sheet_optimizer.RenderMode.OFFSCREEN, screenshotPath=f'{variant}_tubes.mp4',
                                             detActivationThreshold=0.9, detActivationThresholdTubeTri=0.5,
                                             originalDesignMesh=origDesignMesh, fusingCurveSmoothnessConfig=opt_config.FusingCurveSmoothnessParams(0.0, 0.0, 1.0, 1.0), customConfigCallback=config)

sheet_opt.deploy_viewer.setCameraParams(((2.2283129063050544, -2.4424908418018645, 1.480176156096905),
     (-0.20481495134223526, 0.5059968464853937, 0.8378651604247005),
      (-0.10862883029189645, -0.1146986958685825, -0.4968676111207299)))
sheet_opt.flat_viewer.setCameraParams(((-0.005724654779466134, -0.15742800643532207, 3.663779058858404),
     (0.0, 1.0, 0.0),
      (-0.005724654779466134, -0.15742800643532207, 0.0)))
sheet_opt.flat_viewer.showWireframe()

sheet_opt.deploy_viewer.scalarFieldGetter = visualization.ISheetScalarField.TGT_DIST(targetAttractedSheet.sheet(), target_surf)

sheet_opt.optimize()
sheet_opt.save(f'data/helmet_{variant}.opt.pkl.gz')
