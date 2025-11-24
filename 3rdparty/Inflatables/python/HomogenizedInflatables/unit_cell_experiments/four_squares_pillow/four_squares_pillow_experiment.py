import sys
sys.path.append("../")
sys.path.append("../../")
sys.path.append("../../../")
sys.path.append("../../../gmsh")

import experiment_helper
import igl
from periodic_simulation_setup import *
import json

res = 0.1

# type_index = 4
for type_index in range(5):

    isheet, m, marker = pattern_generator_using_gmsh.get_four_square(5, type_index, avg_len_boundary = res, avg_len_embeddings = res)


    visualization.plot_2d_mesh(m, pointList = marker, width = 5, height = 5)

    viewer = TriMeshViewer(isheet, width=768, height=640)


    viewer.showWireframe(False)

    viewer.show()

    # pressure = 0.8
    stiffness_pressure = 0.1


    allowBending = False
    useTFT = True
    disableFusedRegionTFT = False

    from matplotlib import cm

    def cb(i):
        viewer.update(scalarField=utils.getStrains(isheet)[:, 0])


    # Choose strategy for constraining rigid motion
    fixedVars, hessianShift = [], 1e-6

    isheet.setUseTensionFieldEnergy(useTFT)
    isheet.setUseHessianProjectedEnergy(False)
    if (disableFusedRegionTFT):
        isheet.disableFusedRegionTensionFieldTheory(False)
    isheet.pressure = stiffness_pressure

    benchmark.reset()
    print(allowBending, stiffness_pressure, hessianShift, fixedVars)

    opts.niter = 500
    opts.gradTol = 1e-10
    print(opts.factorizer)

    cr = inflation.inflation_newton(isheet, fixedVars, opts, callback=cb, hessianShift = hessianShift)
    viewer.update(scalarField=utils.getStrains(isheet)[:, 0])
    benchmark.report()

    print("max strain: ", max(utils.getStrains(isheet)[:, 0]))


    strains = utils.getStrains(isheet)[:, 0]
    strainField = vis.fields.ScalarField(isheet, strains, colormap = cm.viridis, vmin= 0, vmax = max(utils.getStrains(isheet)[:, 0]))
    viewer.update(scalarField=strainField)
    viewer.saveColorizedObj("four_square_pillow_{}.obj".format(type_index))
