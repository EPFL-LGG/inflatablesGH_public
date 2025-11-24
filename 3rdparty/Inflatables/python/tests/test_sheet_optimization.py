import sys
sys.path.append('..')
import mesh, inflation, fd_validation, wall_generation, numpy as np

# m, fuseMarkers = wall_generation.triangulate_channel_walls(*parametric_pillows.concentricCircles(8, 50), 0.001)
m = mesh.Mesh('../../examples/single_tri.obj')
fuseMarkers = [1, 1, 1]

isheet = inflation.InflatableSheet(m, np.array(fuseMarkers) != 0)
isheet.setUseTensionFieldEnergy(False)
isheet.pressure = 30

targetSurf = mesh.Mesh('data/pringle.obj')
opt = inflation.SheetOptimizer(isheet, targetSurf)

xorig = opt.getVars()
xperturbed = xorig + 0.5 * np.random.uniform(low=-1, high=1, size=xorig.shape)
# Test the collapse barrier by scaling everything down past its activation point.
opt.setVars(xperturbed / np.sqrt(20 + 5))

# fd_perturb = np.zeros_like(xperturbed)
# fd_perturb[opt.numEquilibriumVars()] = 1.0

# print(fd_validation.validateHessian(opt, fd_eps=1e-6, perturb=fd_perturb, etype=opt.EnergyType.CollapseBarrier))
print(fd_validation.validateGrad   (opt, fd_eps=1e-7))
print(fd_validation.validateHessian(opt, fd_eps=1e-6, etype=opt.EnergyType.CollapseBarrier))
