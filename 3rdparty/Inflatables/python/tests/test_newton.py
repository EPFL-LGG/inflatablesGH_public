import parametrization, sparse_matrices, numpy as np
from numpy.linalg import norm

m = parametrization.Mesh("../examples/julius.msh")

lg = parametrization.LocalGlobalParametrizer(m, parametrization.lscm(m))
lg.alphaMin = 1
lg.alphaMax = np.pi / 2
lg.energy()

for i in range(4):
    lg.runIteration()
    print(lg.energy())

rg = parametrization.RegularizedParametrizer(lg)
print(rg.energy())

# Fix a single (u, v) and rotation of the parametric domain to make parametrization unique
fixedVars = [rg.uOffset(), rg.vOffset(), rg.phiOffset()]

#my newton solver makes very slow progress due to many variables entering working set.
#also, the energy appears to be nonconvex even quite close to the optimum...
opts = parametrization.NewtonOptimizerOptions()
opts.useIdentityMetric = True
parametrization.regularized_parametrization_newton(rg, fixedVars)
