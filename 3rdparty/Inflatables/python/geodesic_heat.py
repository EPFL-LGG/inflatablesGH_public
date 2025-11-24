import numpy as np, scipy, scipy.sparse, scipy.sparse.linalg
import MeshFEM
import differential_operators

def computeX(m, sourceVertices, c = 4 / np.sqrt(3)):
    # Choose a timestep proportional to h^2 where h is the average edge length.
    # (As discussed in section 3.2.4 of the paper)
    t = c *  m.volume / m.numElements()

    L = differential_operators.laplacian(m).compressedColumn()
    M = differential_operators.mass(m, lumped=False).compressedColumn()
    A = L + t * M

    mask = np.ones(m.numVertices(), dtype=np.bool)
    mask[sourceVertices] = False

    A_ff = A[:,  mask][mask, :]
    A_fc = A[:, ~mask][mask, :]

    # Solve (M + t L) u = 0 with the constraint u[sourceVertices] = 1
    u = np.ones(m.numVertices())
    u[mask] = scipy.sparse.linalg.spsolve(A_ff, -A_fc @ np.ones(len(sourceVertices)))

    # Compute the heat gradients
    g = differential_operators.gradient(m, u)
    # Normalize the gradients to get an approximate gradient of the distance field
    X = -g / np.linalg.norm(g, axis=1)[:, np.newaxis]
    return X

def poissonSolve(m, sourceVertices, X):
    mask = np.ones(m.numVertices(), dtype=np.bool)
    mask[sourceVertices] = False
    
    L = differential_operators.laplacian(m).compressedColumn()
    divX = differential_operators.divergence(m, X)
    L_ff = L[:, mask][mask, :]
    heatDist = np.zeros(m.numVertices())
    heatDist[mask] = scipy.sparse.linalg.spsolve(L_ff, divX[mask]) 
    return heatDist

def dist(m, sourceVertices):
    return poissonSolve(m, sourceVertices, computeX(m, sourceVertices))
