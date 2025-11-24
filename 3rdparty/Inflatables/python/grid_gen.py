import numpy as np
def triangulatedGrid2D(nx, ny, triangulationType=0):
    if triangulationType not in [0, 1, 2]: raise Exception('Invalid triangulation type')
    X, Y = np.meshgrid(np.linspace(0, 0.5, nx), np.linspace(0, 1, ny))
    V = np.column_stack((X.ravel(), Y.ravel()))
    flatIdx = lambda i, j: j + i * nx
    def triCorners(i, j):
        corners = [[[[i, j], [i, j + 1], [i + 1, j]],         # lower left tri
                   [[i, j + 1], [i + 1, j + 1], [i + 1, j]]], # upper right tri
                   [[[i, j], [i, j + 1], [i + 1, j + 1]],     # upper left tri
                   [[i, j], [i + 1, j + 1], [i + 1, j]]]]     # bottom right tri

        tt = triangulationType if triangulationType in [0, 1] else (i + j) % 2
        return [[flatIdx(c[0], c[1]) for c in cc] for cc in corners[tt]]
    F = [tri for i in range(ny - 1) for j in range(nx - 1) for tri in triCorners(i, j)]
    return V, F
