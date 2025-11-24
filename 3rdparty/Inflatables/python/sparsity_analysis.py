import scipy, sparse_matrices, numpy as np, matplotlib.pyplot as plt
from scipy.sparse import csc_matrix

def to_csc(K): return csc_matrix((K.Ax, K.Ai, K.Ap))

def spy(K, width=6, height=6, markersize=3):
    fig = plt.figure(figsize=(width, height))
    plt.spy(to_csc(K), markersize=markersize)
    plt.show()

# Verify that Ksp holds the sparsity pattern of K.
# Returns "True" if the sparsity patterns match. Otherwise, it returns
# a sparse matrix holding the difference of the two sparsity patterns.
def validate_sparsity_pattern(K, Ksp):
    K_csc   = to_csc(K)
    Ksp_csc = to_csc(Ksp)

    Ksp_csc.sum_duplicates()
    if ((np.min(Ksp_csc.data) != 1.0) or (np.max(Ksp_csc.data) != 1.0)): raise Exception('Invalid sparsity pattern (duplicate/nonunit entries)')

    K_csc.data.fill(1.0)
    diff = Ksp_csc - K_csc
    if (diff.nnz == 0): return True
    return diff
