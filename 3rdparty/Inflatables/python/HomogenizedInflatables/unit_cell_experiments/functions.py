"""
Sample code automatically generated on 2023-02-19 16:53:01

by www.matrixcalculus.org

from input

d/dw w*w'*v - (w'*w)*v = w'*v*eye+w*v'-2*v*w'

where

v is a vector
w is a vector

The generated code is provided "as is" without warranty of any kind.
"""

from __future__ import division, print_function, absolute_import

import numpy as np

def fAndG(v, w):
    assert isinstance(v, np.ndarray)
    dim = v.shape
    assert len(dim) == 1
    v_rows = dim[0]
    assert isinstance(w, np.ndarray)
    dim = w.shape
    assert len(dim) == 1
    w_rows = dim[0]
    assert w_rows == v_rows

    t_0 = (w).dot(v)
    functionValue = ((t_0 * w) - ((w).dot(w) * v))
    gradient = (((t_0 * np.eye(w_rows, w_rows)) + np.outer(w, v)) - (2 * np.outer(v, w)))

    return functionValue, gradient

def checkGradient(v, w):
    # numerical gradient checking
    # f(x + t * delta) - f(x - t * delta) / (2t)
    # should be roughly equal to inner product <g, delta>
    t = 1E-6
    delta = np.random.randn(3)
    f1, _ = fAndG(v, w + t * delta)
    f2, _ = fAndG(v, w - t * delta)
    f, g = fAndG(v, w)
    print('approximation error',
          np.linalg.norm((f1 - f2) / (2*t) - np.tensordot(g, delta, axes=1)))

def generateRandomData():
    v = np.random.randn(3)
    w = np.random.randn(3)

    return v, w

if __name__ == '__main__':
    v, w = generateRandomData()
    functionValue, gradient = fAndG(v, w)
    print('functionValue = ', functionValue)
    print('gradient = ', gradient)

    print('numerical gradient checking ...')
    checkGradient(v, w)
