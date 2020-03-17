import scipy.sparse
import numpy as np
from time import time
# for direct
import scipy.sparse.linalg as sla
# for iterative
from scipy.optimize import minimize

from .progress_bar import printProgressBar


def MakeConsistent(matrix):
    from KratosMultiphysics.EigenSolversApplication import MakeConsistent as MakeConsistentCpp
    from KratosMultiphysics import Vector

    x = Vector()
    a = time()
    MakeConsistentCpp(matrix, x)
    print(time()-a)

    return matrix, x


def MakeConsistentScipy(matrix):
    # TODO precalculate col sum
    # TODO iterate non zeros only
    m = matrix.tocsc(copy=True)
    LHS = matrix.copy()
    for i in range(matrix.shape[0]):
        printProgressBar(i, matrix.shape[1])
        for j in range(matrix.shape[1]):
            m_ij = matrix[i,j]
            if m_ij == 0:
                continue
            v = (m_ij * m.getcol(j)).sum()
            if v != 0:
                LHS[i,j] = v

    x = sla.spsolve(LHS, [1]*matrix.shape[0])

    x = np.sqrt(x)

    _m = matrix.copy()
    return _m.dot(scipy.sparse.diags(x)), x


def MakeConsistentIterative(matrix):
    x = np.array([1.]*matrix.shape[0])

    def ModifyMatrix(_x):
        _m = matrix.copy()
        return _m @ scipy.sparse.diags(_x)

    def f(_x):
        a = ModifyMatrix(_x)
        m = a @ a.transpose()
        v = np.square((m.sum(1) - 1)).sum()
        return v

    result = minimize(f, x)
    return ModifyMatrix(x), x
