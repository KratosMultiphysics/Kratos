# ==============================================================================
#  KratosShapeOptimizationApplication
#
#  License:         BSD License
#                   license: ShapeOptimizationApplication/license.txt
#
#  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
#
# ==============================================================================

# Making KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

# Import additional libraries
import math

# ==============================================================================
def Norm2(_X):
    return math.sqrt(sum(x**2 for x in _X))

# ------------------------------------------------------------------------------
def NormInf3D(_X):
    temp_vec = [_X[3*i]**2 + _X[3*i+1]**2 + _X[3*i+2]**2 for i in range(int(len(_X)/3)) ]
    max_squared_value = max(temp_vec)
    return math.sqrt(max_squared_value)

# ------------------------------------------------------------------------------
def DotProduct(_X,_Y):
    if len(_X) != len(_Y):
        raise RuntimeError("Dot product to be computed but _X and _Y do not have the same dimension!")
    return sum([_X[i]*_Y[i] for i in range(len(_X))])

# ------------------------------------------------------------------------------
def ScalarVectorProduct(scal, _X):
    return [x*scal for x in _X]

# ------------------------------------------------------------------------------
def ProjectOntoInterval(value, bound):
    if bound<0:
        raise ValueError("bound negative")
    return min(max(value,-bound),bound)

# ------------------------------------------------------------------------------
def Zeros(n):
    return [0.0 for i in range(n)]

# ==============================================================================
