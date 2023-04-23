import KratosMultiphysics.OptimizationApplication as KratosOA
import numpy as np
import pybind11
import nlopt

def myfunc(x,grad,z):
    if len(grad) > 0:
       grad[0] = 0.0
       grad[1] = 0.5 / np.sqrt(x[1])
    print("successful")
    return np.sqrt(x[1])


opt = KratosOA.NLOptOptimizer()
opt.set_function(myfunc)
opt.optimize()

