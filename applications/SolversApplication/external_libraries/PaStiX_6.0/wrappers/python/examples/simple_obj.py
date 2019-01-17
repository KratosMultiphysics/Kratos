#!/usr/bin/env python3
"""
 @file simple_obj.py

 PaStiX simple python example with the solver object

 @copyright 2017-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
                      Univ. Bordeaux. All rights reserved.

 @version 6.0.1
 @author Pierre Ramet
 @author Mathieu Faverge
 @author Louis Poirel
 @date 2018-07-16

"""
import pypastix as pastix
import scipy.sparse as sps
import numpy as np

# Set the problem
n = 9
A = sps.spdiags([np.ones(n)*i for i in [4, -1, -1, -1, -1]],
                [0, 1, 3, -1, -3], n, n)
x0 = np.arange(n)
b = A.dot(x0)

# Hack to make sure that the mkl is loaded
tmp = np.eye(2).dot(np.ones(2))

# Factorize
solver = pastix.solver(A, verbose=True)

# Solve
x = solver.solve( b )

# Check the solution
rc = solver.check( x, b, x0=x0 )

solver.finalize()

exit(rc)
