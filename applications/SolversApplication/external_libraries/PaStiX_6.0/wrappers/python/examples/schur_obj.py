#!/usr/bin/env python3
"""
 @file schur_obj.py

 PaStiX Schur python example with an object oriented programing solution.

 @copyright 2017-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
                      Univ. Bordeaux. All rights reserved.

 @version 6.0.1
 @author Pierre Ramet
 @author Mathieu Faverge
 @author Louis Poirel
 @date 2018-07-16

 This example shows how to use pastix solver to solve a system with a
 Schur complement.

"""
import pypastix as pastix
import scipy.sparse as sps
import scipy.linalg as la
import numpy as np

# Hack to make sure that the mkl is loaded
tmp = np.eye(2).dot(np.ones(2))

# Set matrix A from Scipy Sparse matrices
n = 30
A = sps.spdiags([np.ones(n)*i for i in [4., -1, -1]],
                [0, 1, -1], n, n)

x0 = np.arange(n).reshape(n,1)
# Construct b as b = A * x_0
b = A.dot(x0)
x = b.copy()

# Initialize the Schur list as the first third of the elements
nschur = min( int(n / 3), 5000 )
schurlist = np.arange(nschur)

solver = pastix.solver(verbose=2)

solver.schur(A, schurlist)
S = solver.S
f = solver.schur_forward(b)
y = la.solve(S, f, sym_pos=True, lower=True)
x = solver.schur_backward(y, b)

# Check the solution
rc = solver.check( x, b, x0=x0 )

solver.finalize()

exit(rc)
