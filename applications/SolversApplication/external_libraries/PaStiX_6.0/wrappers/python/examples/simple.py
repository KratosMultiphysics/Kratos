#!/usr/bin/env python3
"""
 @file simple.py

 PaStiX simple python example

 @copyright 2017-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
                      Univ. Bordeaux. All rights reserved.

 @version 6.0.1
 @author Pierre Ramet
 @author Mathieu Faverge
 @author Louis Poirel
 @date 2018-07-16

"""
import spm
import pypastix as pastix
import scipy.sparse as sps
import numpy as np

# Hack to make sure that the mkl is loaded
tmp = np.eye(2).dot(np.ones(2))

# Load a sparse matrix from HB driver
spmA = spm.spmatrix( None, driver=spm.driver.Laplacian, filename="10:10:10:2.:1." )
#spmA = spm.spmatrix( None, driver=spm.driver.HB, filename="$PASTIX_DIR/test/matrix/orsirr.rua" )
spmA.printInfo()

# Scale A for low-rank: A / ||A||_f
norm = spmA.norm()
spmA.scale( 1. / norm )

# Generate b and x0 vector such that A * x0 = b
nrhs = 1
x0, b = spmA.genRHS( spm.rhstype.RndX, nrhs )

# Initialize parameters to default values
iparm, dparm = pastix.initParam()

# Startup PaStiX
pastix_data = pastix.init( iparm, dparm )

# Change some parameters
iparm[pastix.iparm.factorization] = pastix.factotype.LU

# Perform analyze
pastix.task_analyze( pastix_data, spmA )

# Perform numerical factorization
pastix.task_numfact( pastix_data, spmA )

# Perform solve
x = b.copy()
pastix.task_solve( pastix_data, spmA, x )

# Refine the solution
pastix.task_refine( pastix_data, spmA, b, x )

# Check solution
spmA.checkAxb( x0, b, x )

pastix.finalize( pastix_data, iparm, dparm )
