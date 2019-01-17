#!/usr/bin/env python3
"""
 @file spm_driver.py

 SPM example to generate a sparse matrix from the spm drivers

 @copyright 2017-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
                      Univ. Bordeaux. All rights reserved.

 @version 6.0.0
 @author Pierre Ramet
 @author Mathieu Faverge
 @author Louis Poirel
 @date 2017-05-04

"""
import spm
import numpy as np

# Hack to make sure that the mkl is loaded
tmp = np.eye(2).dot(np.ones(2))

# Load a sparse matrix from the Laplacian driver
A = spm.spmatrix( None, driver=spm.driver.Laplacian, filename="10:10:10:2.:1." )

# Example from a HB file
#A = spm( None, driver=driver.HB, filename="$PASTIX_DIR/test/matrix/orsirr.rua" )

A.printInfo()

# Scale A for low-rank: A / ||A||_f
norm = A.norm()
A.scale( 1. / norm )

# Generate b and x0 vectors such that A * x0 = b
nrhs = 10
x0, b = A.genRHS( spm.rhstype.RndX, nrhs )

# Check that A * x = b
A.checkAxb( None, b, x0 )

