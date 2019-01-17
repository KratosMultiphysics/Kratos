#!/usr/bin/env python3
"""
 @file schur.py

 PaStiX schur python example

 @copyright 2017-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
                      Univ. Bordeaux. All rights reserved.

 @version 6.0.1
 @author Pierre Ramet
 @author Mathieu Faverge
 @author Louis Poirel
 @date 2018-07-16

"""
import pypastix as pastix
import spm
import scipy.sparse as sps
import scipy.linalg as la
import numpy as np

# Set matrix A
if True:
    # From Scipy Sparse matrices
    n = 30
    A = sps.spdiags([np.ones(n)*i for i in [4., -1, -1]],
                    [0, 1, -1], n, n)

    x0 = np.arange(n).reshape(n,1)
    # Construct b as b = A * x_0
    b = A.dot(x0)
    spmA = spm.spmatrix( A )

else:
    # Load a sparse matrix from RSA driver
    spmA = spm.spmatrix( None, driver=spm.driver.Laplacian, filename="10:10:10:2.:1." )
    n = spmA.spm_c.gN

    # Generate b and x0 vector such that A * x0 = b
    x0, b = spmA.genRHS( spm.rhstype.RndX )

x = b.copy()

# Hack to make sure that the mkl is loaded
tmp = np.eye(2).dot(np.ones(2))

# Convert the scipy sparse matrix to spm storage format
spmA.printInfo()

# Initialize parameters to default values
iparm, dparm = pastix.initParam()

# Startup PaStiX
pastix_data = pastix.init( iparm, dparm )

factotype = pastix.factotype.LLT
iparm[pastix.iparm.factorization]   = factotype
iparm[pastix.iparm.scheduler]       = pastix.scheduler.Sequential
iparm[pastix.iparm.schur_solv_mode] = pastix.solv_mode.Interface

# Initialize the Schur list as the first third of the elements
nschur = min( int(n / 3.), 5000 )
schurlist = np.arange(nschur) + spmA.findBase()
pastix.setSchurUnknownList ( pastix_data, schurlist )

# Perform analyze
pastix.task_analyze( pastix_data, spmA )

# Perform numerical factorization
pastix.task_numfact( pastix_data, spmA )

# Get the Schur complement
S = np.array( np.zeros( (nschur, nschur) ), order='F', dtype=spmA.dtype )
pastix.getSchur( pastix_data, S )

# Store both sides for linalg
if factotype != pastix.factotype.LU:
    S += la.tril(S, -1).T

# 1- Apply P to b
pastix.subtask_applyorder( pastix_data, pastix.dir.Forward, x )

if factotype == pastix.factotype.LU:
    # 2- Forward solve on the non Schur complement part of the system
    pastix.subtask_trsm( pastix_data, pastix.side.Left, pastix.uplo.Lower, pastix.trans.NoTrans, pastix.diag.Unit, x )

    # 3- Solve the Schur complement part
    x[-nschur:] = la.solve(S, x[-nschur:], sym_pos=False)

    # 4- Backward solve on the non Schur complement part of the system
    pastix.subtask_trsm( pastix_data, pastix.side.Left, pastix.uplo.Upper, pastix.trans.NoTrans, pastix.diag.NonUnit, x )

else:
    # 2- Forward solve on the non Schur complement part of the system
    pastix.subtask_trsm( pastix_data, pastix.side.Left, pastix.uplo.Lower, pastix.trans.NoTrans, pastix.diag.NonUnit, x )

    # 3- Solve the Schur complement part
    x[-nschur:] = la.solve(S, x[-nschur:], sym_pos=True, lower=True)

    # 4- Backward solve on the non Schur complement part of the system
    pastix.subtask_trsm( pastix_data, pastix.side.Left, pastix.uplo.Lower, pastix.trans.ConjTrans, pastix.diag.NonUnit, x )

# 5- Apply P^t to x
pastix.subtask_applyorder( pastix_data, pastix.dir.Backward, x )

# Check solution
rc = spmA.checkAxb( x0, b, x )

pastix.finalize( pastix_data, iparm, dparm )

exit(rc)
