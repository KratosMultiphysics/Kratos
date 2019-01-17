"""

 @file __spm__.py

 SPM python wrapper

 @copyright 2017-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
                      Univ. Bordeaux. All rights reserved.

 @version 6.0.0
 @author Pierre Ramet
 @author Mathieu Faverge
 @author Louis Poirel
 @date 2017-05-04

This file has been automatically generated with gen_wrappers.py

"""
from ctypes import *
import numpy as np

from . import libspm
from .enum import __spm_int__

class pyspm_spmatrix_t(Structure):
    _fields_ = [("mtxtype",   c_int               ),
                ("flttype",   c_int               ),
                ("fmttype",   c_int               ),
                ("gN",        __spm_int__         ),
                ("n",         __spm_int__         ),
                ("gnnz",      __spm_int__         ),
                ("nnz",       __spm_int__         ),
                ("gNexp",     __spm_int__         ),
                ("nexp",      __spm_int__         ),
                ("gnnzexp",   __spm_int__         ),
                ("nnzexp",    __spm_int__         ),
                ("dof",       __spm_int__         ),
                ("dofs",      POINTER(__spm_int__)),
                ("layout",    c_int               ),
                ("colptr",    POINTER(__spm_int__)),
                ("rowptr",    POINTER(__spm_int__)),
                ("loc2glob",  POINTER(__spm_int__)),
                ("values",    c_void_p            ) ]

def pyspm_spmInit( spm ):
    libspm.spmInit.argtypes = [ POINTER(pyspm_spmatrix_t) ]
    libspm.spmInit( spm )

def pyspm_spmAlloc( spm ):
    libspm.spmAlloc.argtypes = [ POINTER(pyspm_spmatrix_t) ]
    libspm.spmAlloc( spm )

def pyspm_spmExit( spm ):
    libspm.spmExit.argtypes = [ POINTER(pyspm_spmatrix_t) ]
    libspm.spmExit( spm )

def pyspm_spmCopy( spm ):
    libspm.spmCopy.argtypes = [ POINTER(pyspm_spmatrix_t) ]
    libspm.spmCopy.restype = POINTER(pyspm_spmatrix_t)
    return libspm.spmCopy( spm )

def pyspm_spmBase( spm, baseval ):
    libspm.spmBase.argtypes = [ POINTER(pyspm_spmatrix_t), c_int ]
    libspm.spmBase( spm, baseval )

def pyspm_spmFindBase( spm ):
    libspm.spmFindBase.argtypes = [ POINTER(pyspm_spmatrix_t) ]
    libspm.spmFindBase.restype = __spm_int__
    return libspm.spmFindBase( spm )

def pyspm_spmConvert( ofmttype, ospm ):
    libspm.spmConvert.argtypes = [ c_int, POINTER(pyspm_spmatrix_t) ]
    libspm.spmConvert.restype = c_int
    return libspm.spmConvert( ofmttype, ospm )

def pyspm_spmUpdateComputedFields( spm ):
    libspm.spmUpdateComputedFields.argtypes = [ POINTER(pyspm_spmatrix_t) ]
    libspm.spmUpdateComputedFields( spm )

def pyspm_spmGenFakeValues( spm ):
    libspm.spmGenFakeValues.argtypes = [ POINTER(pyspm_spmatrix_t) ]
    libspm.spmGenFakeValues( spm )

def pyspm_spmNorm( ntype, spm ):
    libspm.spmNorm.argtypes = [ c_int, POINTER(pyspm_spmatrix_t) ]
    libspm.spmNorm.restype = c_double
    return libspm.spmNorm( ntype, spm )

def pyspm_spmMatVec( trans, alpha, spm, x, beta, y ):
    libspm.spmMatVec.argtypes = [ c_int, c_double, POINTER(pyspm_spmatrix_t),
                                  c_void_p, c_double, c_void_p ]
    libspm.spmMatVec.restype = c_int
    return libspm.spmMatVec( trans, alpha, spm, x, beta, y )

def pyspm_spmMatMat( trans, n, alpha, A, B, ldb, beta, C, ldc ):
    libspm.spmMatMat.argtypes = [ c_int, __spm_int__, c_double,
                                  POINTER(pyspm_spmatrix_t), c_void_p,
                                  __spm_int__, c_double, c_void_p, __spm_int__ ]
    libspm.spmMatMat.restype = c_int
    return libspm.spmMatMat( trans, n, alpha, A, B, ldb, beta, C, ldc )

def pyspm_spmScalMatrix( alpha, spm ):
    libspm.spmScalMatrix.argtypes = [ c_double, POINTER(pyspm_spmatrix_t) ]
    libspm.spmScalMatrix( alpha, spm )

def pyspm_spmScalVector( flt, alpha, n, x, incx ):
    libspm.spmScalVector.argtypes = [ c_int, c_double, __spm_int__, c_void_p,
                                      __spm_int__ ]
    libspm.spmScalVector( flt, alpha, n, x, incx )

def pyspm_spmSort( spm ):
    libspm.spmSort.argtypes = [ POINTER(pyspm_spmatrix_t) ]
    libspm.spmSort.restype = c_int
    return libspm.spmSort( spm )

def pyspm_spmMergeDuplicate( spm ):
    libspm.spmMergeDuplicate.argtypes = [ POINTER(pyspm_spmatrix_t) ]
    libspm.spmMergeDuplicate.restype = __spm_int__
    return libspm.spmMergeDuplicate( spm )

def pyspm_spmSymmetrize( spm ):
    libspm.spmSymmetrize.argtypes = [ POINTER(pyspm_spmatrix_t) ]
    libspm.spmSymmetrize.restype = __spm_int__
    return libspm.spmSymmetrize( spm )

def pyspm_spmCheckAndCorrect( spm_in, spm_out ):
    libspm.spmCheckAndCorrect.argtypes = [ POINTER(pyspm_spmatrix_t),
                                           POINTER(pyspm_spmatrix_t) ]
    libspm.spmCheckAndCorrect.restype = c_int
    return libspm.spmCheckAndCorrect( spm_in, spm_out )

def pyspm_spmGenRHS( type, nrhs, spm, x, ldx, b, ldb ):
    libspm.spmGenRHS.argtypes = [ c_int, __spm_int__, POINTER(pyspm_spmatrix_t),
                                  c_void_p, __spm_int__, c_void_p, __spm_int__ ]
    libspm.spmGenRHS.restype = c_int
    return libspm.spmGenRHS( type, nrhs, spm, x, ldx, b, ldb )

def pyspm_spmCheckAxb( eps, nrhs, spm, x0, ldx0, b, ldb, x, ldx ):
    libspm.spmCheckAxb.argtypes = [ c_double, __spm_int__,
                                    POINTER(pyspm_spmatrix_t), c_void_p,
                                    __spm_int__, c_void_p, __spm_int__,
                                    c_void_p, __spm_int__ ]
    libspm.spmCheckAxb.restype = c_int
    return libspm.spmCheckAxb( eps, nrhs, spm, x0, ldx0, b, ldb, x, ldx )

def pyspm_spmIntConvert( n, input ):
    libspm.spmIntConvert.argtypes = [ __spm_int__, c_int_p ]
    libspm.spmIntConvert.restype = POINTER(__spm_int__)
    return libspm.spmIntConvert( n, input )

def pyspm_spmLoad( spm ):
    libspm.spmLoad.argtypes = [ POINTER(pyspm_spmatrix_t), c_void_p ]
    libspm.spmLoad.restype = c_int
    return libspm.spmLoad( spm, None )

def pyspm_spmSave( spm ):
    libspm.spmSave.argtypes = [ POINTER(pyspm_spmatrix_t), c_void_p ]
    libspm.spmSave.restype = c_int
    return libspm.spmSave( spm, None )

def pyspm_spmReadDriver( driver, filename, spm ):
    libspm.spmReadDriver.argtypes = [ c_int, c_char_p,
                                      POINTER(pyspm_spmatrix_t) ]
    libspm.spmReadDriver.restype = c_int
    return libspm.spmReadDriver( driver, filename, spm )

def pyspm_spmParseLaplacianInfo( filename, flttype, dim1, dim2, dim3, alpha,
                                 beta ):
    libspm.spmParseLaplacianInfo.argtypes = [ c_char_p, c_int_p,
                                              POINTER(__spm_int__),
                                              POINTER(__spm_int__),
                                              POINTER(__spm_int__),
                                              POINTER(c_double),
                                              POINTER(c_double) ]
    libspm.spmParseLaplacianInfo.restype = c_int
    return libspm.spmParseLaplacianInfo( filename, flttype, dim1, dim2, dim3,
                                         alpha, beta )

def pyspm_spm2Dense( spm ):
    libspm.spm2Dense.argtypes = [ POINTER(pyspm_spmatrix_t) ]
    libspm.spm2Dense( spm )

def pyspm_spmPrint( spm ):
    libspm.spmPrint.argtypes = [ POINTER(pyspm_spmatrix_t), c_void_p ]
    libspm.spmPrint( spm, None )

def pyspm_spmPrintInfo( spm ):
    libspm.spmPrintInfo.argtypes = [ POINTER(pyspm_spmatrix_t), c_void_p ]
    libspm.spmPrintInfo( spm, None )

def pyspm_spmExpand( spm ):
    libspm.spmExpand.argtypes = [ POINTER(pyspm_spmatrix_t) ]
    libspm.spmExpand.restype = POINTER(pyspm_spmatrix_t)
    return libspm.spmExpand( spm )

def pyspm_spmDofExtend( spm, type, dof ):
    libspm.spmDofExtend.argtypes = [ POINTER(pyspm_spmatrix_t), c_int, c_int ]
    libspm.spmDofExtend.restype = POINTER(pyspm_spmatrix_t)
    return libspm.spmDofExtend( spm, type, dof )


