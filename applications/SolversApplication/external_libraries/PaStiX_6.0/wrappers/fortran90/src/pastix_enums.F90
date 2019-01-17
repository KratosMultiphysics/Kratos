
!
! @file pastix_enums.F90
!
! PaStiX fortran 90 wrapper to define enums and datatypes
!
! @copyright 2017-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
!                      Univ. Bordeaux. All rights reserved.
!
! @version 6.0.1
! @author Mathieu Faverge
! @date 2018-07-16
!
! This file has been automatically generated with gen_wrappers.py
!
module pastix_enums
  use iso_c_binding
  implicit none

  ! enum iparm
  enum, bind(C)
     enumerator :: IPARM_VERBOSE               = 1
     enumerator :: IPARM_IO_STRATEGY           = 2
     enumerator :: IPARM_NNZEROS               = 3
     enumerator :: IPARM_NNZEROS_BLOCK_LOCAL   = 4
     enumerator :: IPARM_ALLOCATED_TERMS       = 5
     enumerator :: IPARM_PRODUCE_STATS         = 6
     enumerator :: IPARM_MC64                  = 7
     enumerator :: IPARM_ORDERING              = 8
     enumerator :: IPARM_ORDERING_DEFAULT      = 9
     enumerator :: IPARM_SCOTCH_SWITCH_LEVEL   = 10
     enumerator :: IPARM_SCOTCH_CMIN           = 11
     enumerator :: IPARM_SCOTCH_CMAX           = 12
     enumerator :: IPARM_SCOTCH_FRAT           = 13
     enumerator :: IPARM_METIS_CTYPE           = 14
     enumerator :: IPARM_METIS_RTYPE           = 15
     enumerator :: IPARM_METIS_NO2HOP          = 16
     enumerator :: IPARM_METIS_NSEPS           = 17
     enumerator :: IPARM_METIS_NITER           = 18
     enumerator :: IPARM_METIS_UFACTOR         = 19
     enumerator :: IPARM_METIS_COMPRESS        = 20
     enumerator :: IPARM_METIS_CCORDER         = 21
     enumerator :: IPARM_METIS_PFACTOR         = 22
     enumerator :: IPARM_METIS_SEED            = 23
     enumerator :: IPARM_METIS_DBGLVL          = 24
     enumerator :: IPARM_SF_KASS               = 25
     enumerator :: IPARM_AMALGAMATION_LVLBLAS  = 26
     enumerator :: IPARM_AMALGAMATION_LVLCBLK  = 27
     enumerator :: IPARM_REORDERING_SPLIT      = 28
     enumerator :: IPARM_REORDERING_STOP       = 29
     enumerator :: IPARM_MIN_BLOCKSIZE         = 30
     enumerator :: IPARM_MAX_BLOCKSIZE         = 31
     enumerator :: IPARM_TASKS2D_LEVEL         = 32
     enumerator :: IPARM_TASKS2D_WIDTH         = 33
     enumerator :: IPARM_ABS                   = 34
     enumerator :: IPARM_ALLCAND               = 35
     enumerator :: IPARM_INCOMPLETE            = 36
     enumerator :: IPARM_LEVEL_OF_FILL         = 37
     enumerator :: IPARM_FACTORIZATION         = 38
     enumerator :: IPARM_STATIC_PIVOTING       = 39
     enumerator :: IPARM_FREE_CSCUSER          = 40
     enumerator :: IPARM_SCHUR_FACT_MODE       = 41
     enumerator :: IPARM_SCHUR_SOLV_MODE       = 42
     enumerator :: IPARM_APPLYPERM_WS          = 43
     enumerator :: IPARM_REFINEMENT            = 44
     enumerator :: IPARM_NBITER                = 45
     enumerator :: IPARM_ITERMAX               = 46
     enumerator :: IPARM_GMRES_IM              = 47
     enumerator :: IPARM_SCHEDULER             = 48
     enumerator :: IPARM_THREAD_NBR            = 49
     enumerator :: IPARM_AUTOSPLIT_COMM        = 50
     enumerator :: IPARM_GPU_NBR               = 51
     enumerator :: IPARM_GPU_MEMORY_PERCENTAGE = 52
     enumerator :: IPARM_GPU_MEMORY_BLOCK_SIZE = 53
     enumerator :: IPARM_COMPRESS_MIN_WIDTH    = 54
     enumerator :: IPARM_COMPRESS_MIN_HEIGHT   = 55
     enumerator :: IPARM_COMPRESS_WHEN         = 56
     enumerator :: IPARM_COMPRESS_METHOD       = 57
     enumerator :: IPARM_COMPRESS_ORTHO        = 58
     enumerator :: IPARM_THREAD_COMM_MODE      = 59
     enumerator :: IPARM_MODIFY_PARAMETER      = 60
     enumerator :: IPARM_START_TASK            = 61
     enumerator :: IPARM_END_TASK              = 62
     enumerator :: IPARM_FLOAT                 = 63
     enumerator :: IPARM_MTX_TYPE              = 64
     enumerator :: IPARM_DOF_NBR               = 65
     enumerator :: IPARM_SIZE                  = 65
  end enum

  ! enum dparm
  enum, bind(C)
     enumerator :: DPARM_FILL_IN            = 1
     enumerator :: DPARM_EPSILON_REFINEMENT = 2
     enumerator :: DPARM_RELATIVE_ERROR     = 3
     enumerator :: DPARM_EPSILON_MAGN_CTRL  = 4
     enumerator :: DPARM_ANALYZE_TIME       = 5
     enumerator :: DPARM_PRED_FACT_TIME     = 6
     enumerator :: DPARM_FACT_TIME          = 7
     enumerator :: DPARM_SOLV_TIME          = 8
     enumerator :: DPARM_FACT_FLOPS         = 9
     enumerator :: DPARM_FACT_THFLOPS       = 10
     enumerator :: DPARM_FACT_RLFLOPS       = 11
     enumerator :: DPARM_SOLV_FLOPS         = 12
     enumerator :: DPARM_SOLV_THFLOPS       = 13
     enumerator :: DPARM_SOLV_RLFLOPS       = 14
     enumerator :: DPARM_REFINE_TIME        = 15
     enumerator :: DPARM_A_NORM             = 16
     enumerator :: DPARM_COMPRESS_TOLERANCE = 17
     enumerator :: DPARM_COMPRESS_MIN_RATIO = 18
     enumerator :: DPARM_SIZE               = 18
  end enum

  ! enum task
  enum, bind(C)
     enumerator :: PastixTaskInit     = 0
     enumerator :: PastixTaskOrdering = 1
     enumerator :: PastixTaskSymbfact = 2
     enumerator :: PastixTaskAnalyze  = 3
     enumerator :: PastixTaskNumfact  = 4
     enumerator :: PastixTaskSolve    = 5
     enumerator :: PastixTaskRefine   = 6
     enumerator :: PastixTaskClean    = 7
  end enum

  ! enum verbose
  enum, bind(C)
     enumerator :: PastixVerboseNot = 0
     enumerator :: PastixVerboseNo  = 1
     enumerator :: PastixVerboseYes = 2
  end enum

  ! enum io
  enum, bind(C)
     enumerator :: PastixIONo        = 0
     enumerator :: PastixIOLoad      = 1
     enumerator :: PastixIOSave      = 2
     enumerator :: PastixIOLoadGraph = 4
     enumerator :: PastixIOSaveGraph = 8
     enumerator :: PastixIOLoadCSC   = 16
     enumerator :: PastixIOSaveCSC   = 32
  end enum

  ! enum fact_mode
  enum, bind(C)
     enumerator :: PastixFactModeLocal = 0
     enumerator :: PastixFactModeSchur = 1
     enumerator :: PastixFactModeBoth  = 2
  end enum

  ! enum solv_mode
  enum, bind(C)
     enumerator :: PastixSolvModeLocal     = 0
     enumerator :: PastixSolvModeInterface = 1
     enumerator :: PastixSolvModeSchur     = 2
  end enum

  ! enum refine
  enum, bind(C)
     enumerator :: PastixRefineGMRES    = 0
     enumerator :: PastixRefineCG       = 1
     enumerator :: PastixRefineSR       = 2
     enumerator :: PastixRefineBiCGSTAB = 3
  end enum

  ! enum factotype
  enum, bind(C)
     enumerator :: PastixFactPOTRF = 0
     enumerator :: PastixFactSYTRF = 1
     enumerator :: PastixFactGETRF = 2
     enumerator :: PastixFactPXTRF = 3
     enumerator :: PastixFactHETRF = 4
     enumerator :: PastixFactLLH   = 0
     enumerator :: PastixFactLDLT  = 1
     enumerator :: PastixFactLU    = 2
     enumerator :: PastixFactLLT   = 3
     enumerator :: PastixFactLDLH  = 4
  end enum

  ! enum scheduler
  enum, bind(C)
     enumerator :: PastixSchedSequential = 0
     enumerator :: PastixSchedStatic     = 1
     enumerator :: PastixSchedParsec     = 2
     enumerator :: PastixSchedStarPU     = 3
     enumerator :: PastixSchedDynamic    = 4
  end enum

  ! enum order
  enum, bind(C)
     enumerator :: PastixOrderScotch   = 0
     enumerator :: PastixOrderMetis    = 1
     enumerator :: PastixOrderPersonal = 2
     enumerator :: PastixOrderPtScotch = 3
     enumerator :: PastixOrderParMetis = 4
  end enum

  ! enum threadmode
  enum, bind(C)
     enumerator :: PastixThreadMultiple = 1
     enumerator :: PastixThreadFunneled = 2
  end enum

  ! enum error
  enum, bind(C)
     enumerator :: PASTIX_SUCCESS            = 0
     enumerator :: PASTIX_ERR_UNKNOWN        = 1
     enumerator :: PASTIX_ERR_ALLOC          = 2
     enumerator :: PASTIX_ERR_NOTIMPLEMENTED = 3
     enumerator :: PASTIX_ERR_OUTOFMEMORY    = 4
     enumerator :: PASTIX_ERR_THREAD         = 5
     enumerator :: PASTIX_ERR_INTERNAL       = 6
     enumerator :: PASTIX_ERR_BADPARAMETER   = 7
     enumerator :: PASTIX_ERR_FILE           = 8
     enumerator :: PASTIX_ERR_INTEGER_TYPE   = 9
     enumerator :: PASTIX_ERR_IO             = 10
     enumerator :: PASTIX_ERR_MPI            = 11
  end enum

  ! enum compress_when
  enum, bind(C)
     enumerator :: PastixCompressNever      = 0
     enumerator :: PastixCompressWhenBegin  = 1
     enumerator :: PastixCompressWhenEnd    = 2
     enumerator :: PastixCompressWhenDuring = 3
  end enum

  ! enum compress_method
  enum, bind(C)
     enumerator :: PastixCompressMethodSVD  = 0
     enumerator :: PastixCompressMethodRRQR = 1
  end enum

  ! enum compress_ortho
  enum, bind(C)
     enumerator :: PastixCompressOrthoCGS       = 0
     enumerator :: PastixCompressOrthoQR        = 1
     enumerator :: PastixCompressOrthoPartialQR = 2
  end enum

  ! enum layout
  enum, bind(C)
     enumerator :: PastixRowMajor = 101
     enumerator :: PastixColMajor = 102
  end enum

  ! enum trans
  enum, bind(C)
     enumerator :: PastixNoTrans   = 111
     enumerator :: PastixTrans     = 112
     enumerator :: PastixConjTrans = 113
  end enum

  ! enum mtxtype
  enum, bind(C)
     enumerator :: PastixGeneral   = PastixNoTrans
     enumerator :: PastixSymmetric = PastixTrans
     enumerator :: PastixHermitian = PastixConjTrans
  end enum

  ! enum uplo
  enum, bind(C)
     enumerator :: PastixUpper      = 121
     enumerator :: PastixLower      = 122
     enumerator :: PastixUpperLower = 123
  end enum

  ! enum coefside
  enum, bind(C)
     enumerator :: PastixLCoef  = 0
     enumerator :: PastixUCoef  = 1
     enumerator :: PastixLUCoef = 2
  end enum

  ! enum diag
  enum, bind(C)
     enumerator :: PastixNonUnit = 131
     enumerator :: PastixUnit    = 132
  end enum

  ! enum side
  enum, bind(C)
     enumerator :: PastixLeft  = 141
     enumerator :: PastixRight = 142
  end enum

  ! enum normtype
  enum, bind(C)
     enumerator :: PastixOneNorm       = 171
     enumerator :: PastixFrobeniusNorm = 174
     enumerator :: PastixInfNorm       = 175
     enumerator :: PastixMaxNorm       = 177
  end enum

  ! enum dir
  enum, bind(C)
     enumerator :: PastixDirForward  = 391
     enumerator :: PastixDirBackward = 392
  end enum

  integer, parameter :: pastix_int_t = PASTIX_INT_KIND

contains

  function pastix_getintsize()
    integer :: pastix_getintsize
    pastix_getintsize = kind(PASTIX_INT_KIND)
    return
  end function pastix_getintsize

end module pastix_enums
