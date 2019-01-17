/**
 *
 * @file old_api.h
 *
 * @copyright 2013-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * file to keep compatibility with older API.
 *
 * @version 6.0.1
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @date 2018-07-16
 *
 */
#ifndef _old_api_h_
#define _old_api_h_

#define pastix_float_t void

/* Error numbers, need to conserve it MURGE compliant */
#define NO_ERR             PASTIX_SUCCESS
#define UNKNOWN_ERR        PASTIX_ERR_UNKNOWN
#define ALLOC_ERR          PASTIX_ERR_ALLOC
#define NOTIMPLEMENTED_ERR PASTIX_ERR_NOTIMPLEMENTED
#define OUTOFMEMORY_ERR    PASTIX_ERR_OUTOFMEMORY
#define THREAD_ERR         PASTIX_ERR_THREAD
#define INTERNAL_ERR       PASTIX_ERR_INTERNAL
#define BADPARAMETER_ERR   PASTIX_ERR_BADPARAMETER
#define FILE_ERR           PASTIX_ERR_FILE
#define INTEGER_TYPE_ERR   PASTIX_ERR_INTEGER_TYPE
#define IO_ERR             PASTIX_ERR_IO
#define MPI_ERR            PASTIX_ERR_MPI

/* Removed from PaStiX 6.0.0 */
#define ASSERT_ERR         -1
#define BAD_DEFINE_ERR     -1
#define FLOAT_TYPE_ERR     -1
#define MATRIX_ERR         -1
#define STEP_ORDER_ERR     -1

#define API_NO  0
#define API_YES 1

/* Former IPARM values */
/*
 * Backward compatibility
 */
enum IPARM_ACCESS_DEPRECATED {
    IPARM_DEFAULT_ORDERING      = IPARM_ORDERING_DEFAULT,
    IPARM_ORDERING_SWITCH_LEVEL = IPARM_SCOTCH_SWITCH_LEVEL,
    IPARM_ORDERING_CMIN         = IPARM_SCOTCH_CMIN,
    IPARM_ORDERING_CMAX         = IPARM_SCOTCH_CMAX,
    IPARM_ORDERING_FRAT         = IPARM_SCOTCH_FRAT,
    IPARM_AMALGAMATION_LEVEL    = IPARM_AMALGAMATION_LVLCBLK,
    IPARM_ONLY_REFINE           = 0,
    IPARM_RHS_MAKING            = -1,
    IPARM_BINDTHRD              = -1,
    IPARM_SCHUR                 = -1,
    IPARM_ISOLATE_ZEROS         = -1
};

/* Former DPARM values */
#define DPARM_RAFF_TIME     DPARM_REFINE_TIME
#define DPARM_FACT_FLOPS    DPARM_FACT_THFLOPS

/* Former API values */

/* _POS_ 1 */
#define API_TASK_INIT       PastixTaskInit
#define API_TASK_SCOTCH     PastixTaskOrdering
#define API_TASK_ORDERING   PastixTaskOrdering
#define API_TASK_FAX        PastixTaskSymbfact
#define API_TASK_SYMBFACT   PastixTaskSymbfact
#define API_TASK_BLEND      PastixTaskAnalyze
#define API_TASK_ANALYSE    PastixTaskAnalyze
#define API_TASK_SOPALIN    PastixTaskNumfact
#define API_TASK_NUMFACT    PastixTaskNumfact
#define API_TASK_UPDOWN     PastixTaskSolve
#define API_TASK_SOLVE      PastixTaskSolve
#define API_TASK_REFINE     PastixTaskRefine
#define API_TASK_REFINEMENT PastixTaskRefine
#define API_TASK_CLEAN      PastixTaskClean

/* _POS_ 4 */
#define API_FACT_LLT        PastixFactPOTRF
#define API_FACT_LDLT       PastixFactSYTRF
#define API_FACT_LU         PastixFactGETRF
#define API_FACT_LDLH       PastixFactSYTRF

/* _POS_ 5 */
#define API_VERBOSE_NOT     PastixVerboseNot
#define API_VERBOSE_NO      PastixVerboseNo
#define API_VERBOSE_YES     PastixVerboseYes

/* _POS_ 6 */
#define API_IO_NO           PastixIONo
#define API_IO_LOAD         PastixIOLoad
#define API_IO_SAVE         PastixIOSave
#define API_IO_LOAD_GRAPH   PastixIOLoadGraph
#define API_IO_SAVE_GRAPH   PastixIOSaveGraph
#define API_IO_LOAD_CSC     PastixIOLoadCSC
#define API_IO_SAVE_CSC     PastixIOSaveCSC

/* _POS_ 7 */
#define API_RHS_B           SpmRhsRndB
#define API_RHS_1           SpmRhsOne
#define API_RHS_I           SpmRhsI
#define API_RHS_0           SpmRhsRndX

/* _POS_ 8 */
#define API_RAFF_GMRES      PastixGMRES
#define API_RAFF_GRAD       PastixCG
#define API_RAFF_PIVOT      PastixSR
#define API_RAFF_BICGSTAB   PastixBiGSTAB

/* _POS_ 11 */
#define API_ORDER_SCOTCH    PastixOrderScotch
#define API_ORDER_METIS     PastixOrderMetis
#define API_ORDER_PERSONAL  PastixOrderPersonal
#define API_ORDER_LOAD      PastixOrderLoad
#define API_ORDER_PTSCOTCH  PastixOrderPtScotch

/* _POS_ 61 */
#define API_REALSINGLE      PastixFloat
#define API_REALDOUBLE      PastixDouble
#define API_COMPLEXSINGLE   PastixComplex32
#define API_COMPLEXDOUBLE   PastixComplex64

/**
 * Some define for old pastix compatibility
 */
#define API_SYM_YES         PastixSymmetric
#define API_SYM_HER         PastixHermitian
#define API_SYM_NO          PastixGeneral

/**
 * @brief The list of matrix driver readers and generators from Spm
 */
#define pastix_driver_t spm_driver_t
#define PastixDriverRSA        SpmDriverRSA
#define PastixDriverHB         SpmDriverHB
#define PastixDriverIJV        SpmDriverIJV
#define PastixDriverMM         SpmDriverMM
#define PastixDriverLaplacian  SpmDriverLaplacian
#define PastixDriverXLaplacian SpmDriverXLaplacian
#define PastixDriverGraph      SpmDriverGraph
#define PastixDriverSPM        SpmDriverSPM

#endif /* _old_api_h_ */
