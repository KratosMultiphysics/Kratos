/**
 *
 * @file pastix/api.h
 *
 * PaStiX API enums parameters.
 *
 * @copyright 2004-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.1
 * @author Xavier Lacoste
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @date 2018-07-16
 *
 * @addtogroup pastix_api
 * @{
 *
 **/
#ifndef _pastix_api_h_
#define _pastix_api_h_

/**
 * @brief Integer parameters
 */
typedef enum pastix_iparm_e {
    IPARM_VERBOSE,               /**< Verbose mode (@see pastix_verbose_t)                           Default: PastixVerboseNo          IN    */
    IPARM_IO_STRATEGY,           /**< IO strategy  (@see pastix_io_t)                                Default: PastixIONo               IN    */

    /* Stats */
    IPARM_NNZEROS,               /**< Number of nonzero entries in the factorized matrix             Default: -                        OUT   */
    IPARM_NNZEROS_BLOCK_LOCAL,   /**< Number of nonzero entries in the local block factorized matrix Default: -                        OUT   */
    IPARM_ALLOCATED_TERMS,       /**< Maximum memory allocated for matrix terms                      Default: -                        OUT   */
    IPARM_PRODUCE_STATS,         /**< Compute some statistiques (such as precision error)            Default: 0                        IN    */

    /* Scaling */
    IPARM_MC64,                  /**< MC64 operation                                                 Default: 0                        IN    */

    /* Ordering */
    IPARM_ORDERING,              /**< Choose ordering                                                Default: PastixOrderScotch        IN    */
    IPARM_ORDERING_DEFAULT,      /**< Use default ordering parameters with Scotch or Metis           Default: 1                        IN    */

    /* Subset for Scotch */
    IPARM_SCOTCH_SWITCH_LEVEL,   /**< Ordering switch level    (see Scotch Manual)                   Default: 120                      IN    */
    IPARM_SCOTCH_CMIN,           /**< Ordering cmin parameter  (see Scotch Manual)                   Default: 0                        IN    */
    IPARM_SCOTCH_CMAX,           /**< Ordering cmax parameter  (see Scotch Manual)                   Default: 100000                   IN    */
    IPARM_SCOTCH_FRAT,           /**< Ordering frat parameter  (see Scotch Manual)                   Default: 8                        IN    */

    /* Subset for Metis */
    IPARM_METIS_CTYPE,           /**< Metis parameters (see Metis Manual)                            Default: METIS_CTYPE_SHEM         IN    */
    IPARM_METIS_RTYPE,           /**< Metis parameters (see Metis Manual)                            Default: METIS_RTYPE_SEP1SIDED    IN    */
    IPARM_METIS_NO2HOP,          /**< Metis parameters (see Metis Manual)                            Default: 0                        IN    */
    IPARM_METIS_NSEPS,           /**< Metis parameters (see Metis Manual)                            Default: 1                        IN    */
    IPARM_METIS_NITER,           /**< Metis parameters (see Metis Manual)                            Default: 10                       IN    */
    IPARM_METIS_UFACTOR,         /**< Metis parameters (see Metis Manual)                            Default: 200                      IN    */
    IPARM_METIS_COMPRESS,        /**< Metis parameters (see Metis Manual)                            Default: 1                        IN    */
    IPARM_METIS_CCORDER,         /**< Metis parameters (see Metis Manual)                            Default: 0                        IN    */
    IPARM_METIS_PFACTOR,         /**< Metis parameters (see Metis Manual)                            Default: 0                        IN    */
    IPARM_METIS_SEED,            /**< Metis parameters (see Metis Manual)                            Default: 3452                     IN    */
    IPARM_METIS_DBGLVL,          /**< Metis parameters (see Metis Manual)                            Default: 0                        IN    */

    /* Symbolic Factorization */
    IPARM_SF_KASS,               /**< Force KASS instead of Fax to perform symbolic factorization    Default: 0                        IN    */
    IPARM_AMALGAMATION_LVLBLAS,  /**< Amalgamation level                                             Default: 5                        IN    */
    IPARM_AMALGAMATION_LVLCBLK,  /**< Amalgamation level                                             Default: 5                        IN    */

    /* Reordering */
    IPARM_REORDERING_SPLIT,      /**< Reordering split level                                         Default: 0                        IN    */
    IPARM_REORDERING_STOP,       /**< Reordering stop criteria                                       Default: PASTIX_INT_MAX           IN    */

    /* Analyze */
    IPARM_MIN_BLOCKSIZE,         /**< Minimum block size                                             Default: 160                      IN    */
    IPARM_MAX_BLOCKSIZE,         /**< Maximum block size                                             Default: 320                      IN    */
    IPARM_TASKS2D_LEVEL,         /**< 2D Distribution level (-1 for autolevel, 0 for 1D)             Default: -1                       IN    */
    IPARM_TASKS2D_WIDTH,         /**< Minimal width for 2D tasks with autolevel                      Default: IPARM_MIN_BLOCKSIZE      IN    */
    IPARM_ABS,                   /**< ABS level (Automatic Blocksize Splitting)                      Default: 0                        IN    */
    IPARM_ALLCAND,               /**< Allow all threads to be candidate in the proportional mapping  Default: 0                        IN    */

    /* Incomplete */
    IPARM_INCOMPLETE,            /**< Incomplete factorization                                       Default: 0                        IN    */
    IPARM_LEVEL_OF_FILL,         /**< Level of fill for incomplete factorization                     Default: 0                        IN    */

    /* Factorization */
    IPARM_FACTORIZATION,         /**< Factorization mode                                             Default: PastixFactLU             IN    */
    IPARM_STATIC_PIVOTING,       /**< Static pivoting                                                Default: -                        OUT   */
    IPARM_FREE_CSCUSER,          /**< Free user CSC                                                  Default: 0                        IN    */
    IPARM_SCHUR_FACT_MODE,       /**< Specify if the Schur is factorized (@see pastix_fact_mode_t)   Default: PastixFactModeLocal      IN    */

    /* Solve */
    IPARM_SCHUR_SOLV_MODE,       /**< Specify the solve parts to apply (@see pastix_solv_mode_t)     Default: PastixSolvModeLocal      IN    */
    IPARM_APPLYPERM_WS,          /**< Enable/disable extra workspace for a thread-safe swap          Default: 1                        IN    */

    /* Refinement */
    IPARM_REFINEMENT,            /**< Refinement mode                                                Default: PastixRefineGMRES        IN    */
    IPARM_NBITER,                /**< Number of iterations performed in refinement                   Default: -                        OUT   */
    IPARM_ITERMAX,               /**< Maximum iteration number for refinement                        Default: 250                      IN    */
    IPARM_GMRES_IM,              /**< GMRES restart parameter                                        Default: 25                       IN    */

    /* Context */
    IPARM_SCHEDULER,             /**< Scheduler mode                                                 Default: PastixSchedStatic        IN    */
    IPARM_THREAD_NBR,            /**< Number of threads per process (-1 for auto detect)             Default: -1                       IN    */
    IPARM_AUTOSPLIT_COMM,        /**< Automaticaly split communicator to have one MPI task by node   Default: 0                        IN    */

    /* GPU */
    IPARM_GPU_NBR,               /**< Number of GPU devices                                          Default: 0                        IN    */
    IPARM_GPU_MEMORY_PERCENTAGE, /**< Maximum percentage of the GPU memory used by the solver        Default: 95                       IN    */
    IPARM_GPU_MEMORY_BLOCK_SIZE, /**< Size of GPU memory pages (for PaRSEC runtime)                  Default: 32 * 1024                IN    */

    /* Compression */
    IPARM_COMPRESS_MIN_WIDTH,    /**< Minimum width to compress a supernode                          Default: 120                      IN    */
    IPARM_COMPRESS_MIN_HEIGHT,   /**< Minimum height to compress an off-diagonal block               Default: 20                       IN    */
    IPARM_COMPRESS_WHEN,         /**< When to compress a supernode                                   Default: PastixCompressNever      IN    */
    IPARM_COMPRESS_METHOD,       /**< Compression method (SVD/RRQR)                                  Default: PastixCompressMethodRRQR IN    */
    IPARM_COMPRESS_ORTHO,        /**< Orthogonalization method                                       Default: PastixCompressOrthoCGS   IN    */

    /* MPI modes */
    IPARM_THREAD_COMM_MODE,      /**< Threaded communication mode                                    Default: PastixThreadMultiple     IN    */

    /* Subset for old interface */
    IPARM_MODIFY_PARAMETER,      /**< Indicate if parameters have been set by user                   Default: 1                        IN    */
    IPARM_START_TASK,            /**< Indicate the first step to execute                             Default: PastixTaskOrdering       IN    */
    IPARM_END_TASK,              /**< Indicate the last step to execute                              Default: PastixTaskClean          IN    */
    IPARM_FLOAT,                 /**< Indicate the arithmetics                                       Default: PastixDouble             IN    */
    IPARM_MTX_TYPE,              /**< Indicate matrix format                                         Default: -1                       IN    */
    IPARM_DOF_NBR,               /**< Degree of freedom per node                                     Default: 1                        IN    */
    IPARM_SIZE
} pastix_iparm_t;


/**
 * @brief Float parameters
 */
typedef enum pastix_dparm_e {
    DPARM_FILL_IN,               /**< Maximum memory (-DMEMORY_USAGE)                   Default: -                OUT */
    DPARM_EPSILON_REFINEMENT,    /**< Epsilon for refinement                            Default: -1.              IN  */
    DPARM_RELATIVE_ERROR,        /**< Relative backward error                           Default: -                OUT */
    DPARM_EPSILON_MAGN_CTRL,     /**< Epsilon for magnitude control                     Default: 0.               IN  */
    DPARM_ANALYZE_TIME,          /**< Time for Analyse step (wallclock)                 Default: -                OUT */
    DPARM_PRED_FACT_TIME,        /**< Predicted factorization time                      Default: -                OUT */
    DPARM_FACT_TIME,             /**< Time for Numerical Factorization step (wallclock) Default: -                OUT */
    DPARM_SOLV_TIME,             /**< Time for Solve step (wallclock)                   Default: -                OUT */
    DPARM_FACT_FLOPS,            /**< Factorization GFlops/s                            Default: -                OUT */
    DPARM_FACT_THFLOPS,          /**< Factorization theoretical Flops                   Default: -                OUT */
    DPARM_FACT_RLFLOPS,          /**< Factorization performed Flops                     Default: -                OUT */
    DPARM_SOLV_FLOPS,            /**< Solve GFlops/s                                    Default: -                OUT */
    DPARM_SOLV_THFLOPS,          /**< Solve theoretical Flops                           Default: -                OUT */
    DPARM_SOLV_RLFLOPS,          /**< Solve performed Flops                             Default: -                OUT */
    DPARM_REFINE_TIME,           /**< Time for Refinement step (wallclock)              Default: -                OUT */
    DPARM_A_NORM,                /**< ||A||_f norm                                      Default: -                OUT */
    DPARM_COMPRESS_TOLERANCE,    /**< Tolerance for low-rank kernels                    Default: 0.01             IN  */
    DPARM_COMPRESS_MIN_RATIO,    /**< Min ratio for rank w.r.t. strict rank             Default: 1.0              IN  */
    DPARM_SIZE
} pastix_dparm_t;

/**
 * @brief Main steps for the pastix() interface.
 *
 * Those enums are used of the IPARM_START_TASK and IPARM_END_TASK parameters
 * that configure the pastix() call.
 */
typedef enum pastix_task_e {
    PastixTaskInit       = 0, /**< Startup the library          */
    PastixTaskOrdering   = 1, /**< Ordering                     */
    PastixTaskSymbfact   = 2, /**< Symbolic factorization       */
    PastixTaskAnalyze    = 3, /**< Tasks mapping and scheduling */
    PastixTaskNumfact    = 4, /**< Numerical factorization      */
    PastixTaskSolve      = 5, /**< Numerical solve              */
    PastixTaskRefine     = 6, /**< Numerical refinement         */
    PastixTaskClean      = 7  /**< Clean                        */
} pastix_task_t;

/**
 * @brief Verbose modes
 */
typedef enum pastix_verbose_e {
    PastixVerboseNot = 0, /**< Nothing  */
    PastixVerboseNo  = 1, /**< Default  */
    PastixVerboseYes = 2  /**< Extended */
} pastix_verbose_t;

/**
 * @brief IO strategy for graph and ordering
 */
typedef enum pastix_io_e {
    PastixIONo         = 0, /**< No output or input */
    PastixIOLoad       = 1, /**< Load ordering and symbol matrix instead of applying symbolic factorisation step */
    PastixIOSave       = 2, /**< Save ordering and symbol matrix after symbolic factorisation step */
    PastixIOLoadGraph  = 4, /**< Load graph  during ordering step */
    PastixIOSaveGraph  = 8, /**< Save graph  during ordering step */
    PastixIOLoadCSC    = 16,/**< Load CSC(d) during ordering step */
    PastixIOSaveCSC    = 32 /**< Save CSC(d) during ordering step */
} pastix_io_t;

/**
 * @brief Factorization Schur modes
 *
 * Describe which part of the matrix is factorized or not
 *
 */
typedef enum pastix_fact_mode_e {
    PastixFactModeLocal   = 0,
    PastixFactModeSchur   = 1,
    PastixFactModeBoth    = 2
} pastix_fact_mode_t;

/**
 * @brief Solve Schur modes
 *
 * Describe which part of the solve is applied with the matrix
 *
 * \f[ A = \left( \begin{array}{cc}
 *             L_{11}U_{11} & U_{12} \\
 *             L_{21}       & S_{22} \end{array} \right) \f]
 *
 * For the lower part (and symmetrically for upper part):
 *   -# Solve \f[ L_{11} * x_{11} = b_{11} \f]
 *   -# Apply the update \f[ b_{22} = b_{22} - L_{21} * b_{11} \f]
 *   -# Solve the lower part of \f[ S_{22} * x_{22} = b_{22} \f] if S22 has been previously factorized.
 *
 * PastixSolvModeLocal applies only the step 1.
 * PastixSolvModeInterface applies steps 1 and 2.
 * PastixSolvModeSchur applies all steps.
 *
 */
typedef enum pastix_solv_mode_e {
    PastixSolvModeLocal     = 0,
    PastixSolvModeInterface = 1,
    PastixSolvModeSchur     = 2
} pastix_solv_mode_t;

/**
 * @brief Iterative refinement algorithms
 */
typedef enum pastix_refine_e {
    PastixRefineGMRES,   /**< GMRES              */
    PastixRefineCG,      /**< Conjugate Gradiant */
    PastixRefineSR,      /**< Simple refinement  */
    PastixRefineBiCGSTAB /**< BiCGStab           */
} pastix_refine_t;

/**
 * @brief Arithmetic types.
 *
 * This describes the different arithmetics that can be stored in a sparse matrix.
 * @remark The values start at 2 for compatibility purpose with PLASMA and
 * DPLASMA libraries, and they match the ones used in spm.
 *
 * @sa spm_coeftype_t
 *
 * @{
 */
#define pastix_coeftype_t spm_coeftype_t
#define PastixPattern   SpmPattern
#define PastixFloat     SpmFloat
#define PastixDouble    SpmDouble
#define PastixComplex32 SpmComplex32
#define PastixComplex64 SpmComplex64
/**
 * @}
 */

/**
 * @brief Factorization algorithms available for IPARM_FACTORIZATION parameter
 */
typedef enum pastix_factotype_e {
    PastixFactPOTRF = 0, /**< Cholesky factorization                   */
    PastixFactSYTRF = 1, /**< LDL^t factorization                      */
    PastixFactGETRF = 2, /**< LU factorization                         */
    PastixFactPXTRF = 3, /**< LL^t factorization for complex matrices  */
    PastixFactHETRF = 4, /**< LDL^h factorization for complex matrices */

    PastixFactLLH  = 0, /**< LL^h factorization for complex matrices  */
    PastixFactLDLT = 1, /**< LDL^t factorization                      */
    PastixFactLU   = 2, /**< LU factorization                         */
    PastixFactLLT  = 3, /**< LL^t factorization                       */
    PastixFactLDLH = 4, /**< LDL^h factorization for complex matrices */
} pastix_factotype_t;

/**
 * @brief Scheduler
 */
typedef enum pastix_scheduler_e {
    PastixSchedSequential = 0, /**< Sequential                           */
    PastixSchedStatic     = 1, /**< Shared memory with static scheduler  */
    PastixSchedParsec     = 2, /**< PaRSEC scheduler                     */
    PastixSchedStarPU     = 3, /**< StarPU scheduler                     */
    PastixSchedDynamic    = 4, /**< Shared memory with dynamic scheduler */
} pastix_scheduler_t;

/**
 * @brief Ordering strategy
 */
enum pastix_order_e {
    PastixOrderScotch,   /**< Use Scotch ordering                         */
    PastixOrderMetis,    /**< Use Metis ordering                          */
    PastixOrderPersonal, /**< Apply user's permutation, or load from file */
    PastixOrderPtScotch, /**< Use Pt-Scotch ordering                      */
    PastixOrderParMetis  /**< Use ParMetis ordering                       */
};

#if defined(PASTIX_WITH_MPI)
/**
 * @brief MPI thread mode
 */
typedef enum pastix_threadmode_e {
    PastixThreadMultiple = 1, /**< All threads communicate              */
    PastixThreadFunneled = 2  /**< One thread perform all the MPI Calls */
} pastix_threadmode_t;
#endif /* defined(PASTIX_WITH_MPI) */

/**
 * @brief Error codes
 */
typedef enum pastix_error_e {
    PASTIX_SUCCESS            = 0,  /**< No error                     */
    PASTIX_ERR_UNKNOWN        = 1,  /**< Unknown error                */
    PASTIX_ERR_ALLOC          = 2,  /**< Allocation error             */
    PASTIX_ERR_NOTIMPLEMENTED = 3,  /**< Not implemented feature      */
    PASTIX_ERR_OUTOFMEMORY    = 4,  /**< Not enough memory            */
    PASTIX_ERR_THREAD         = 5,  /**< Error with threads           */
    PASTIX_ERR_INTERNAL       = 6,  /**< Internal error               */
    PASTIX_ERR_BADPARAMETER   = 7,  /**< Bad parameters given         */
    PASTIX_ERR_FILE           = 8,  /**< Error in In/Out operations   */
    PASTIX_ERR_INTEGER_TYPE   = 9,  /**< Error with integer types     */
    PASTIX_ERR_IO             = 10, /**< Error with input/output      */
    PASTIX_ERR_MPI            = 11  /**< Error with MPI calls         */
} pastix_error_t;

/**
 * @brief Compression strategy available for IPARM_COMPRESS_WHEN parameter
 */
typedef enum pastix_compress_when_e {
    PastixCompressNever,
    PastixCompressWhenBegin,
    PastixCompressWhenEnd,
    PastixCompressWhenDuring
} pastix_compress_when_t;

/**
 * @brief Compression method available for IPARM_COMPRESS_METHOD parameter
 */
typedef enum pastix_compress_method_e {
    PastixCompressMethodSVD,
    PastixCompressMethodRRQR
} pastix_compress_method_t;

/**
 * @brief Orthogonalization method available for IPARM_COMPRESS_ORTHO parameter
 */
typedef enum pastix_compress_ortho_e {
    PastixCompressOrthoCGS,
    PastixCompressOrthoQR,
    PastixCompressOrthoPartialQR,
} pastix_compress_ortho_t;

/**
 *
 * @name Constants compatible with CBLAS & LAPACK & PLASMA
 * @{
 *    The naming and numbering of the following constants is consistent with:
 *
 *       - CBLAS from Netlib (http://www.netlib.org/blas/blast-forum/cblas.tgz)
 *       - C Interface to LAPACK from Netlib (http://www.netlib.org/lapack/lapwrapc/)
 *       - Plasma (http://icl.cs.utk.edu/plasma/index.html)
 *
 */

/**
 * @brief Direction of the matrix storage
 */
typedef enum pastix_layout_e {
    PastixRowMajor  = 101, /**< Storage in row major order    */
    PastixColMajor  = 102  /**< Storage in column major order */
} pastix_layout_t;

/**
 * @brief Transpostion
 */
typedef enum pastix_trans_e {
    PastixNoTrans   = 111, /**< Use A         */
    PastixTrans     = 112, /**< Use A^t       */
    PastixConjTrans = 113  /**< Use conj(A^t) */
} pastix_trans_t;

/**
 * @brief Matrix symmetry type property.
 * @remark Must match transposition.
 */
typedef enum pastix_mtxtype_e {
    PastixGeneral   = PastixNoTrans,    /**< The matrix is general   */
    PastixSymmetric = PastixTrans,      /**< The matrix is symmetric */
    PastixHermitian = PastixConjTrans   /**< The matrix is hermitian */
} pastix_mtxtype_t;

/**
 * @brief Upper/Lower part
 */
typedef enum pastix_uplo_e {
    PastixUpper      = 121, /**< Use lower triangle of A */
    PastixLower      = 122, /**< Use upper triangle of A */
    PastixUpperLower = 123  /**< Use the full A          */
} pastix_uplo_t;

/**
 * @brief Data blocks used in the kernel
 */
typedef enum pastix_coefside_e {
    PastixLCoef  = 0, /**< Coefficients of the lower triangular L are used         */
    PastixUCoef  = 1, /**< Coefficients of the upper triangular U are used         */
    PastixLUCoef = 2  /**< Coefficients of the upper/lower triangular U/L are used */
} pastix_coefside_t;

/**
 * @brief Diagonal
 */
typedef enum pastix_diag_e {
    PastixNonUnit = 131, /**< Diagonal is non unitary */
    PastixUnit    = 132  /**< Diagonal is unitary     */
} pastix_diag_t;

/**
 * @brief Side of the operation
 */
typedef enum pastix_side_e {
    PastixLeft  = 141, /**< Apply operator on the left  */
    PastixRight = 142  /**< Apply operator on the right */
} pastix_side_t;

/**
 * @brief Norms
 */
typedef enum pastix_normtype_e {
    PastixOneNorm       = 171, /**< One norm:       max_j( sum_i( |a_{ij}| ) )   */
    PastixFrobeniusNorm = 174, /**< Frobenius norm: sqrt( sum_{i,j} (a_{ij}^2) ) */
    PastixInfNorm       = 175, /**< Inifinite norm: max_i( sum_j( |a_{ij}| ) )   */
    PastixMaxNorm       = 177  /**< Inifinite norm: max_{i,j}( | a_{ij} | )      */
} pastix_normtype_t;

/**
 * @brief Direction
 */
typedef enum pastix_dir_e {
    PastixDirForward  = 391, /**< Forward direction   */
    PastixDirBackward = 392, /**< Backward direction  */
} pastix_dir_t;

/**
 * @}
 */

#endif /* _pastix_api_h_ */

/**
 * @}
 */
