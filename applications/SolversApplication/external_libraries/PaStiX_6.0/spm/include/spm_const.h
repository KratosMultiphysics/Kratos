/**
 *
 * @file spm/api.h
 *
 * Spm API enums parameters.
 *
 * @copyright 2004-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.0
 * @author Xavier Lacoste
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @author Matthieu Kuhn
 * @date 2013-06-24
 *
 * @addtogroup spm_api
 * @{
 *
 **/
#ifndef _spm_api_h_
#define _spm_api_h_

/********************************************************************
 * CBLAS value address
 */
#ifndef CBLAS_SADDR
#define CBLAS_SADDR( a_ ) (&(a_))
#endif

/**
 * @brief Verbose modes
 */
typedef enum spm_verbose_e {
    SpmVerboseNot = 0, /**< Nothing  */
    SpmVerboseNo  = 1, /**< Default  */
    SpmVerboseYes = 2  /**< Extended */
} spm_verbose_t;

/**
 * @brief Arithmetic types.
 *
 * This describes the different arithmetics that can be stored in a sparse matrix.
 * @remark The values start at 2 for compatibility purpose with PLASMA and
 * DPLASMA libraries.
 */
typedef enum spm_coeftype_e {
    SpmPattern   = 0, /**< Pattern only, no values are stored */
    SpmFloat     = 2, /**< Single precision real              */
    SpmDouble    = 3, /**< Double precision real              */
    SpmComplex32 = 4, /**< Single precision complex           */
    SpmComplex64 = 5  /**< Double precision complex           */
} spm_coeftype_t;

/**
 * @brief Sparse matrix format
 */
typedef enum spm_fmttype_e {
    SpmCSC, /**< Compressed sparse column */
    SpmCSR, /**< Compressed sparse row    */
    SpmIJV  /**< Coordinates              */
} spm_fmttype_t;

/**
 * @brief Error codes
 */
typedef enum spm_error_e {
    SPM_SUCCESS            = 0,  /**< No error                     */
    SPM_ERR_UNKNOWN        = 1,  /**< Unknown error                */
    SPM_ERR_ALLOC          = 2,  /**< Allocation error             */
    SPM_ERR_NOTIMPLEMENTED = 3,  /**< Not implemented feature      */
    SPM_ERR_OUTOFMEMORY    = 4,  /**< Not enough memory            */
    SPM_ERR_THREAD         = 5,  /**< Error with threads           */
    SPM_ERR_INTERNAL       = 6,  /**< Internal error               */
    SPM_ERR_BADPARAMETER   = 7,  /**< Bad parameters given         */
    SPM_ERR_FILE           = 8,  /**< Error in In/Out operations   */
    SPM_ERR_INTEGER_TYPE   = 9,  /**< Error with integer types     */
    SPM_ERR_IO             = 10, /**< Error with input/output      */
    SPM_ERR_MPI            = 11  /**< Error with MPI calls         */
} spm_error_t;

/**
 * @brief The list of matrix driver readers and generators
 */
typedef enum spm_driver_e {
    SpmDriverRSA,        /**< RSA Fortran driver (deprecated)                 */
    SpmDriverHB,         /**< Harwell Boeing driver                           */
    SpmDriverIJV,        /**< IJV Coordinate driver                           */
    SpmDriverMM,         /**< Matrix Market C driver                          */
    SpmDriverLaplacian,  /**< 3, 5, or 7 points Laplacian stencil generator   */
    SpmDriverXLaplacian, /**< 15-points Laplacian stencil generator           */
    SpmDriverGraph,      /**< Scotch Graph driver                             */
    SpmDriverSPM,        /**< SPM matrix driver                               */
    /* SpmDriverDMM,        /\**< Distributed Matrix Market driver                *\/ */
    /* SpmDriverCSCD,       /\**< CSC distributed driver                          *\/ */
    /* SpmDriverPetscS,     /\**< Petsc Symmetric driver                          *\/ */
    /* SpmDriverPetscU,     /\**< Pestc Unssymmetric driver                       *\/ */
    /* SpmDriverPetscH,     /\**< Pestc Hermitian driver                          *\/ */
    /* SpmDriverCCC,        /\**< Not supported yet *\/ */
    /* SpmDriverRCC,        /\**< Not supported yet *\/ */
    /* SpmDriverOlaf,       /\**< Not supported yet *\/ */
    /* SpmDriverPeer,       /\**< Not supported yet *\/ */
    /* SpmDriverBRGM,       /\**< Not supported yet *\/ */
    /* SpmDriverBRGMD,      /\**< Not supported yet *\/ */
} spm_driver_t;

/**
 * @brief How to generate RHS
 */
typedef enum spm_rhstype_e {
    SpmRhsOne,
    SpmRhsI,
    SpmRhsRndX,
    SpmRhsRndB
} spm_rhstype_t;

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
typedef enum spm_layout_e {
    SpmRowMajor  = 101, /**< Storage in row major order    */
    SpmColMajor  = 102  /**< Storage in column major order */
} spm_layout_t;

/**
 * @brief Transpostion
 */
typedef enum spm_trans_e {
    SpmNoTrans   = 111, /**< Use A         */
    SpmTrans     = 112, /**< Use A^t       */
    SpmConjTrans = 113  /**< Use conj(A^t) */
} spm_trans_t;

/**
 * @brief Matrix symmetry type property.
 * @remark Must match transposition.
 */
typedef enum spm_mtxtype_e {
    SpmGeneral   = SpmNoTrans,    /**< The matrix is general   */
    SpmSymmetric = SpmTrans,      /**< The matrix is symmetric */
    SpmHermitian = SpmConjTrans   /**< The matrix is hermitian */
} spm_mtxtype_t;

/**
 * @brief Upper/Lower part
 */
typedef enum spm_uplo_e {
    SpmUpper      = 121, /**< Use lower triangle of A */
    SpmLower      = 122, /**< Use upper triangle of A */
    SpmUpperLower = 123  /**< Use the full A          */
} spm_uplo_t;

/**
 * @brief Diagonal
 */
typedef enum spm_diag_e {
    SpmNonUnit = 131, /**< Diagonal is non unitary */
    SpmUnit    = 132  /**< Diagonal is unitary     */
} spm_diag_t;

/**
 * @brief Side of the operation
 */
typedef enum spm_side_e {
    SpmLeft  = 141, /**< Apply operator on the left  */
    SpmRight = 142  /**< Apply operator on the right */
} spm_side_t;

/**
 * @brief Norms
 */
typedef enum spm_normtype_e {
    SpmOneNorm       = 171, /**< One norm:       max_j( sum_i( |a_{ij}| ) )   */
    SpmFrobeniusNorm = 174, /**< Frobenius norm: sqrt( sum_{i,j} (a_{ij}^2) ) */
    SpmInfNorm       = 175, /**< Inifinite norm: max_i( sum_j( |a_{ij}| ) )   */
    SpmMaxNorm       = 177  /**< Inifinite norm: max_{i,j}( | a_{ij} | )      */
} spm_normtype_t;

/**
 * @brief Direction
 */
typedef enum spm_dir_e {
    SpmDirForward  = 391, /**< Forward direction   */
    SpmDirBackward = 392, /**< Backward direction  */
} spm_dir_t;

/**
 * @}
 */

#endif /* _spm_api_h_ */

/**
 * @}
 */
