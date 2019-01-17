/**
 *
 * @file symbol_cost_perfs.c
 *
 * PaStiX symbol functions to compute the computational time induced by the chosen
 * symbolic structure with a given performance model.
 *
 * @copyright 1999-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.1
 * @author David Goudin
 * @author Francois Pellegrini
 * @author Mathieu Faverge
 * @author Pascal Henon
 * @author Pierre Ramet
 * @date 2018-07-16
 *
 * @addtogroup symbol_dev_cost
 * @{
 *
 **/
#include "common.h"
#include "symbol_cost.h"
#include "perf.h"

/**
 * @name Performance functions set
 * @{
 *
 */

/**
 * @brief Time model of the computation of the diagonal block
 * @param[in] N Size of the matrix block
 * @return Returns the time cost for the given arithmetic and factorization
 * @details Cholesky complex case
 */
static inline double
perfs_zpotrf_diag( pastix_int_t N ) {
    double total = PERF_POTRF(N);
    assert( N > 0 );
    return (total > 0.) ? total : 0.;
}

/**
 * @copydoc perfs_zpotrf_diag
 * @details Cholesky real case
 */
static inline double
perfs_dpotrf_diag( pastix_int_t N ) {
    double total = PERF_POTRF(N);
    assert( N > 0 );
    return (total > 0.) ? total : 0.;
}

/**
 * @copydoc perfs_zpotrf_diag
 * @details LU complex case
 */
static inline double
perfs_zgetrf_diag( pastix_int_t N ) {
    /* Approximate to 2 times potrf */
    return 2. * perfs_zpotrf_diag( N );
}

/**
 * @copydoc perfs_zpotrf_diag
 * @details LU real case
 */
static inline double
perfs_dgetrf_diag( pastix_int_t N ) {
    /* Approximate to 2 times potrf */
    return 2. * perfs_dpotrf_diag( N );
}

/**
 * @copydoc perfs_zpotrf_diag
 * @details LDL^t complex case
 */
static inline double
perfs_zsytrf_diag( pastix_int_t N ) {
    double total = PERF_SYTRF(N) + (double)N * PERF_COPY(N);
    assert( N > 0 );
    return (total > 0.) ? total : 0.;
}

/**
 * @copydoc perfs_zpotrf_diag
 * @details LDL^t real case
 */
static inline double
perfs_dsytrf_diag( pastix_int_t N ) {
    double total = PERF_SYTRF(N) + (double)N * PERF_COPY(N);
    assert( N > 0 );
    return (total > 0.) ? total : 0.;
}

/**
 * @brief Time performance model of the solve step
 * @param[in] M Number of rows of the B matrix in the TRSM, and size of the matrix A
 * @param[in] N Number of columns of the B matrix in the TRSM
 * @return Returns the time cost for the given arithmetic and factorization
 * @details Cholesky complex case
 */
static inline double
perfs_zpotrf_trsm( pastix_int_t M, pastix_int_t N ) {
    double total = PERF_TRSM( N, M );
    assert( (M > 0) && (N > 0) );
    return (total > 0.) ? total : 0.;
}

/**
 * @copydoc perfs_zpotrf_trsm
 * @details Cholesky real case
 */
static inline double
perfs_dpotrf_trsm( pastix_int_t M, pastix_int_t N ) {
    double total = PERF_TRSM( N, M );
    assert( (M > 0) && (N > 0) );
    return (total > 0.) ? total : 0.;
}

/**
 * @copydoc perfs_zpotrf_trsm
 * @details LU complex case
 */
static inline double
perfs_zgetrf_trsm( pastix_int_t M, pastix_int_t N ) {
    double total = 2. * PERF_TRSM( N, M );
    assert( (M > 0) && (N > 0) );
    return (total > 0.) ? total : 0.;
}

/**
 * @copydoc perfs_zpotrf_trsm
 * @details LU real case
 */
static inline double
perfs_dgetrf_trsm( pastix_int_t M, pastix_int_t N ) {
    double total = 2. * PERF_TRSM( N, M );
    assert( (M > 0) && (N > 0) );
    return (total > 0.) ? total : 0.;
}

/**
 * @copydoc perfs_zpotrf_trsm
 * @details LDL^t complex case
 */
static inline double
perfs_zsytrf_trsm( pastix_int_t M, pastix_int_t N ) {
    double total = PERF_TRSM( N, M )
        +  (double)N * (PERF_SCAL(M) + PERF_COPY(M));
    assert( (M > 0) && (N > 0) );
    return (total > 0.) ? total : 0.;
}

/**
 * @copydoc perfs_zpotrf_trsm
 * @details LDL^t real case
 */
static inline double
perfs_dsytrf_trsm( pastix_int_t M, pastix_int_t N ) {
    double total = PERF_TRSM( N, M )
        +  (double)N * (PERF_SCAL(M) + PERF_COPY(M));
    assert( (M > 0) && (N > 0) );
    return (total > 0.) ? total : 0.;
}

/**
 * @brief Time performance model of the update step per block (see sum2d())
 * @param[in] M Number of rows of the A and C matrices in the GEMM
 * @param[in] N Number of columns of the B and C matrices in the GEMM
 * @param[in] K Number of columns of the A matrix, and rows of the B matrix in the GEMM
 * @return Returns the time cost for the given arithmetic and factorization
 * @details Cholesky complex case
 */
static inline double
perfs_zpotrf_blkupdate( pastix_int_t M, pastix_int_t N, pastix_int_t K )
{
    double total = PERF_GEMM( M, N, K ) + PERF_GEAM( M, N );
    assert( (M > 0) && (N > 0) && (K > 0) );
    return (total > 0.) ? total : 0.;
}

/**
 * @copydoc flops_zpotrf_blkupdate
 * @details Cholesky real case
 */
static inline double
perfs_dpotrf_blkupdate( pastix_int_t M, pastix_int_t N, pastix_int_t K )
{
    double total = PERF_GEMM( M, N, K ) + PERF_GEAM( M, N );
    assert( (M > 0) && (N > 0) && (K > 0) );
    return (total > 0.) ? total : 0.;
}

/**
 * @copydoc flops_zpotrf_blkupdate
 * @details LU complex case
 */
static inline double
perfs_zgetrf_blkupdate( pastix_int_t M, pastix_int_t N, pastix_int_t K )
{
    double total = PERF_GEMM(     M, N, K ) + PERF_GEAM(     M, N )
        +          PERF_GEMM( (M-N), N, K ) + PERF_GEAM( (M-N), N );
    assert( (M > 0) && (N > 0) && (K > 0) );
    return (total > 0.) ? total : 0.;
}

/**
 * @copydoc flops_zpotrf_blkupdate
 * @details LU real case
 */
static inline double
perfs_dgetrf_blkupdate( pastix_int_t M, pastix_int_t N, pastix_int_t K )
{
    double total = PERF_GEMM(     M, N, K ) + PERF_GEAM(     M, N )
        +          PERF_GEMM( (M-N), N, K ) + PERF_GEAM( (M-N), N );
    assert( (M > 0) && (N > 0) && (K > 0) );
    return (total > 0.) ? total : 0.;
}

/**
 * @copydoc flops_zpotrf_blkupdate
 * @details LDL^t complex case
 */
static inline double
perfs_zsytrf_blkupdate( pastix_int_t M, pastix_int_t N, pastix_int_t K )
{
    double total = PERF_GEMM( M, N, K ) + PERF_GEAM( M, N );
    assert( (M > 0) && (N > 0) && (K > 0) );
    return (total > 0.) ? total : 0.;
}

/**
 * @copydoc flops_zpotrf_blkupdate
 * @details LDL^t real case
 */
static inline double
perfs_dsytrf_blkupdate( pastix_int_t M, pastix_int_t N, pastix_int_t K )
{
    double total = PERF_GEMM( M, N, K ) + PERF_GEAM( M, N );
    assert( (M > 0) && (N > 0) && (K > 0) );
    return (total > 0.) ? total : 0.;
}

/**
 * @}
 */

symbol_function_t perfstable[2][5] = {
    {
        {perfs_dpotrf_diag, perfs_dpotrf_trsm, NULL, perfs_dpotrf_blkupdate },
        {perfs_dsytrf_diag, perfs_dsytrf_trsm, NULL, perfs_dsytrf_blkupdate },
        {perfs_dgetrf_diag, perfs_dgetrf_trsm, NULL, perfs_dgetrf_blkupdate },
        {perfs_dpotrf_diag, perfs_dpotrf_trsm, NULL, perfs_dpotrf_blkupdate },
        {perfs_dsytrf_diag, perfs_dsytrf_trsm, NULL, perfs_dsytrf_blkupdate }
    },
    {
        {perfs_zpotrf_diag, perfs_zpotrf_trsm, NULL, perfs_zpotrf_blkupdate },
        {perfs_zsytrf_diag, perfs_zsytrf_trsm, NULL, perfs_zsytrf_blkupdate },
        {perfs_zgetrf_diag, perfs_zgetrf_trsm, NULL, perfs_zgetrf_blkupdate },
        {perfs_zpotrf_diag, perfs_zpotrf_trsm, NULL, perfs_zpotrf_blkupdate },
        {perfs_zsytrf_diag, perfs_zsytrf_trsm, NULL, perfs_zsytrf_blkupdate }
    }
};

/**
 * @}
 */
