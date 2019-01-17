/**
 *
 * @file symbol_cost_flops.c
 *
 * PaStiX symbol functions to compute the number of flops induced by the chosen
 * symbolic structure.
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
#include "flops.h"

/**
 * @name Flops functions set
 * @{
 *
 */

/**
 * @brief Computations flops of diagonal blocks
 * @param[in] N Size of the matrix block
 * @return Returns the number of flops for the given arithmetic and factorization
 * @details Cholesky complex case
 */
static inline double
flops_zpotrf_diag( pastix_int_t N ) {
    return FLOPS_ZPOTRF( N );
}

/**
 * @copydoc flops_zpotrf_diag
 * @details Cholesky real case
 */
static inline double
flops_dpotrf_diag( pastix_int_t N ) {
    return FLOPS_DPOTRF( N );
}

/**
 * @copydoc flops_zpotrf_diag
 * @details LU complex case
 */
static double
flops_zgetrf_diag( pastix_int_t N ) {
    return FLOPS_ZGETRF( N, N );
}

/**
 * @copydoc flops_zpotrf_diag
 * @details LU real case
 */
static inline double
flops_dgetrf_diag( pastix_int_t N ) {
    return FLOPS_DGETRF( N, N );
}

/**
 * @copydoc flops_zpotrf_diag
 * @details LDL^t complex case
 */
static inline double
flops_zsytrf_diag( pastix_int_t N ) {
    return FLOPS_ZSYTRF( N );
}

/**
 * @copydoc flops_zpotrf_diag
 * @details LDL^t real case
 */
static inline double
flops_dsytrf_diag( pastix_int_t N ) {
    return FLOPS_DSYTRF( N );
}

/**
 * @brief Computations flops of the solve step
 * @param[in] M Number of rows of the B matrix in the TRSM, and size of the matrix A
 * @param[in] N Number of columns of the B matrix in the TRSM
 * @return Returns the number of flops for the given arithmetic and factorization
 * @details Cholesky complex case
 */
static inline double
flops_zpotrf_trsm( pastix_int_t M, pastix_int_t N ) {
    return FLOPS_ZTRSM( PastixRight, M, N );
}

/**
 * @copydoc flops_zpotrf_trsm
 * @details Cholesky real case
 */
static inline double
flops_dpotrf_trsm( pastix_int_t M, pastix_int_t N ) {
    return FLOPS_DTRSM( PastixRight, M, N );
}

/**
 * @copydoc flops_zpotrf_trsm
 * @details LU complex case
 */
static double
flops_zgetrf_trsm( pastix_int_t M, pastix_int_t N ) {
    return 2. * FLOPS_ZTRSM( PastixRight, M, N );
}

/**
 * @copydoc flops_zpotrf_trsm
 * @details LU real case
 */
static inline double
flops_dgetrf_trsm( pastix_int_t M, pastix_int_t N ) {
    return 2. * FLOPS_DTRSM( PastixRight, M, N );
}

/**
 * @copydoc flops_zpotrf_trsm
 * @details LDL^t complex case
 */
static inline double
flops_zsytrf_trsm( pastix_int_t M, pastix_int_t N ) {
    return FLOPS_ZTRSM( PastixRight, M, N ) + 6. * (double)N * (double)M;
}

/**
 * @copydoc flops_zpotrf_trsm
 * @details LDL^t real case
 */
static inline double
flops_dsytrf_trsm( pastix_int_t M, pastix_int_t N ) {
    return FLOPS_DTRSM( PastixRight, M, N ) + (double)N * (double)M;
}

/**
 * @brief Theroretical computation flops of the update step per coumn block (see sum1d())
 * @param[in] K Number of columns of A, and rows of B in the GEMM operation
 * @param[in] M Dimension of all other sizes in the GEMM
 * @return Returns the number of flops for the given arithmetic and factorization
 * @details Cholesky complex case
 */
static inline double
flops_zpotrf_update( pastix_int_t K, pastix_int_t M ) {
    return FLOPS_ZHERK( K, M );
}

/**
 * @copydoc flops_zpotrf_update
 * @details Cholesky real case
 */
static inline double
flops_dpotrf_update( pastix_int_t K, pastix_int_t M ) {
    return FLOPS_DSYRK( K, M );
}

/**
 * @copydoc flops_zpotrf_update
 * @details LU complex case
 */
static double
flops_zgetrf_update( pastix_int_t K, pastix_int_t M ) {
    return FLOPS_ZGEMM( M, M, K );
}

/**
 * @copydoc flops_zpotrf_update
 * @details LU real case
 */
static inline double
flops_dgetrf_update( pastix_int_t K, pastix_int_t M ) {
    return FLOPS_DGEMM( M, M, K );
}

/**
 * @copydoc flops_zpotrf_update
 * @details LDL^t complex case
 */
static inline double
flops_zsytrf_update( pastix_int_t K, pastix_int_t M ) {
    return FLOPS_ZSYRK( K, M ) + 6. * (double)M * (double)M;
}

/**
 * @copydoc flops_zpotrf_update
 * @details LDL^t real case
 */
static inline double
flops_dsytrf_update( pastix_int_t K, pastix_int_t M ) {
    return FLOPS_DSYRK( K, M ) + (double)M * (double)M;
}

/**
 * @brief Computes the theoretical number of flops of the update step per block (see sum2d())
 * @param[in] M Number of rows of the A and C matrices in the GEMM
 * @param[in] N Number of columns of the B and C matrices in the GEMM
 * @param[in] K Number of columns of the A matrix, and rows of the B matrix in the GEMM
 * @return Returns the number of flops for the given arithmetic and factorization
 * @details Cholesky complex case
 */
static inline double
flops_zpotrf_blkupdate( pastix_int_t M, pastix_int_t N, pastix_int_t K )
{
    return FLOPS_ZGEMM( M, N, K ) + 2. * (double)M * (double)N;
}

/**
 * @copydoc flops_zpotrf_blkupdate
 * @details Cholesky real case
 */
static inline double
flops_dpotrf_blkupdate( pastix_int_t M, pastix_int_t N, pastix_int_t K )
{
    return FLOPS_DGEMM( M, N, K ) + (double)M * (double)N;
}

/**
 * @copydoc flops_zpotrf_blkupdate
 * @details LU complex case
 */
static double
flops_zgetrf_blkupdate( pastix_int_t M, pastix_int_t N, pastix_int_t K )
{
    return FLOPS_ZGEMM( M, N, K ) + FLOPS_ZGEMM( M-N, N, K )
        + 2. * (double)M * (double)N + 2. * (double)(M-N) * (double)(N); /* Add step */
}

/**
 * @copydoc flops_zpotrf_blkupdate
 * @details LU real case
 */
static inline double
flops_dgetrf_blkupdate( pastix_int_t M, pastix_int_t N, pastix_int_t K )
{
    return FLOPS_DGEMM( M, N, K ) + FLOPS_DGEMM( M-N, N, K )
        + (double)M * (double)N + (double)(M-N) * (double)(N); /* Add step */
}

/**
 * @copydoc flops_zpotrf_blkupdate
 * @details LDL^t complex case
 */
static inline double
flops_zsytrf_blkupdate( pastix_int_t M, pastix_int_t N, pastix_int_t K )
{
    /* If we consider that we stored the D * A^t somewhere */
#if 0
    return FLOPS_ZGEMM( M, N, K )
        + 2. * (double)M * (double)N;
#else
    /* If not, as it is the case in the runtime */
    return FLOPS_ZGEMM( M, N, K )
        + 2. * (double)M * (double)N   /* Add step   */
        + 6. * (double)M * (double)N;  /* Scale step */
#endif
}

/**
 * @copydoc flops_zpotrf_blkupdate
 * @details LDL^t real case
 */
static inline double
flops_dsytrf_blkupdate( pastix_int_t M, pastix_int_t N, pastix_int_t K )
{
    /* If we consider that we stored the D * A^t somewhere */
#if 0
    return FLOPS_DGEMM( M, N, K )
        + (double)M * (double)N;
#else
    /* If not, as it is the case in the runtime */
    return FLOPS_DGEMM( M, N, K )
        + (double)M * (double)N   /* Add step   */
        + (double)M * (double)N;  /* Scale step */
#endif
}

/**
 * @}
 */

symbol_function_t flopstable[2][5] = {
    {
        {flops_dpotrf_diag, flops_dpotrf_trsm, flops_dpotrf_update, flops_dpotrf_blkupdate },
        {flops_dsytrf_diag, flops_dsytrf_trsm, flops_dsytrf_update, flops_dsytrf_blkupdate },
        {flops_dgetrf_diag, flops_dgetrf_trsm, flops_dgetrf_update, flops_dgetrf_blkupdate },
        {flops_dpotrf_diag, flops_dpotrf_trsm, flops_dpotrf_update, flops_dpotrf_blkupdate },
        {flops_dsytrf_diag, flops_dsytrf_trsm, flops_dsytrf_update, flops_dsytrf_blkupdate }
    },
    {
        {flops_zpotrf_diag, flops_zpotrf_trsm, flops_zpotrf_update, flops_zpotrf_blkupdate },
        {flops_zsytrf_diag, flops_zsytrf_trsm, flops_zsytrf_update, flops_zsytrf_blkupdate },
        {flops_zgetrf_diag, flops_zgetrf_trsm, flops_zgetrf_update, flops_zgetrf_blkupdate },
        {flops_zpotrf_diag, flops_zpotrf_trsm, flops_zpotrf_update, flops_zpotrf_blkupdate },
        {flops_zsytrf_diag, flops_zsytrf_trsm, flops_zsytrf_update, flops_zsytrf_blkupdate }
    }
};

/**
 * @}
 */
