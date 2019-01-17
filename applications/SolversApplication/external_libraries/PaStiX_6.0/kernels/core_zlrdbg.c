/**
 *
 * @file core_zlrdbg.c
 *
 * PaStiX low-rank kernel debug routines that may be call within gdb.
 *
 * @copyright 2016-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.1
 * @author Gregoire Pichon
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @date 2018-07-16
 * @precisions normal z -> c d s
 *
 **/
#include "common.h"
#include <cblas.h>
#include <lapacke.h>
#include "pastix_zlrcores.h"
#include "z_nan_check.h"

/**
 *******************************************************************************
 *
 * @brief Print the svd values of the given matrix.
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the matrix A.
 *
 * @param[in] N
 *          The number of columns of the matrix A.
 *
 * @param[in] A
 *          The matrix A to study of size lda-by-N
 *
 * @param[in] lda
 *          The leading dimension of the matrix A. lda = max( 1, M )
 *
 *******************************************************************************/
void
core_zlrdbg_printsvd( pastix_int_t              M,
                      pastix_int_t              N,
                      const pastix_complex64_t *A,
                      pastix_int_t              lda )
{
    pastix_int_t i;
    pastix_int_t minMN = pastix_imin( M, N );
    size_t lrwork =  2 * minMN;
    size_t lzwork =  M * N;
    pastix_complex64_t *W;
    double *s, *superb;

    W = malloc( lzwork * sizeof(pastix_complex64_t) + lrwork * sizeof(double) );
    s = (double*)(W + M*N);
    superb = s + minMN;

    LAPACKE_zlacpy_work( LAPACK_COL_MAJOR, 'A', M, N, A, lda, W, M );
    LAPACKE_zgesvd(LAPACK_COL_MAJOR, 'N', 'N', M, N, W, M, s, NULL, 1, NULL, 1, superb );

    for(i=0; i<minMN; i++) {
        fprintf( stderr, "%e ", s[i] );
    }
    fprintf(stderr, "\n");

    free(W);
}

/**
 *******************************************************************************
 *
 * @brief Check the orthogonality of the matrix A
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the matrix A.
 *
 * @param[in] N
 *          The number of columns of the matrix A.
 *
 * @param[in] A
 *          The matrix A to study of size lda-by-N
 *
 * @param[in] lda
 *          The leading dimension of the matrix A. lda = max( 1, M )
 *
 *******************************************************************************
 *
 * @retval 0 if the matrix A is orthogonal
 * @retval 1 if the matrix A is not orthogonal
 *
 *******************************************************************************/
int
core_zlrdbg_check_orthogonality( pastix_int_t              M,
                                 pastix_int_t              N,
                                 const pastix_complex64_t *A,
                                 pastix_int_t              lda )
{
    pastix_complex64_t *Id;
    double alpha, beta;
    double normQ, res;
    pastix_int_t info_ortho;
    pastix_int_t minMN = pastix_imin(M, N);
    pastix_int_t maxMN = pastix_imax(M, N);
    double eps = LAPACKE_dlamch_work('e');

    alpha = 1.0;
    beta  = -1.0;

    /* Build the identity matrix */
    Id = malloc( minMN * minMN * sizeof(pastix_complex64_t) );
    LAPACKE_zlaset_work( LAPACK_COL_MAJOR, 'A', minMN, minMN,
                         0., 1., Id, minMN );

    if (M > N) {
        /* Perform Id - Q'Q */
        cblas_zherk(CblasColMajor, CblasUpper, CblasConjTrans, N, M, alpha, A, lda, beta, Id, minMN);
    }
    else {
        /* Perform Id - QQ' */
        cblas_zherk(CblasColMajor, CblasUpper, CblasNoTrans,   M, N, alpha, A, lda, beta, Id, minMN);
    }

    normQ = LAPACKE_zlanhe_work( LAPACK_COL_MAJOR, 'M', 'U', minMN, Id, minMN, NULL );
    res = normQ / (maxMN * eps);

    if ( isnan(res) || isinf(res) || (res > 60.0) ) {
        fprintf(stderr, "Check Orthogonality: || I - Q*Q' || = %e, ||Id-Q'*Q||_oo / (N*eps) = %e : \n",
                normQ, res );
        info_ortho = 1;
    }
    else {
        info_ortho = 0;
    }

    free(Id);
    return info_ortho;
}

/**
 *******************************************************************************
 *
 * @brief Check the orthogonality of the matrix A relatively to the matrix B
 *
 * Check that A^t B = 0
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the matrix A.
 *
 * @param[in] NA
 *          The number of columns of the matrix A.
 *
 * @param[in] NB
 *          The number of columns of the matrix B.
 *
 * @param[in] A
 *          The matrix A to study of size lda-by-NA
 *
 * @param[in] lda
 *          The leading dimension of the matrix A. lda = max( 1, M )
 *
 * @param[in] B
 *          The matrix B to study of size ldb-by-NB
 *
 * @param[in] ldb
 *          The leading dimension of the matrix B. ldb = max( 1, M )
 *
 *******************************************************************************
 *
 * @retval 0 if the matrices A and B are orthogonal
 * @retval 1 if the matrices A anb B are not orthogonal
 *
 *******************************************************************************/
int
core_zlrdbg_check_orthogonality_AB( pastix_int_t M, pastix_int_t NA, pastix_int_t NB,
                                    const pastix_complex64_t *A, pastix_int_t lda,
                                    const pastix_complex64_t *B, pastix_int_t ldb )
{
    pastix_complex64_t *Zero;
    double norm, res;
    pastix_int_t info_ortho;
    double eps = LAPACKE_dlamch_work('e');
    pastix_complex64_t zone = 1.0;
    pastix_complex64_t zzero = 0.0;

    /* Build the null matrix */
    Zero = malloc( NA * NB * sizeof(pastix_complex64_t) );
    LAPACKE_zlaset_work( LAPACK_COL_MAJOR, 'A', NA, NB,
                         0., 0., Zero, NA );

    cblas_zgemm(CblasColMajor, CblasConjTrans, CblasNoTrans,
                NA, NB, M,
                CBLAS_SADDR(zone),  A, lda,
                                    B, ldb,
                CBLAS_SADDR(zzero), Zero, NA);

    norm = LAPACKE_zlange_work( LAPACK_COL_MAJOR, 'M', NA, NB, Zero, NA, NULL );
    res = norm / (M * eps);

    if ( isnan(res) || isinf(res) || (res > 60.0) ) {
        fprintf(stderr, "Check Orthogonality: || A' B || = %e, || A' B ||_oo / (M*eps) = %e : \n",
                norm, res );
        info_ortho = 1;
    }
    else {
        info_ortho = 0;
    }

    free(Zero);
    return info_ortho;
}
