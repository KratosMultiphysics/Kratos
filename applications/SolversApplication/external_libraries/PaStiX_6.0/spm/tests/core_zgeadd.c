/**
 *
 * @file core_zgeadd.c
 *
 * PaStiX kernel routines
 *
 * @copyright 2010-2015 Univ. of Tennessee, Univ. of California Berkeley and
 *                      Univ. of Colorado Denver. All rights reserved.
 * @copyright 2012-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.0
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 **/
#include <spm_tests.h>
#include "cblas.h"

/**
 ******************************************************************************
 *
 * @brief Add two matrices together.
 *
 * Perform the operation:  B <- alpha * op(A) + B
 *
 *******************************************************************************
 *
 * @param[in] trans
 *         @arg SpmNoTrans:   No transpose, op( A ) = A;
 *         @arg SpmTrans:     Transpose, op( A ) = A';
 *         @arg SpmConjTrans: Conjugate Transpose, op( A ) = conj(A').
 *
 * @param[in] M
 *          Number of rows of the matrix B.
 *          Number of rows of the matrix A, if trans == SpmNoTrans, number of
 *          columns of A otherwise.
 *
 * @param[in] N
 *          Number of columns of the matrix B.
 *          Number of columns of the matrix A, if trans == SpmNoTrans, number
 *          of rows of A otherwise.
 *
 * @param[in] alpha
 *          Scalar factor of A.
 *
 * @param[in] A
 *          Matrix of size LDA-by-N, if trans == SpmNoTrans, LDA-by-M,
 *          otherwise.
 *
 * @param[in] LDA
 *          Leading dimension of the array A. LDA >= max(1,K).
 *          K = M if trans == SpmNoTrans, K = N otherwise.
 *
 * @param[in] beta
 *          Scalar factor of B.
 *
 * @param[inout] B
 *          Matrix of size LDB-by-N.
 *
 * @param[in] LDB
 *          Leading dimension of the array B. LDB >= max(1,M)
 *
 *******************************************************************************
 *
 * @retval SPM_SUCCESS successful exit
 * @retval <0 if -i, the i-th argument had an illegal value
 * @retval 1, not yet implemented
 *
 ******************************************************************************/
int
core_zgeadd( spm_trans_t            trans,
             spm_int_t              M,
             spm_int_t              N,
             spm_complex64_t        alpha,
             const spm_complex64_t *A,
             spm_int_t              LDA,
             spm_complex64_t        beta,
             spm_complex64_t       *B,
             spm_int_t              LDB)
{
    int i, j;

#if !defined(NDEBUG)
    if ((trans < SpmNoTrans)   ||
        (trans > SpmConjTrans))
    {
        return -1;
    }

    if (M < 0) {
        return -2;
    }
    if (N < 0) {
        return -3;
    }
    if ( ((trans == SpmNoTrans) && (LDA < spm_imax(1,M)) && (M > 0)) ||
         ((trans != SpmNoTrans) && (LDA < spm_imax(1,N)) && (N > 0)) )
    {
        return -6;
    }
    if ( (LDB < spm_imax(1,M)) && (M > 0) ) {
        return -8;
    }
#endif

    switch( trans ) {
#if defined(PRECISION_z) || defined(PRECISION_c)
    case SpmConjTrans:
        if ( alpha == 0.0 ) {
            for (j=0; j<N; j++) {
                for(i=0; i<M; i++, B++) {
                    *B = beta * (*B);
                }
                B += LDB-M;
            }
        }
        else if ( beta == 0.0 ) {
            for (j=0; j<N; j++, A++) {
                for(i=0; i<M; i++, B++) {
                    *B = alpha * conj(A[LDA*i]);
                }
                B += LDB-M;
            }
        }
        else {
            for (j=0; j<N; j++, A++) {
                for(i=0; i<M; i++, B++) {
                    *B = beta * (*B) + alpha * conj(A[LDA*i]);
                }
                B += LDB-M;
            }
        }
        break;
#endif /* defined(PRECISION_z) || defined(PRECISION_c) */

    case SpmTrans:
        if ( alpha == 0.0 ) {
            for (j=0; j<N; j++) {
                for(i=0; i<M; i++, B++) {
                    *B = beta * (*B);
                }
                B += LDB-M;
            }
        }
        else if ( beta == 0.0 ) {
            for (j=0; j<N; j++, A++) {
                for(i=0; i<M; i++, B++) {
                    *B = alpha * A[LDA*i];
                }
                B += LDB-M;
            }
        }
        else {
            for (j=0; j<N; j++, A++) {
                for(i=0; i<M; i++, B++) {
                    *B = beta * (*B) + alpha * A[LDA*i];
                }
                B += LDB-M;
            }
        }
        break;

    case SpmNoTrans:
    default:
        if ( alpha == 0.0 ) {
            for (j=0; j<N; j++) {
                for(i=0; i<M; i++, B++) {
                    *B = beta * (*B);
                }
                B += LDB-M;
            }
        }
        else if ( beta == 0.0 ) {
            for (j=0; j<N; j++) {
                for(i=0; i<M; i++, B++, A++) {
                    *B = alpha * (*A);
                }
                A += LDA-M;
                B += LDB-M;
            }
        }
        else {
            for (j=0; j<N; j++) {
                for(i=0; i<M; i++, B++, A++) {
                    *B = beta * (*B) + alpha * (*A);
                }
                A += LDA-M;
                B += LDB-M;
            }
        }
    }
    return SPM_SUCCESS;
}
