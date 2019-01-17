/**
 *
 * @file core_ztradd.c
 *
 * @copyright 2012-2014 Univ. of Tennessee, Univ. of California Berkeley and
 *                      Univ. of Colorado Denver. All rights reserved.
 * @copyright 2012-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.1
 * @author Mathieu Faverge
 * @date 2018-07-16
 * @precisions normal z -> c d s
 *
 **/
#include "common.h"
#include "blend/solver.h"
#include "pastix_zcores.h"
#include "cblas.h"

/**
 ******************************************************************************
 *
 * @brief Add two triangular matrices together as in PBLAS pztradd.
 *
 *    B <- alpha * op(A)  + beta * B,
 *
 * where op(X) = X, X', or conj(X')
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          Specifies the shape of A and B matrices:
 *          @arg PastixUpperLower: A and B are general matrices.
 *          @arg PastixUpper: op(A) and B are upper trapezoidal matrices.
 *          @arg PastixLower: op(A) and B are lower trapezoidal matrices.
 *
 * @param[in] trans
 *          Specifies whether the matrix A is non-transposed, transposed, or
 *          conjugate transposed
 *          @arg PastixNoTrans:   op(A) = A
 *          @arg PastixTrans:     op(A) = A'
 *          @arg PastixConjTrans: op(A) = conj(A')
 *
 * @param[in] M
 *          Number of rows of the matrices op(A) and B.
 *
 * @param[in] N
 *          Number of columns of the matrices op(A) and B.
 *
 * @param[in] alpha
 *          Scalar factor of A.
 *
 * @param[in] A
 *          Matrix of size LDA-by-N, if trans = PastixNoTrans, LDA-by-M
 *          otherwise.
 *
 * @param[in] LDA
 *          Leading dimension of the array A. LDA >= max(1,M) if trans =
 *          PastixNoTrans, LDA >= max(1,N) otherwise.
 *
 * @param[in] beta
 *          Scalar factor of B.
 *
 * @param[inout] B
 *          Matrix of size LDB-by-N.
 *          On exit, B = alpha * op(A) + beta * B
 *
 * @param[in] LDB
 *          Leading dimension of the array B. LDB >= max(1,M)
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS successful exit
 * @retval <0 if -i, the i-th argument had an illegal value
 *
 ******************************************************************************/
int
core_ztradd( pastix_uplo_t uplo, pastix_trans_t trans, pastix_int_t M, pastix_int_t N,
             pastix_complex64_t alpha, const pastix_complex64_t *A, pastix_int_t LDA,
             pastix_complex64_t beta,        pastix_complex64_t *B, pastix_int_t LDB)
{
    int i, j;

    if (uplo == PastixUpperLower){
        int rc = core_zgeadd( trans, M, N, alpha, A, LDA, beta, B, LDB );
        if (rc != PASTIX_SUCCESS)
            return rc-1;
        else
            return rc;
    }

#if !defined(NDEBUG)
    if ((uplo != PastixUpper) &&
        (uplo != PastixLower))
    {
        return -1;
    }

    if ((trans < PastixNoTrans)   ||
        (trans > PastixConjTrans))
    {
        return -2;
    }

    if (M < 0) {
        return -3;
    }
    if (N < 0) {
        return -4;
    }
    if ( ((trans == PastixNoTrans) && (LDA < pastix_imax(1,M)) && (M > 0)) ||
         ((trans != PastixNoTrans) && (LDA < pastix_imax(1,N)) && (N > 0)) )
    {
        return -7;
    }
    if ( (LDB < pastix_imax(1,M)) && (M > 0) ) {
        return -9;
    }
#endif

    /**
     * PastixLower
     */
    if (uplo == PastixLower) {
        switch( trans ) {
#if defined(PRECISION_z) || defined(PRECISION_c)
        case PastixConjTrans:
            for (j=0; j<N; j++, A++) {
                for(i=j; i<M; i++, B++) {
                    *B = beta * (*B) + alpha * conj(A[LDA*i]);
                }
                B += LDB-M+j+1;
            }
            break;
#endif /* defined(PRECISION_z) || defined(PRECISION_c) */

        case PastixTrans:
            for (j=0; j<N; j++, A++) {
                for(i=j; i<M; i++, B++) {
                    *B = beta * (*B) + alpha * A[LDA*i];
                }
                B += LDB-M+j+1;
            }
            break;

        case PastixNoTrans:
        default:
            for (j=0; j<N; j++) {
                for(i=j; i<M; i++, B++, A++) {
                    *B = beta * (*B) + alpha * (*A);
                }
                B += LDB-M+j+1;
                A += LDA-M+j+1;
            }
        }
    }
    /**
     * PastixUpper
     */
    else {
        switch( trans ) {
#if defined(PRECISION_z) || defined(PRECISION_c)
        case PastixConjTrans:
            for (j=0; j<N; j++, A++) {
                int mm = pastix_imin( j+1, M );
                for(i=0; i<mm; i++, B++) {
                    *B = beta * (*B) + alpha * conj(A[LDA*i]);
                }
                B += LDB-mm;
            }
            break;
#endif /* defined(PRECISION_z) || defined(PRECISION_c) */

        case PastixTrans:
            for (j=0; j<N; j++, A++) {
                int mm = pastix_imin( j+1, M );
                for(i=0; i<mm; i++, B++) {
                    *B = beta * (*B) + alpha * (A[LDA*i]);
                }
                B += LDB-mm;
            }
            break;

        case PastixNoTrans:
        default:
            for (j=0; j<N; j++) {
                int mm = pastix_imin( j+1, M );
                for(i=0; i<mm; i++, B++, A++) {
                    *B = beta * (*B) + alpha * (*A);
                }
                B += LDB-mm;
                A += LDA-mm;
            }
        }
    }
    return PASTIX_SUCCESS;
}
