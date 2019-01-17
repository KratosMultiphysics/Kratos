/**
 *
 * @file core_zgemdm.c
 *
 * PaStiX kernel routines.
 *
 * @copyright 2010-2015 Univ. of Tennessee, Univ. of California Berkeley and
 *                      Univ. of Colorado Denver. All rights reserved.
 * @copyright 2015-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.1
 * @author Dulceneia Becker
 * @author Mathieu Faverge
 * @date 2018-07-16
 * @precisions normal z -> c d s
 *
 **/
#include "common.h"
#include "cblas.h"
#include "lapacke.h"

/**
 ******************************************************************************
 *
 * @brief Perform one of the following matrix-matrix operations
 *
 *       C := alpha*op( A )*D*op( B ) + beta*C,
 *
 * where op( X ) is one of
 *
 *       op( X ) = X   or   op( X ) = X',
 *
 * alpha and beta are scalars, and A, B, C and D are matrices, with
 *
 *       op( A ) an m by k matrix,
 *       op( B ) an k by n matrix,
 *       C an m by n matrix and
 *       D an k by k matrix.
 *
 *******************************************************************************
 *
 * @param[in] transA
 *         @arg PastixNoTrans   : No transpose, op( A ) = A;
 *         @arg PastixConjTrans : Transpose, op( A ) = A'.
 *
 * @param[in] transB
 *         @arg PastixNoTrans   : No transpose, op( B ) = B;
 *         @arg PastixConjTrans : Transpose, op( B ) = B'.
 *
 * @param[in] M
 *         The number of rows of the matrix op( A ) and of the
 *         matrix C. M must be at least zero.
 *
 * @param[in] N
 *         The number of columns of the matrix op( B ) and the
 *         number of columns of the matrix C. N must be at least zero.
 *
 * @param[in] K
 *         The number of columns of the matrix op( A ), the number of
 *         rows of the matrix op( B ), and the number of rows and columns
 *         of matrix D. K must be at least zero.
 *
 * @param[in] alpha
 *         On entry, alpha specifies the scalar alpha.
 *
 * @param[in] A
 *         pastix_complex64_t array of DIMENSION ( LDA, ka ), where ka is
 *         k when transA = PastixTrans, and is m otherwise.
 *         Before entry with transA = PastixTrans, the leading m by k
 *         part of the array A must contain the matrix A, otherwise
 *         the leading k by m part of the array A must contain the
 *         matrix A.
 *
 * @param[in] LDA
 *        On entry, LDA specifies the first dimension of A as declared
 *        in the calling (sub) program. When transA = PastixTrans then
 *        LDA must be at least  max( 1, m ), otherwise LDA must be at
 *        least max( 1, k ).
 *
 * @param[in] B
 *        pastix_complex64_t array of DIMENSION ( LDB, kb ), where kb is
 *        n  when transB = PastixTrans, and is k otherwise.
 *        Before entry with transB = PastixTrans, the leading k by n
 *        part of the array B must contain the matrix B, otherwise
 *        the leading n by k part of the array B must contain the
 *        matrix B.
 *
 * @param[in] LDB
 *       On entry, LDB specifies the first dimension of B as declared
 *       in the calling (sub) program. When transB = PastixTrans then
 *       LDB must be at least  max( 1, k ), otherwise LDB must be at
 *       least max( 1, n ).
 *
 * @param[in] beta
 *       On entry, beta specifies the scalar beta. When beta is
 *       supplied as zero then C need not be set on input.
 *
 * @param[in] C
 *       pastix_complex64_t array of DIMENSION ( LDC, n ).
 *       Before entry, the leading m by n part of the array C must
 *       contain the matrix C,  except when beta is zero, in which
 *       case C need not be set on entry.
 *       On exit, the array C is overwritten by the m by n matrix
 *       ( alpha*op( A )*D*op( B ) + beta*C ).
 *
 * @param[in] LDC
 *       On entry, LDC specifies the first dimension of C as declared
 *       in the calling (sub) program. LDC must be at least
 *       max( 1, m ).
 *
 * @param[in] D
 *        pastix_complex64_t array of DIMENSION ( LDD, k ).
 *        Before entry, the leading k by k part of the array D
 *        must contain the matrix D.
 *
 * @param[in] incD
 *       On entry, LDD specifies the first dimension of D as declared
 *       in the calling (sub) program. LDD must be at least
 *       max( 1, k ).
 *
 * @param[in] WORK
 *       pastix_complex64_t array, dimension (MAX(1,LWORK))
 *
 * @param[in] LWORK
 *       The length of WORK.
 *       On entry, if transA = PastixTrans and transB = PastixTrans then
 *       LWORK >= max(1, K*N). Otherwise LWORK >= max(1, M*K).
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS successful exit
 * @retval <0 if -i, the i-th argument had an illegal value
 *
 ******************************************************************************/
int
core_zgemdm( pastix_trans_t transA, pastix_trans_t transB,
             int M, int N, int K,
                   pastix_complex64_t alpha,
             const pastix_complex64_t *A,    int LDA,
             const pastix_complex64_t *B,    int LDB,
                   pastix_complex64_t beta,
                   pastix_complex64_t *C,    int LDC,
             const pastix_complex64_t *D,    int incD,
                   pastix_complex64_t *WORK, int LWORK )
{
    int j, Am, Bm;
    pastix_complex64_t delta;
    pastix_complex64_t *wD2, *w;
    const pastix_complex64_t *wD;

    Am = (transA == PastixNoTrans ) ? M : K;
    Bm = (transB == PastixNoTrans ) ? K : N;

    /* Check input arguments */
    if ((transA < PastixNoTrans) || (transA > PastixConjTrans)) {
        return -1;
    }
    if ((transB < PastixNoTrans) || (transB > PastixConjTrans)) {
        return -2;
    }
    if (M < 0) {
        return -3;
    }
    if (N < 0) {
        return -4;
    }
    if (K < 0) {
        return -5;
    }
    if ((LDA < pastix_imax(1,Am)) && (Am > 0)) {
        return -8;
    }
    if ((LDB < pastix_imax(1,Bm)) && (Bm > 0)) {
        return -10;
    }
    if ((LDC < pastix_imax(1,M)) && (M > 0)) {
        return -13;
    }
    if ( incD < 0 ) {
        return -15;
    }
    if ( ( ( transA == PastixNoTrans ) && ( LWORK < (M+1)*K) ) ||
         ( ( transA != PastixNoTrans ) && ( LWORK < (N+1)*K) ) ){
        errorPrint("CORE_gemdm: Illegal value of LWORK\n");
        if (transA == PastixNoTrans ) {
            errorPrint("LWORK %d < (M=%d+1)*K=%d ", LWORK, M, K);
        }
        if (transA == PastixNoTrans ) {
            errorPrint("LWORK %d < (N=%d+1)*K=%d ", LWORK, N, K);
        }
        return -17;
    }

    /* Quick return */
    if (M == 0 || N == 0 ||
        ((alpha == 0.0 || K == 0) && beta == 1.0) ) {
        return PASTIX_SUCCESS;
    }

    if ( incD == 1 ) {
        wD = D;
    } else {
        wD2 = WORK;
        cblas_zcopy(K, D, incD, wD2, 1);
        wD = wD2;
    }
    w = WORK + K;

    /*
     * transA == PastixNoTrans
     */
    if ( transA == PastixNoTrans )
    {
        /* WORK = A * D */
        for (j=0; j<K; j++, wD++) {
            delta = *wD;
            cblas_zcopy(M, &A[LDA*j], 1,       &w[M*j], 1);
            cblas_zscal(M, CBLAS_SADDR(delta), &w[M*j], 1);
        }

        /* C = alpha * WORK * op(B) + beta * C */
        cblas_zgemm(CblasColMajor, CblasNoTrans, (CBLAS_TRANSPOSE)transB,
                    M, N, K,
                    CBLAS_SADDR(alpha), w, M,
                                        B, LDB,
                    CBLAS_SADDR(beta),  C, LDC);
    }
    else
    {
        if ( transB == PastixNoTrans ) /* Worst case*/
        {
            /* WORK = (D * B)' */
            for (j=0; j<K; j++, wD++) {
                delta = *wD;
                cblas_zcopy(N, &B[j],     LDB,     &w[N*j], 1);
                cblas_zscal(N, CBLAS_SADDR(delta), &w[N*j], 1);
            }

            /* C = alpha * op(A) * WORK' + beta * C */
            cblas_zgemm(CblasColMajor, (CBLAS_TRANSPOSE)transA, CblasTrans,
                        M, N, K,
                        CBLAS_SADDR(alpha), A, LDA,
                                            w, N,
                        CBLAS_SADDR(beta),  C, LDC);
        }
        else
        {
#if defined(PRECISION_z) || defined(PRECISION_c)
            if ( transB == PastixConjTrans )
            {
                /* WORK = D * B' */
                for (j=0; j<K; j++, wD++) {
                    delta = *wD;
                    cblas_zcopy(N, &B[LDB*j], 1,       &w[N*j], 1);
                    LAPACKE_zlacgv_work(N,             &w[N*j], 1);
                    cblas_zscal(N, CBLAS_SADDR(delta), &w[N*j], 1);
                }
            }
            else
#endif
            {
                /* WORK = D * B' */
              for (j=0; j<K; j++, wD++) {
                    delta = *wD;
                    cblas_zcopy(N, &B[LDB*j], 1,       &w[N*j], 1);
                    cblas_zscal(N, CBLAS_SADDR(delta), &w[N*j], 1);
                }
            }

            /* C = alpha * op(A) * WORK + beta * C */
            cblas_zgemm(CblasColMajor, (CBLAS_TRANSPOSE)transA, CblasNoTrans,
                        M, N, K,
                        CBLAS_SADDR(alpha), A, LDA,
                                            w, N,
                        CBLAS_SADDR(beta),  C, LDC);
        }
    }
    return PASTIX_SUCCESS;
}
