/**
 *
 * @file core_zgeadd.c
 *
 * PaStiX kernel routines
 *
 * @copyright 2010-2015 Univ. of Tennessee, Univ. of California Berkeley and
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
 * @brief Add two matrices together.
 *
 * Perform the operation:  B <- alpha * op(A) + B
 *
 *******************************************************************************
 *
 * @param[in] trans
 *         @arg PastixNoTrans:   No transpose, op( A ) = A;
 *         @arg PastixTrans:     Transpose, op( A ) = A';
 *         @arg PastixConjTrans: Conjugate Transpose, op( A ) = conj(A').
 *
 * @param[in] M
 *          Number of rows of the matrix B.
 *          Number of rows of the matrix A, if trans == PastixNoTrans, number of
 *          columns of A otherwise.
 *
 * @param[in] N
 *          Number of columns of the matrix B.
 *          Number of columns of the matrix A, if trans == PastixNoTrans, number
 *          of rows of A otherwise.
 *
 * @param[in] alpha
 *          Scalar factor of A.
 *
 * @param[in] A
 *          Matrix of size LDA-by-N, if trans == PastixNoTrans, LDA-by-M,
 *          otherwise.
 *
 * @param[in] LDA
 *          Leading dimension of the array A. LDA >= max(1,K).
 *          K = M if trans == PastixNoTrans, K = N otherwise.
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
 * @retval PASTIX_SUCCESS successful exit
 * @retval <0 if -i, the i-th argument had an illegal value
 * @retval 1, not yet implemented
 *
 ******************************************************************************/
int
core_zgeadd( pastix_trans_t            trans,
             pastix_int_t              M,
             pastix_int_t              N,
             pastix_complex64_t        alpha,
             const pastix_complex64_t *A,
             pastix_int_t              LDA,
             pastix_complex64_t        beta,
             pastix_complex64_t       *B,
             pastix_int_t              LDB)
{
    int i, j;

#if !defined(NDEBUG)
    if ((trans < PastixNoTrans)   ||
        (trans > PastixConjTrans))
    {
        return -1;
    }

    if (M < 0) {
        return -2;
    }
    if (N < 0) {
        return -3;
    }
    if ( ((trans == PastixNoTrans) && (LDA < pastix_imax(1,M)) && (M > 0)) ||
         ((trans != PastixNoTrans) && (LDA < pastix_imax(1,N)) && (N > 0)) )
    {
        return -6;
    }
    if ( (LDB < pastix_imax(1,M)) && (M > 0) ) {
        return -8;
    }
#endif

    switch( trans ) {
#if defined(PRECISION_z) || defined(PRECISION_c)
    case PastixConjTrans:
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

    case PastixTrans:
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

    case PastixNoTrans:
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
    return PASTIX_SUCCESS;
}

/**
 ******************************************************************************
 *
 * @brief Add two column blocks together.
 *
 * Perform the operation: cblk2 <- cblk1 + cblk2 using core_zgeadd().
 *
 *******************************************************************************
 *
 * @param[in] cblk1
 *          The pointer to the data structure that describes the panel to add.
 *          Next column blok must be accessible through cblk1[1].
 *
 * @param[in] cblk2
 *          The pointer to the data structure that describes the panel in which
 *          we add.
 *          Next column blok must be accessible through cblk2[1].
 *
 * @param[in] L
 *          The pointer to the lower matrix storing the coefficients of the
 *          panel. Must be of size cblk.stride -by- cblk.width
 *
 * @param[inout] Cl
 *          The pointer to the lower matrix storing the coefficients of the
 *          updated panel. Must be of size cblk.stride -by- cblk.width
 *          On exit Cl = Cl + L.
 *
 * @param[in] U
 *          The pointer to the upper matrix storing the coefficients of the
 *          panel. Must be of size cblk.stride -by- cblk.width. Ignored if
 *          NULL.
 *
 * @param[inout] Cu
 *          The pointer to the upper matrix storing the coefficients of the
 *          updated panel. Must be of size cblk.stride -by- cblk.width
 *          On exit Cu = Cu + U.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS on success.
 *
 ******************************************************************************/
int
cpucblk_zgeaddsp1d( const SolverCblk         *cblk1,
                    SolverCblk               *cblk2,
                    const pastix_complex64_t *L,
                    pastix_complex64_t       *Cl,
                    const pastix_complex64_t *U,
                    pastix_complex64_t       *Cu )
{
    SolverBlok *iterblok;
    SolverBlok *firstblok;
    SolverBlok *lastblok;
    SolverBlok *fblok;
    pastix_int_t ncol1 = cblk1->lcolnum - cblk1->fcolnum + 1;
    const pastix_complex64_t *ga;
    pastix_complex64_t *gb;

    assert(0 /* Outdated */);

    firstblok = cblk1->fblokptr;
    lastblok  = cblk1[1].fblokptr;
    fblok = cblk2->fblokptr;

    assert(cblk1->fcolnum >= cblk2->fcolnum);

    for (iterblok = firstblok; iterblok < lastblok; iterblok++) {
        pastix_int_t nrow;
        /* Find facing bloknum */
        while (!is_block_inside_fblock( iterblok, fblok )) {
            fblok++;
            assert( fblok < cblk2[1].fblokptr );
        }
        ga = L + iterblok->coefind;
        gb = Cl + cblk2->stride*(cblk1->fcolnum-cblk2->fcolnum) +
            fblok->coefind +
            iterblok->frownum - fblok->frownum;
        nrow = iterblok->lrownum - iterblok->frownum + 1;
        core_zgeadd( PastixNoTrans,
                     nrow, ncol1,
                     1.0, ga, cblk1->stride,
		     1.0, gb, cblk2->stride );
        if (U != NULL) {
            ga = U + iterblok->coefind;
            gb = Cu + cblk2->stride*(cblk1->fcolnum-cblk2->fcolnum) +
                fblok->coefind +
                iterblok->frownum - fblok->frownum;
            core_zgeadd( PastixNoTrans,
                         nrow, ncol1,
                         1.0, ga, cblk1->stride,
                         1.0, gb, cblk2->stride );
        }
    }
    return PASTIX_SUCCESS;
}
