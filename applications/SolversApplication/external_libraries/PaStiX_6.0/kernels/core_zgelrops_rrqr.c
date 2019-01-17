/**
 *
 * @file core_zgelrops_rrqr.c
 *
 * PaStiX low-rank kernel routines using Rank-revealing QR based on Lapack GEQP3.
 *
 * @copyright 2016-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.1
 * @author Alfredo Buttari
 * @author Gregoire Pichon
 * @date 2018-07-16
 * @precisions normal z -> c d s
 *
 **/
#include "common.h"
#include <cblas.h>
#include <lapacke.h>
#include "flops.h"
#include "kernels_trace.h"
#include "blend/solver.h"
#include "pastix_zcores.h"
#include "pastix_zlrcores.h"
#include "z_nan_check.h"
#include "pastix_lowrank.h"

#ifndef DOXYGEN_SHOULD_SKIP_THIS
static pastix_complex64_t mzone = -1.0;
static pastix_complex64_t zone  =  1.0;
static pastix_complex64_t zzero =  0.0;
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/**
 *******************************************************************************
 *
 * @brief Try to orthognalize the u part of the low-rank form, and update the v
 * part accordingly using full QR.
 *
 * This function considers a low-rank matrix resulting from the addition of two
 * matrices B += A, with A of smaller or equal size to B.
 * The product has the form: U * V^t
 *
 * The U part of the low-rank form must be orthognalized to get the smaller
 * possible rank during the rradd operation. This function perfoms this by
 * applying a full QR factorization on the U part.
 *
 *  U = Q R, then U' = Q, and V' = R * V
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the u1u2 matrix.
 *
 * @param[in] N
 *          The number of columns of the v1v2 matrix.
 *
 * @param[in] rank
 *          The number of columns of the U matrix, and the number of rows of the
 *          V part in the v1v2 matrix.
 *
 * @param[inout] U
 *          The U matrix of size ldu -by- rank. On exit, Q from U = Q R.
 *
 * @param[in] ldu
 *          The leading dimension of the U matrix. ldu >= max(1, M)
 *
 * @param[inout] V
 *          The V matrix of size ldv -by- N.
 *          On exit, R * V, with R from U = Q R.
 *
 * @param[in] ldv
 *          The leading dimension of the V matrix. ldv >= max(1, rank)
 *
 *******************************************************************************
 *
 * @return The number of flops required to perform the operation.
 *
 *******************************************************************************/
static inline pastix_fixdbl_t
core_zlrorthu_fullqr( pastix_int_t M,  pastix_int_t N, pastix_int_t rank,
                      pastix_complex64_t *U, pastix_int_t ldu,
                      pastix_complex64_t *V, pastix_int_t ldv )
{
    pastix_int_t minMK = pastix_imin( M, rank );
    pastix_int_t lwork = M * 32 + minMK;
    pastix_int_t ret;
    pastix_complex64_t *W = malloc( lwork * sizeof(pastix_complex64_t) );
    pastix_complex64_t *tau, *work;
    pastix_fixdbl_t flops = 0.;

    tau  = W;
    work = W + minMK;
    lwork -= minMK;

    assert( M >= rank );

    /* Compute U = Q * R */
    ret = LAPACKE_zgeqrf_work( LAPACK_COL_MAJOR, M, rank,
                               U, ldu, tau, work, lwork );
    assert( ret == 0 );
    flops += FLOPS_ZGEQRF( M, rank );

    /* Compute V' = R * V' */
    cblas_ztrmm( CblasColMajor,
                 CblasLeft, CblasUpper,
                 CblasNoTrans, CblasNonUnit,
                 rank, N, CBLAS_SADDR(zone),
                 U, ldu, V, ldv );
    flops += FLOPS_ZTRMM( PastixLeft, rank, N );

    /* Generate the Q */
    ret = LAPACKE_zungqr_work( LAPACK_COL_MAJOR, M, rank, rank,
                               U, ldu, tau, work, lwork );
    assert( ret == 0 );
    flops += FLOPS_ZUNGQR( M, rank, rank );

    free(W);

    (void)ret;
    return flops;
}

/**
 *******************************************************************************
 *
 * @brief Try to orthognalize the U part of the low-rank form, and update the V
 * part accordingly using partial QR.
 *
 * This function considers a low-rank matrix resulting from the addition of two
 * matrices B += A, with A of smaller or equal size to B.
 * The product has the form: U * V^t
 *
 * The U part of the low-rank form must be orthognalized to get the smaller
 * possible rank during the rradd operation. This function perfoms this by
 * applying a full QR factorization on the U part.
 *
 * In that case, it takes benefit from the fact that U = [ u1, u2 ], and V = [
 * v1, v2 ] with u2 and v2 wich are matrices of respective size M2-by-r2, and
 * r2-by-N2, offset by offx and offy
 *
 * The steps are:
 *    - Scaling of u2 with removal of the null columns
 *    - Orthogonalization of u2 relatively to u1
 *    - Application of the update to v2
 *    - orthogonalization through QR of u2
 *    - Update of V
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the u1u2 matrix.
 *
 * @param[in] N
 *          The number of columns of the v1v2 matrix.
 *
 * @param[in] r1
 *          The number of columns of the U matrix in the u1 part, and the number
 *          of rows of the V part in the v1 part.
 *
 * @param[inout] r2ptr
 *          The number of columns of the U matrix in the u2 part, and the number
 *          of rows of the V part in the v2 part. On exit, this rank is reduced
 *          y the number of null columns found in U.
 *
 * @param[in] offx
 *          The row offset of the matrix u2 in U.
 *
 * @param[in] offy
 *          The column offset of the matrix v2 in V.
 *
 * @param[inout] U
 *          The U matrix of size ldu -by- rank. On exit, the orthogonalized U.
 *
 * @param[in] ldu
 *          The leading dimension of the U matrix. ldu >= max(1, M)
 *
 * @param[inout] V
 *          The V matrix of size ldv -by- N.
 *          On exit, the updated V matrix.
 *
 * @param[in] ldv
 *          The leading dimension of the V matrix. ldv >= max(1, rank)
 *
 *******************************************************************************
 *
 * @return The number of flops required to perform the operation.
 *
 *******************************************************************************/
static inline pastix_fixdbl_t
core_zlrorthu_partialqr( pastix_int_t M,  pastix_int_t N,
                         pastix_int_t r1, pastix_int_t *r2ptr,
                         pastix_int_t offx, pastix_int_t offy,
                         pastix_complex64_t *U, pastix_int_t ldu,
                         pastix_complex64_t *V, pastix_int_t ldv )
{
    pastix_int_t r2 = *r2ptr;
    pastix_int_t minMN = pastix_imin( M, r2 );
    pastix_int_t ldwork = pastix_imax( r1 * r2, M * 32 + minMN );
    pastix_int_t ret, i;
    pastix_complex64_t *u1 = U;
    pastix_complex64_t *u2 = U + r1 * ldu;
    pastix_complex64_t *v1 = V;
    pastix_complex64_t *v2 = V + r1;
    pastix_complex64_t *W = malloc( ldwork * sizeof(pastix_complex64_t) );
    pastix_complex64_t *tau, *work;
    pastix_fixdbl_t flops = 0.;
    double norm, eps;

    tau = W;
    work = W + minMN;
    ldwork -= minMN;

    eps = LAPACKE_dlamch_work('e');

    /* Scaling */
    for (i=0; i<r2; i++, u2 += ldu, v2++) {
        norm = cblas_dznrm2( M, u2, 1 );
        if ( norm > (M * eps) ) {
            cblas_zdscal( M, 1. / norm, u2, 1   );
            cblas_zdscal( N, norm,      v2, ldv );
        }
        else {
            if ( i < (r2-1) ) {
                cblas_zswap( M, u2, 1, U + (r1+r2-1) * ldu, 1 );
                memset( U + (r1+r2-1) * ldu, 0,  M * sizeof(pastix_complex64_t) );

                cblas_zswap( N, v2, ldv, V + (r1+r2-1),     ldv );
                LAPACKE_zlaset_work( LAPACK_COL_MAJOR, 'A', 1, N,
                                     0., 0., V + (r1+r2-1), ldv );
                r2--;
                i--;
                u2-= ldu;
                v2--;
            }
            else {
                memset( u2, 0,  M * sizeof(pastix_complex64_t) );
                LAPACKE_zlaset_work( LAPACK_COL_MAJOR, 'A', 1, N,
                                     0., 0., v2, ldv );
                r2--;
            }
        }
    }
    u2 = U + r1 * ldu;
    v2 = V + r1;

    *r2ptr = r2;

    if ( r2 == 0 ) {
        free( W );
        return 0.;
    }

    /* Compute W = u1^t u2 */
    cblas_zgemm( CblasColMajor, CblasConjTrans, CblasNoTrans,
                 r1, r2, M,
                 CBLAS_SADDR(zone),  u1, ldu,
                                     u2, ldu,
                 CBLAS_SADDR(zzero), W,  r1 );
    flops += FLOPS_ZGEMM( r1, r2, M );

    /* Compute u2 = u2 - u1 ( u1^t u2 ) = u2 - u1 * W */
    cblas_zgemm( CblasColMajor, CblasNoTrans, CblasNoTrans,
                 M, r2, r1,
                 CBLAS_SADDR(mzone), u1, ldu,
                                     W,  r1,
                 CBLAS_SADDR(zone),  u2, ldu );
    flops += FLOPS_ZGEMM( M, r2, r1 );

    /* Update v1 = v1 + ( u1^t u2 ) v2 = v1 + W * v2 */
    cblas_zgemm( CblasColMajor, CblasNoTrans, CblasNoTrans,
                 r1, N, r2,
                 CBLAS_SADDR(zone), W,  r1,
                                    v2, ldv,
                 CBLAS_SADDR(zone), v1, ldv );
    flops += FLOPS_ZGEMM( r1, N, r2 );

#if !defined(PASTIX_LR_CGS1)
    /* Compute W = u1^t u2 */
    cblas_zgemm( CblasColMajor, CblasConjTrans, CblasNoTrans,
                 r1, r2, M,
                 CBLAS_SADDR(zone),  u1, ldu,
                                     u2, ldu,
                 CBLAS_SADDR(zzero), W,  r1 );
    flops += FLOPS_ZGEMM( r1, r2, M );

    /* Compute u2 = u2 - u1 ( u1^t u2 ) = u2 - u1 * W */
    cblas_zgemm( CblasColMajor, CblasNoTrans, CblasNoTrans,
                 M, r2, r1,
                 CBLAS_SADDR(mzone), u1, ldu,
                                     W,  r1,
                 CBLAS_SADDR(zone),  u2, ldu );
    flops += FLOPS_ZGEMM( M, r2, r1 );

    /* Update v1 = v1 + ( u1^t u2 ) v2 = v1 + W * v2 */
    cblas_zgemm( CblasColMajor, CblasNoTrans, CblasNoTrans,
                 r1, N, r2,
                 CBLAS_SADDR(zone), W,  r1,
                                    v2, ldv,
                 CBLAS_SADDR(zone), v1, ldv );
    flops += FLOPS_ZGEMM( r1, N, r2 );
#endif

#if defined(PASTIX_DEBUG_LR)
    if ( core_zlrdbg_check_orthogonality_AB( M, r1, r2, u1, ldu, u2, ldu ) != 0 ) {
        fprintf(stderr, "partialQR: u2 not correctly projected with u1\n" );
    }
#endif

    /* Compute u2 = Q * R */
    ret = LAPACKE_zgeqrf_work( LAPACK_COL_MAJOR, M, r2,
                               u2, ldu, tau, work, ldwork );
    assert( ret == 0 );
    flops += FLOPS_ZGEQRF( M, r2 );

    /* Compute v2' = R * v2 */
    cblas_ztrmm( CblasColMajor,
                 CblasLeft, CblasUpper,
                 CblasNoTrans, CblasNonUnit,
                 r2, N, CBLAS_SADDR(zone),
                 u2, ldu, v2, ldv);
    flops += FLOPS_ZTRMM( PastixLeft, r2, N );

    /* Generate the Q */
    ret = LAPACKE_zungqr_work( LAPACK_COL_MAJOR, M, r2, r2,
                               u2, ldu, tau, work, ldwork );
    assert( ret == 0 );
    flops += FLOPS_ZUNGQR( M, r2, r2 );

#if defined(PASTIX_DEBUG_LR)
    if ( core_zlrdbg_check_orthogonality_AB( M, r1, r2, u1, ldu, u2, ldu ) != 0 ) {
        fprintf(stderr, "partialQR: Final u2 not orthogonal to u1\n" );
    }
#endif

    free( W );

    (void)ret;
    (void)offx;
    (void)offy;

    return flops;
}

/**
 *******************************************************************************
 *
 * @brief Try to orthognalize the U part of the low-rank form, and update the V
 * part accordingly using CGS.
 *
 * This function considers a low-rank matrix resulting from the addition of two
 * matrices B += A, with A of smaller or equal size to B.
 * The product has the form: U * V^t
 *
 * The U part of the low-rank form must be orthognalized to get the smaller
 * possible rank during the rradd operation. This function perfoms this by
 * applying a full QR factorization on the U part.
 *
 * In that case, it takes benefit from the fact that U = [ u1, u2 ], and V = [
 * v1, v2 ] with u2 and v2 wich are matrices of respective size M2-by-r2, and
 * r2-by-N2, offset by offx and offy
 *
 * The steps are:
 *    - for each column of u2
 *       - Scaling of u2 with removal of the null columns
 *       - Orthogonalization of u2 relatively to u1
 *       - Remove the column if null
 *
 *******************************************************************************
 *
 * @param[in] M1
 *          The number of rows of the U matrix.
 *
 * @param[in] N1
 *          The number of columns of the U matrix.
 *
 * @param[in] M2
 *          The number of rows of the u2 part of the U matrix.
 *
 * @param[in] N2
 *          The number of columns of the v2 part of the V matrix.
 *
 * @param[in] r1
 *          The number of columns of the U matrix in the u1 part, and the number
 *          of rows of the V part in the v1 part.
 *
 * @param[inout] r2ptr
 *          The number of columns of the U matrix in the u2 part, and the number
 *          of rows of the V part in the v2 part. On exit, this rank is reduced
 *          y the number of null columns found in U.
 *
 * @param[in] offx
 *          The row offset of the matrix u2 in U.
 *
 * @param[in] offy
 *          The column offset of the matrix v2 in V.
 *
 * @param[inout] U
 *          The U matrix of size ldu -by- rank. On exit, the orthogonalized U.
 *
 * @param[in] ldu
 *          The leading dimension of the U matrix. ldu >= max(1, M)
 *
 * @param[inout] V
 *          The V matrix of size ldv -by- N.
 *          On exit, the updated V matrix.
 *
 * @param[in] ldv
 *          The leading dimension of the V matrix. ldv >= max(1, rank)
 *
 *******************************************************************************
 *
 * @return The number of flops required to perform the operation.
 *
 *******************************************************************************/
static inline pastix_fixdbl_t
core_zlrorthu_cgs( pastix_int_t M1,  pastix_int_t N1,
                   pastix_int_t M2,  pastix_int_t N2,
                   pastix_int_t r1, pastix_int_t *r2ptr,
                   pastix_int_t offx, pastix_int_t offy,
                   pastix_complex64_t *U, pastix_int_t ldu,
                   pastix_complex64_t *V, pastix_int_t ldv )
{
    pastix_int_t r2 = *r2ptr;
    pastix_complex64_t *u1 = U;
    pastix_complex64_t *u2 = U + r1 * ldu;
    pastix_complex64_t *v1 = V;
    pastix_complex64_t *v2 = V + r1;
    pastix_complex64_t *W;
    pastix_fixdbl_t flops = 0.0;
    pastix_int_t i, rank = r1 + r2;
    pastix_int_t ldwork = rank;
    double eps, norm;
    double norm_before, alpha;

    assert( M1 >= (M2 + offx) );
    assert( N1 >= (N2 + offy) );

    W     = malloc(ldwork * sizeof(pastix_complex64_t));
    eps   = LAPACKE_dlamch( 'e' );
    alpha = 1. / sqrt(2);

    /* Classical Gram-Schmidt */
    for (i=r1; i<rank; i++, u2 += ldu, v2++) {

        norm = cblas_dznrm2( M2, u2 + offx, 1 );
        if ( norm > ( M2 * eps ) ) {
            cblas_zdscal( M2, 1. / norm, u2 + offx,       1   );
            cblas_zdscal( N2, norm,      v2 + offy * ldv, ldv );
        }
        else {
            rank--; r2--;
            if ( i < rank ) {
                cblas_zswap( M2, u2 + offx, 1, U + rank * ldu + offx, 1 );
#if !defined(NDEBUG)
                memset( U + rank * ldu, 0,  M1 * sizeof(pastix_complex64_t) );
#endif

                cblas_zswap( N2, v2 + offy * ldv, ldv, V + offy * ldv + rank, ldv );

#if !defined(NDEBUG)
                LAPACKE_zlaset_work( LAPACK_COL_MAJOR, 'A', 1, N1,
                                     0., 0., V + rank, ldv );
#endif
                i--;
                u2-= ldu;
                v2--;
            }
#if !defined(NDEBUG)
            else {
                memset( u2, 0,  M1 * sizeof(pastix_complex64_t) );
                LAPACKE_zlaset_work( LAPACK_COL_MAJOR, 'A', 1, N1,
                                     0., 0., v2, ldv );
            }
#endif
            continue;
        }

        /* Compute W = u1^t u2 */
        cblas_zgemv( CblasColMajor, CblasConjTrans,
                     M2, i,
                     CBLAS_SADDR(zone),  u1+offx, ldu,
                                         u2+offx, 1,
                     CBLAS_SADDR(zzero), W,       1 );
        flops += FLOPS_ZGEMM( M2, i, 1 );

        /* Compute u2 = u2 - u1 ( u1^t u2 ) = u2 - u1 * W */
        cblas_zgemv( CblasColMajor, CblasNoTrans,
                     M1, i,
                     CBLAS_SADDR(mzone), u1, ldu,
                                         W,  1,
                     CBLAS_SADDR(zone),  u2, 1 );
        flops += FLOPS_ZGEMM( M1, i, 1 );

        /* Update v1 = v1 + ( u1^t u2 ) v2 = v1 + W * v2 */
        cblas_zgemm( CblasColMajor, CblasNoTrans, CblasNoTrans,
                     i, N1, 1,
                     CBLAS_SADDR(zone), W,  i,
                                        v2, ldv,
                     CBLAS_SADDR(zone), v1, ldv );
        flops += FLOPS_ZGEMM( i, N1, 1 );

        norm_before = cblas_dznrm2( i, W,  1 );
        norm        = cblas_dznrm2( M1, u2, 1 );

#if !defined(PASTIX_LR_CGS1)
        if ( norm <= (alpha * norm_before) ){
            /* Compute W = u1^t u2 */
            cblas_zgemv( CblasColMajor, CblasConjTrans,
                         M1, i,
                         CBLAS_SADDR(zone),  u1, ldu,
                                             u2, 1,
                         CBLAS_SADDR(zzero), W,  1 );
            flops += FLOPS_ZGEMM( M1, i, 1 );

            /* Compute u2 = u2 - u1 ( u1^t u2 ) = u2 - u1 * W */
            cblas_zgemv( CblasColMajor, CblasNoTrans,
                         M1, i,
                         CBLAS_SADDR(mzone), u1, ldu,
                                             W,  1,
                         CBLAS_SADDR(zone),  u2, 1 );
            flops += FLOPS_ZGEMM( M1, i, 1 );

            /* Update v1 = v1 + ( u1^t u2 ) v2 = v1 + W * v2 */
            cblas_zgemm( CblasColMajor, CblasNoTrans, CblasNoTrans,
                         i, N1, 1,
                         CBLAS_SADDR(zone), W,  i,
                                            v2, ldv,
                         CBLAS_SADDR(zone), v1, ldv );
            flops += FLOPS_ZGEMM( i, N1, 1 );

            norm = cblas_dznrm2( M1, u2, 1 );
        }
#endif

        if ( norm > M1 * eps ) {
            cblas_zdscal( M1, 1. / norm, u2, 1   );
            cblas_zdscal( N1, norm,      v2, ldv );
        }
        else {
            rank--; r2--;
            if ( i < rank ) {
                cblas_zswap( M1, u2, 1, U + rank * ldu, 1 );
                memset( U + rank * ldu, 0,  M1 * sizeof(pastix_complex64_t) );

                cblas_zswap( N1, v2, ldv, V + rank,     ldv );
                LAPACKE_zlaset_work( LAPACK_COL_MAJOR, 'A', 1, N1,
                                     0., 0., V + rank, ldv );
                i--;
                u2-= ldu;
                v2--;
            }
            else {
                memset( u2, 0,  M1 * sizeof(pastix_complex64_t) );
                LAPACKE_zlaset_work( LAPACK_COL_MAJOR, 'A', 1, N1,
                                     0., 0., v2, ldv );
            }
        }
    }
    free(W);

#if defined(PASTIX_DEBUG_LR)
    {
        u2 = U + r1 * ldu;
        if ( core_zlrdbg_check_orthogonality_AB( M1, r1, r2, u1, ldu, u2, ldu ) != 0 ) {
            fprintf(stderr, "cgs: Final u2 not orthogonal to u1\n" );
        }
    }
#endif

    *r2ptr = r2;

    (void)offy;
    (void)N2;
    return flops;
}

/**
 *******************************************************************************
 *
 * @brief Compute a rank-reavealing QR factorization.
 *
 * This routine is originated from the LAPACK kernels zgeqp3/zlaqps and was
 * modified by A. Buttari for MUMPS-BLR.
 *
 *******************************************************************************
 *
 * @param[in] m
 *          Number of rows of the matrix A.
 *
 * @param[in] n
 *          Number of columns of the matrix A.
 *
 * @param[in] A
 *          The matrix of dimension lda-by-n that need to be compressed
 *
 * @param[in] lda
 *          The leading dimension of the matrix A. lda >= max(1, m)
 *
 * @param[out] jpvt
 *          The array that describes the permutation of A
 *
 * @param[out] tau
 *          Contains scalar factors of the elementary reflectors for the matrix
 *          Q
 *
 * @param[in] work
 *          Workspace array
 *
 * @param[in] ldwork
 *          The leading dimension of work
 *
 * @param[in] rwork
 *          Workspace array used to store partial and exact column norms
 *
 * @param[in] tol
 *          The relative tolerance criteria. Computations are stopped when the
 *          norm of the residual matrix is lower than tol.
 *
 * @param[in] nb
 *          Blocking size for GEMM
 *
 * @param[in] maxrank
 *         Maximum rank of A. Computations are stopped when the rank exceeds
 *         maxrank
 *
 *******************************************************************************
 *
 * @return  This routine will return the rank of A (>=0) or an error (<0)
 *
 *******************************************************************************/
int
core_zrrqr( pastix_int_t m, pastix_int_t n,
            pastix_complex64_t *A, pastix_int_t lda,
            pastix_int_t *jpvt, pastix_complex64_t *tau,
            pastix_complex64_t *work, pastix_int_t ldwork,
            double *rwork,
            double tol, pastix_int_t nb, pastix_int_t maxrank )
{
    pastix_int_t minMN, ldf;
    pastix_int_t j, k, jb, itemp, lsticc, pvt;
    double temp, temp2, machine_prec;
    pastix_complex64_t akk;

    pastix_complex64_t *auxv, *f;

    /* Partial (VN1) and exact (VN2) column norms */
    double *VN1, *VN2;

    /* Number or rows of A that have been factorized */
    pastix_int_t offset = 0;

    /* Rank */
    pastix_int_t rk = 0;

#if !defined(NDEBUG)
    if (m < 0) {
        return -1;
    }
    if (n < 0) {
        return -2;
    }
    if (lda < pastix_imax(1, m)) {
        return -4;
    }
    if( ldwork < n) {
        return -8;
    }
#endif

    minMN = pastix_imin(m, n);
    if ((minMN == 0) || (tol <= 0.0)) {
        return 0;
    }

    VN1 = rwork;
    VN2 = rwork + n;

    auxv = work;
    f    = work + ldwork;
    ldf  = ldwork;

    /*
     * Initialize partial column norms. The first N elements of work store the
     * exact column norms.
     */
    /* TODO: call PLASMA/internal kernel */
    for (j=0; j<n; j++){
        VN1[j]  = cblas_dznrm2(m, A + j*lda, 1);
        VN2[j]  = VN1[j];
        jpvt[j] = j;
    }

    offset = 0;
    machine_prec = sqrt(LAPACKE_dlamch_work('e'));
    rk = 0;

    while ( rk < minMN ) {
        /* jb equivalent to kb in LAPACK xLAQPS: number of columns actually factorized */
        jb     = pastix_imin(nb, minMN-offset+1);
        lsticc = 0;

        /* column being factorized among jb */
        k = 0;

        while(k < jb) {

            rk = offset+k;

            /*
             * Rank is too large for compression
             */
            if (rk > maxrank) {
                return -1;
            }
            if (rk == minMN) {
                return rk;
            }

            pvt = rk + cblas_idamax(n-rk, VN1 + rk, 1);

            if (VN1[pvt] < tol) {
                double residual = cblas_dnrm2( n-rk, VN1 + rk, 1 );
                if (residual < tol) {
                    return rk;
                }
            }

            /*
             * Pivot is not within the current column: we swap
             */
            if (pvt != rk){
                assert( pvt < n );
                cblas_zswap( m, A + pvt * lda, 1,
                                A + rk  * lda, 1 );
                cblas_zswap( k, f + (pvt-offset), ldf,
                                f + k,            ldf );

                itemp     = jpvt[pvt];
                jpvt[pvt] = jpvt[rk];
                jpvt[rk]  = itemp;
                VN1[pvt]  = VN1[rk];
                VN2[pvt]  = VN2[rk];
            }

            /*
             * Apply previous Householder reflectors to column K
             * A(RK:M,RK) := A(RK:M,RK) - A(RK:M,OFFSET+1:RK-1)*F(K,1:K-1)**H
             */
            if (k > 0) {
#if defined(PRECISION_c) || defined(PRECISION_z)
                for (j=0; j<k; j++){
                    f[j * ldf + k] = conj(f[j * ldf + k]);
                }
#endif

                assert( (rk < n) && (rk < m) );
                cblas_zgemv( CblasColMajor, CblasNoTrans, m-rk, k,
                             CBLAS_SADDR(mzone), A + offset * lda + rk, lda,
                                                 f +                k,  ldf,
                             CBLAS_SADDR(zone),  A + rk     * lda + rk, 1 );

#if defined(PRECISION_c) || defined(PRECISION_z)
                for (j=0; j<k; j++){
                    f[j * ldf + k] = conj(f[j * ldf + k]);
                }
#endif
            }

            /*
             * Generate elementary reflector H(k).
             */
            if ((rk+1) < m) {
                LAPACKE_zlarfg(m-rk, A + rk * lda + rk, A + rk * lda + (rk+1), 1, tau + rk);
            }
            else{
                LAPACKE_zlarfg(1,    A + rk * lda + rk, A + rk * lda + rk,     1, tau + rk);
            }

            akk = A[rk * lda + rk];
            A[rk * lda + rk] = zone;

            /*
             * Compute Kth column of F:
             * F(K+1:N,K) := tau(K)*A(RK:M,K+1:N)**H*A(RK:M,K).
             */
            if ((rk+1) < n) {
                pastix_complex64_t alpha = tau[rk];
                cblas_zgemv( CblasColMajor, CblasConjTrans, m-rk, n-rk-1,
                             CBLAS_SADDR(alpha), A + (rk+1) * lda + rk,    lda,
                                                 A +  rk    * lda + rk,    1,
                             CBLAS_SADDR(zzero), f +  k     * ldf + k + 1, 1 );
            }

            /*
             * Padding F(1:K,K) with zeros.
             */
            memset( f + k * ldf, 0, k * sizeof( pastix_complex64_t ) );

            /*
             * Incremental updating of F:
             * F(1:N,K) := F(1:N-OFFSET,K) - tau(RK)*F(1:N,1:K-1)*A(RK:M,OFFSET+1:RK-1)**H*A(RK:M,RK).
             */
            if (k > 0) {
                pastix_complex64_t alpha = -tau[rk];
                cblas_zgemv( CblasColMajor, CblasConjTrans, m-rk, k,
                             CBLAS_SADDR(alpha), A + offset * lda + rk, lda,
                                                 A + rk     * lda + rk, 1,
                             CBLAS_SADDR(zzero), auxv,                  1 );

                cblas_zgemv( CblasColMajor, CblasNoTrans, n-offset, k,
                             CBLAS_SADDR(zone), f,           ldf,
                                                auxv,        1,
                             CBLAS_SADDR(zone), f + k * ldf, 1);
            }

            /*
             * Update the current row of A:
             * A(RK,RK+1:N) := A(RK,RK+1:N) - A(RK,OFFSET+1:RK)*F(K+1:N,1:K)**H.
             */
            if ((rk+1) < n) {
#if defined(PRECISION_c) || defined(PRECISION_z)
                cblas_zgemm( CblasColMajor, CblasNoTrans, CblasConjTrans,
                             1, n-rk-1, k+1,
                             CBLAS_SADDR(mzone), A + (offset) * lda + rk,    lda,
                                                 f +                  (k+1), ldf,
                             CBLAS_SADDR(zone),  A + (rk + 1) * lda + rk,    lda );
#else
                cblas_zgemv( CblasColMajor, CblasNoTrans, n-rk-1, k+1,
                             CBLAS_SADDR(mzone), f +                  (k+1), ldf,
                                                 A + (offset) * lda + rk,    lda,
                             CBLAS_SADDR(zone),  A + (rk + 1) * lda + rk,    lda );
#endif
            }

            /*
             * Update partial column norms.
             */
            if ((rk+1) < minMN) {
                for (j=rk+1; j<n; j++) {
                    if (VN1[j] != 0.0) {
                        /*
                         * NOTE: The following 4 lines follow from the analysis in
                         * Lapack Working Note 176.
                         */
                        temp  = cabs( A[j * lda + rk] ) / VN1[j];
                        temp2 = (1.0 + temp) * (1.0 - temp);
                        temp  = (temp2 > 0.0) ? temp2 : 0.0;

                        temp2 = temp * ((VN1[j] / VN2[j]) * ( VN1[j] / VN2[j]));
                        if (temp2 < machine_prec){
                            VN2[j] = (double)lsticc;
                            lsticc = j;
                        }
                        else{
                            VN1[j] = VN1[j] * sqrt(temp);
                        }
                    }
                }
            }

            A[rk * lda + rk] = akk;
            k++;

            if (lsticc != 0) break;
        }

        /*
         * Apply the block reflector to the rest of the matrix:
         * A(RK+1:M,RK+1:N) := A(RK+1:M,RK+1:N) -
         * A(RK+1:M,OFFSET+1:RK)*F(K+1:N-OFFSET,1:K)**H.
         */
        if ((rk+1) < minMN) {
            cblas_zgemm( CblasColMajor, CblasNoTrans, CblasConjTrans,
                         m-rk-1, n-rk-1, k,
                         CBLAS_SADDR(mzone), A + (offset) * lda + rk + 1, lda,
                                             f +                  k,      ldf,
                         CBLAS_SADDR(zone),  A + (rk+1)   * lda + rk + 1, lda );
        }

        /* Recomputation of difficult columns. */
        while (lsticc > 0) {
            assert(lsticc < n);
            itemp = (pastix_int_t) (VN2[lsticc]);

            VN1[lsticc] = cblas_dznrm2(m-rk-1, A + lsticc * lda + rk+1, 1 );

            /*
             * NOTE: The computation of VN1( LSTICC ) relies on the fact that
             * SNRM2 does not fail on vectors with norm below the value of
             * SQRT(DLAMCH('S'))
             */
            VN2[lsticc] = VN1[lsticc];
            lsticc = itemp;
        }

        rk++;
        offset = rk;
    }

    return rk;
}

/**
 *******************************************************************************
 *
 * @brief Convert a full rank matrix in a low rank matrix, using RRQR.
 *
 *******************************************************************************
 *
 * @param[in] tol
 *          The tolerance used as a criterai to eliminate information from the
 *          full rank matrix
 *
 * @param[in] rklimit
 *          The maximum rank to store the matrix in low-rank format. If
 *          -1, set to min(m, n) / PASTIX_LR_MINRATIO.
 *
 * @param[in] m
 *          Number of rows of the matrix A, and of the low rank matrix Alr.
 *
 * @param[in] n
 *          Number of columns of the matrix A, and of the low rank matrix Alr.
 *
 * @param[in] A
 *          The matrix of dimension lda-by-n that need to be compressed
 *
 * @param[in] lda
 *          The leading dimension of the matrix A. lda >= max(1, m)
 *
 * @param[out] Alr
 *          The low rank matrix structure that will store the low rank
 *          representation of A
 *
 *******************************************************************************/
pastix_fixdbl_t
core_zge2lr_rrqr( pastix_fixdbl_t tol, pastix_int_t rklimit,
                  pastix_int_t m, pastix_int_t n,
                  const void *Avoid, pastix_int_t lda,
                  pastix_lrblock_t *Alr )
{
    int                 ret;
    pastix_int_t        nb     = 32;
    pastix_int_t        ldwork = pastix_imax(m, n);
    pastix_complex64_t *work, *Acpy, *tau;
    double             *rwork;
    pastix_int_t       *jpvt;
    pastix_int_t        zsize, rsize;
    pastix_complex64_t *zwork;
    pastix_complex64_t *A = (pastix_complex64_t*)Avoid;
    pastix_fixdbl_t     flops;

    double norm = LAPACKE_zlange_work( LAPACK_COL_MAJOR, 'f', m, n,
                                       A, lda, NULL );

    /* work */
    zsize = (2 * nb + 1) * ldwork;
    /* Acpy */
    zsize += m * n;
    /* tau */
    zsize += n;

    /* rwork */
    rsize = 2* n;

    zwork = malloc( zsize * sizeof(pastix_complex64_t) + rsize * sizeof(double) );
    rwork = (double*)(zwork + zsize);

    work = zwork;
    Acpy = zwork + (2 * nb + 1) * ldwork;
    tau  = Acpy + m * n;

    MALLOC_INTERN( jpvt, n, pastix_int_t );

    /**
     * Allocate a temorary Low rank matrix
     */
    core_zlralloc( m, n, pastix_imin( m, n ), Alr );

    ret = LAPACKE_zlacpy_work( LAPACK_COL_MAJOR, 'A', m, n,
                               A, lda, Acpy, m );
    assert(ret == 0);

    rklimit = (rklimit == -1) ? core_get_rklimit( m, n ) : rklimit;
    rklimit = pastix_imin( rklimit, pastix_imin(m, n) );
    ret = core_zrrqr( m, n,
                      Acpy, m,
                      jpvt, tau,
                      work, ldwork,
                      rwork,
                      tol * norm, nb, rklimit );
    if (ret == -1) {
        flops = FLOPS_ZGEQRF( m, n );
    }
    else {
        flops = FLOPS_ZGEQRF( m, ret ) + FLOPS_ZUNMQR( m, n-ret, ret, PastixLeft );
    }

    /**
     * Resize the space used by the low rank matrix
     */
    core_zlrsze( 0, m, n, Alr, ret, -1, rklimit );

    /**
     * It was not interesting to compress, so we store the dense version in Alr
     */
    if ( Alr->rk == -1 ) {
        memFree_null(Alr->u);
        memFree_null(Alr->v);
        Alr->u     = malloc( m * n * sizeof(pastix_complex64_t) );
        Alr->rk    = -1;
        Alr->rkmax = m;

        ret = LAPACKE_zlacpy_work( LAPACK_COL_MAJOR, 'A', m, n,
                                   A, lda, Alr->u, Alr->rkmax );
        assert(ret == 0);
    }
    /**
     * We compute Q/R obtained thanks to core_zrrqr
     */
    else if (Alr->rk != 0) {
        pastix_int_t i;

        /* Temporary space to permute Alr->v */
        pastix_complex64_t *permQR;
        pastix_complex64_t *U = Alr->u;
        pastix_complex64_t *V = Alr->v;

        MALLOC_INTERN( permQR, n * Alr->rk, pastix_complex64_t );
        memset( permQR, 0, n * Alr->rk * sizeof(pastix_complex64_t) );

        ret = LAPACKE_zlacpy_work( LAPACK_COL_MAJOR, 'U', Alr->rk, n,
                                   Acpy, m,
                                   permQR, Alr->rk );
        assert(ret == 0);

        /* Permute V */
        for (i=0; i<n; i++){
            memcpy(V + jpvt[i] * Alr->rk,
                   permQR + i * Alr->rk,
                   Alr->rk * sizeof(pastix_complex64_t));
        }

        memFree_null( permQR );

        /* Compute Q factor on u */
        ret = LAPACKE_zlacpy_work( LAPACK_COL_MAJOR, 'A', m, Alr->rk,
                                   Acpy, m, U, m );
        assert(ret == 0);

        ret = LAPACKE_zungqr( LAPACK_COL_MAJOR, m, Alr->rk, Alr->rk,
                              U , m, tau );
        assert(ret == 0);
        flops += FLOPS_ZUNGQR( m, Alr->rk, Alr->rk );
    }

#if defined(PASTIX_DEBUG_LR)
    if ( Alr->rk > 0 ) {
        int rc = core_zlrdbg_check_orthogonality( m, Alr->rk, Alr->u, m );
        if (rc == 1) {
            fprintf(stderr, "Failed to compress a matrix and generate an othrogonal u\n" );
        }
    }
#endif
    memFree_null( zwork );
    memFree_null( jpvt );
    return flops;
}

/**
 *******************************************************************************
 *
 * @brief Add two LR structures A=(-u1) v1^T and B=u2 v2^T into u2 v2^T
 *
 *    u2v2^T - u1v1^T = (u2 u1) (v2 v1)^T
 *    Orthogonalize (u2 u1) = (u2, u1 - u2(u2^T u1)) * (I u2^T u1)
 *                                                     (0    I   )
 *    Compute RRQR decomposition of (I u2^T u1) * (v2 v1)^T
 *                                  (0    I   )
 *
 *******************************************************************************
 *
 * @param[in] lowrank
 *          The structure with low-rank parameters.
 *
 * @param[in] transA1
 *         @arg PastixNoTrans:  No transpose, op( A ) = A;
 *         @arg PastixTrans:  Transpose, op( A ) = A';
 *
 * @param[in] alpha
 *          alpha * A is add to B
 *
 * @param[in] M1
 *          The number of rows of the matrix A.
 *
 * @param[in] N1
 *          The number of columns of the matrix A.
 *
 * @param[in] A
 *          The low-rank representation of the matrix A.
 *
 * @param[in] M2
 *          The number of rows of the matrix B.
 *
 * @param[in] N2
 *          The number of columns of the matrix B.
 *
 * @param[in] B
 *          The low-rank representation of the matrix B.
 *
 * @param[in] offx
 *          The horizontal offset of A with respect to B.
 *
 * @param[in] offy
 *          The vertical offset of A with respect to B.
 *
 *******************************************************************************
 *
 * @return  The new rank of u2 v2^T or -1 if ranks are too large for
 *          recompression
 *
 *******************************************************************************/
pastix_fixdbl_t
core_zrradd_rrqr( const pastix_lr_t *lowrank, pastix_trans_t transA1, const void *alphaptr,
                  pastix_int_t M1, pastix_int_t N1, const pastix_lrblock_t *A,
                  pastix_int_t M2, pastix_int_t N2,       pastix_lrblock_t *B,
                  pastix_int_t offx, pastix_int_t offy)
{
    pastix_int_t rankA, rank, M, N, minV;
    pastix_int_t i, ret, new_rank, rklimit;
    pastix_int_t ldau, ldav, ldbu, ldbv, ldu, ldv;
    pastix_complex64_t *u1u2, *v1v2, *u;
    pastix_complex64_t *zbuf, *tauV;
    size_t wzsize;
    double norm;
    double tol = lowrank->tolerance;

    /* RRQR parameters / workspace */
    pastix_int_t        nb = 32;
    pastix_int_t        ldwork;
    pastix_int_t       *jpvt;
    pastix_complex64_t *zwork;
    double             *rwork;
    pastix_complex64_t  alpha = *((pastix_complex64_t*)alphaptr);
    pastix_fixdbl_t     flops, total_flops = 0.;

#if defined(PASTIX_DEBUG_LR)
    if ( B->rk > 0 ) {
        int rc = core_zlrdbg_check_orthogonality( M2, B->rk, B->u, M2 );
        if (rc == 1) {
            fprintf(stderr, "Failed to have B->u orthogonal in entry of rradd\n" );
        }
    }
#endif

    rankA = (A->rk == -1) ? pastix_imin(M1, N1) : A->rk;
    rank  = rankA + B->rk;
    M = pastix_imax(M2, M1);
    N = pastix_imax(N2, N1);

    minV = pastix_imin(N, rank);

    assert(M2 == M && N2 == N);
    assert(B->rk != -1);

    assert( A->rk <= A->rkmax);
    assert( B->rk <= B->rkmax);

    if ( ((M1 + offx) > M2) ||
         ((N1 + offy) > N2) )
    {
        errorPrint("Dimensions are not correct");
        assert(0 /* Incorrect dimensions */);
        return total_flops;
    }

    /*
     * A is rank null, nothing to do
     */
    if ( A->rk == 0 ) {
        return total_flops;
    }

    /*
     * Let's handle case where B is a null matrix
     *   B = alpha A
     */
    if ( B->rk == 0 ) {
        core_zlrcpy( lowrank, transA1, alpha,
                     M1, N1, A, M2, N2, B,
                     offx, offy );
        return total_flops;
    }

    /*
     * The rank is too big, let's try to compress
     */
    if ( rank > pastix_imin( M, N ) ) {
        assert(0);
    }

    /*
     * Let's define leading dimensions
     */
    ldau = (A->rk == -1) ? A->rkmax : M1;
    ldav = (transA1 == PastixNoTrans) ? A->rkmax : N1;
    ldbu = M;
    ldbv = B->rkmax;
    ldu = M;
    ldv = rank;

    /*
     * Let's compute the size of the workspace
     */
    /* u1u2 and v1v2 */
    wzsize = (M+N) * rank;
    /* tauV */
    wzsize += minV;

    /* RRQR workspaces */
    ldwork = pastix_imax(rank, N);
    wzsize += (2 * nb + 1) * ldwork;

    zbuf = malloc( wzsize * sizeof(pastix_complex64_t) + 2 * pastix_imax(rank, N) * sizeof(double) );

    u1u2  = zbuf;
    v1v2  = u1u2 + ldu * rank;
    tauV  = v1v2 + N * ldv;
    zwork = tauV + rank;

    rwork = (double*)(zbuf + wzsize);

    /*
     * Concatenate U2 and U1 in u1u2
     *  [ u2  0  ]
     *  [ u2  u1 ]
     *  [ u2  0  ]
     */
    core_zlrconcatenate_u( alpha,
                           M1, N1, A,
                           M2,     B,
                           offx, u1u2 );

    /*
     * Concatenate V2 and V1 in v1v2
     *  [ v2^h v2^h v2^h ]
     *  [ 0    v1^h 0    ]
     */
    core_zlrconcatenate_v( transA1, alpha,
                           M1, N1, A,
                               N2, B,
                           offy, v1v2 );

    /*
     * Orthogonalize [u2, u1]
     * [u2, u1] = [u2, u1 - u2(u2Tu1)] * (I u2Tu1)
     *                                   (0   I  )
     */

    /* We do not care is A was integrated into v1v2 */
    if (rankA != 0) {

        kernel_trace_start_lvl2( PastixKernelLvl2_LR_add2C_rradd_orthogonalize );

        if ( pastix_lr_ortho == PastixCompressOrthoQR ){
            flops = core_zlrorthu_fullqr( M, N, B->rk + rankA,
                                          u1u2, ldu, v1v2, ldv );
        }
        else if ( pastix_lr_ortho == PastixCompressOrthoPartialQR ){
            flops = core_zlrorthu_partialqr( M, N, B->rk, &rankA, offx, offy,
                                             u1u2, ldu, v1v2, ldv );
        }
        else{
            flops = core_zlrorthu_cgs( M2, N2, M1, N1, B->rk, &rankA, offx, offy,
                                       u1u2, ldu, v1v2, ldv );
        }
        kernel_trace_stop_lvl2( flops );

        total_flops += flops;
    }

    rank = B->rk + rankA;

    if (rankA == 0) {
        /*
         * The final B += A fit in B
         * Lets copy and return
         */
        LAPACKE_zlacpy_work( LAPACK_COL_MAJOR, 'A', M, B->rk, u1u2, ldu, B->u, ldbu );
        LAPACKE_zlacpy_work( LAPACK_COL_MAJOR, 'A', B->rk, N, v1v2, ldv, B->v, ldbv );

        free(zbuf);
        return total_flops;
    }

    MALLOC_INTERN( jpvt, pastix_imax(rank, N), pastix_int_t );

    norm = LAPACKE_zlange_work( LAPACK_COL_MAJOR, 'f', rank, N,
                                v1v2, ldv, NULL );

    /*
     * Perform RRQR factorization on (I u2Tu1) v1v2 = (Q2 R2)
     *                               (0   I  )
     */
    rklimit = core_get_rklimit( M, N );

    kernel_trace_start_lvl2( PastixKernelLvl2_LR_add2C_rradd_recompression );
    new_rank = core_zrrqr( rank, N,
                           v1v2, ldv,
                           jpvt, tauV,
                           zwork, ldwork,
                           rwork,
                           tol * norm, nb, pastix_imin(rklimit, rank) );
    flops = (new_rank == -1) ? FLOPS_ZGEQRF( rank, N )
        :                     (FLOPS_ZGEQRF( rank, new_rank ) +
                               FLOPS_ZUNMQR( rank, N-new_rank, new_rank, PastixLeft ));
    kernel_trace_stop_lvl2_rank( flops, new_rank );
    total_flops += flops;

    /*
     * First case: The rank is too big, so we decide to uncompress the result
     */
    if ( (new_rank > rklimit) ||
         (new_rank == -1) )
    {
        pastix_lrblock_t Bbackup = *B;

        core_zlralloc( M, N, -1, B );
        u = B->u;

        /* Uncompress B */
        flops = FLOPS_ZGEMM( M, N, Bbackup.rk );
        kernel_trace_start_lvl2( PastixKernelLvl2_LR_add2C_uncompress );
        cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                    M, N, Bbackup.rk,
                    CBLAS_SADDR(zone),  Bbackup.u, ldbu,
                                        Bbackup.v, ldbv,
                    CBLAS_SADDR(zzero), u, M );
        kernel_trace_stop_lvl2( flops );
        total_flops += flops;

        /* Add A into it */
        if ( A->rk == -1 ) {
            flops = 2 * M1 * N1;
            kernel_trace_start_lvl2( PastixKernelLvl2_FR_GEMM );
            core_zgeadd( transA1, M1, N1,
                         alpha, A->u, ldau,
                         zone, u + offy * M + offx, M);
            kernel_trace_stop_lvl2( flops );
        }
        else {
            flops = FLOPS_ZGEMM( M1, N1, A->rk );
            kernel_trace_start_lvl2( PastixKernelLvl2_FR_GEMM );
            cblas_zgemm(CblasColMajor, CblasNoTrans, (CBLAS_TRANSPOSE)transA1,
                        M1, N1, A->rk,
                        CBLAS_SADDR(alpha), A->u, ldau,
                                            A->v, ldav,
                        CBLAS_SADDR(zone), u + offy * M + offx, M);
            kernel_trace_stop_lvl2( flops );
        }
        total_flops += flops;
        core_zlrfree(&Bbackup);
        memFree_null(zbuf);
        memFree_null(jpvt);
        return total_flops;
    }
    else if ( new_rank == 0 ) {
        core_zlrfree(B);
        memFree_null(zbuf);
        memFree_null(jpvt);
        return total_flops;
    }

    /*
     * We need to reallocate the buffer to store the new compressed version of B
     * because it wasn't big enough
     */
    ret = core_zlrsze( 0, M, N, B, new_rank, -1, -1 );
    assert( ret != -1 );
    assert( B->rkmax >= new_rank );
    assert( B->rkmax >= B->rk    );

    ldbv = B->rkmax;

    /* B->v = P v1v2 */
    {
        pastix_complex64_t *tmpV;
        pastix_int_t lm;

        memset(B->v, 0, N * ldbv * sizeof(pastix_complex64_t));
        tmpV = B->v;
        for (i=0; i<N; i++){
            lm = pastix_imin( new_rank, i+1 );
            memcpy(tmpV + jpvt[i] * ldbv,
                   v1v2 + i       * ldv,
                   lm * sizeof(pastix_complex64_t));
        }
    }

    /* Compute Q2 factor */
    {
        flops = FLOPS_ZUNGQR( rank, new_rank, new_rank )
            +   FLOPS_ZGEMM( M, new_rank, rank );

        kernel_trace_start_lvl2( PastixKernelLvl2_LR_add2C_rradd_computeNewU );
        ret = LAPACKE_zungqr( LAPACK_COL_MAJOR, rank, new_rank, new_rank,
                              v1v2, ldv, tauV );
        assert(ret == 0);

        cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                    M, new_rank, rank,
                    CBLAS_SADDR(zone),  u1u2, ldu,
                                        v1v2, ldv,
                    CBLAS_SADDR(zzero), B->u, ldbu);
        kernel_trace_stop_lvl2( flops );
        total_flops += flops;
    }

    memFree_null(zbuf);
    memFree_null(jpvt);

    (void)ret;
    return total_flops;
}
