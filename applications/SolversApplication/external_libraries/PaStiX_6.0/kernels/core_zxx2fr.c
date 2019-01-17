/**
 *
 * @file core_zxx2fr.c
 *
 * PaStiX low-rank kernel routines that form the product of two matrices A and B
 * into a low-rank form for an update on a full rank matrix.
 *
 * @copyright 2016-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.1
 * @author Mathieu Faverge
 * @author Gregoire Pichon
 * @author Pierre Ramet
 * @date 2018-07-16
 * @precisions normal z -> c d s
 *
 **/
#include "common.h"
#include <cblas.h>
#include "pastix_zlrcores.h"
#include "kernels_trace.h"
#ifndef DOXYGEN_SHOULD_SKIP_THIS
static pastix_complex64_t zone  =  1.0;
static pastix_complex64_t zzero =  0.0;
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/**
 *******************************************************************************
 *
 * @brief Perform the full-rank operation C = alpha * op(A) * op(B) + beta C
 *
 *******************************************************************************
 *
 * @param[inout] params
 *          The LRMM structure that stores all the parameters used in the LRMM
 *          functions family.
 *          On exit, the C matrix contains the product AB aligned with its own
 *          dimensions.
 *          @sa core_zlrmm_t
 *
 *******************************************************************************
 *
 * @return The number of flops required to perform the operation.
 *
 *******************************************************************************/
pastix_fixdbl_t
core_zfrfr2fr( core_zlrmm_t *params )
{
    pastix_int_t ldau, ldbu, ldcu;
    pastix_complex64_t *Cptr;
    pastix_fixdbl_t flops;
    PASTE_CORE_ZLRMM_PARAMS( params );
    ldau = (transA == PastixNoTrans) ? M : K;
    ldbu = (transB == PastixNoTrans) ? K : N;
    ldcu = Cm;

    Cptr  = C->u;
    Cptr += ldcu * offy + offx;

    pastix_atomic_lock( lock );
    assert( C->rk == -1 ); /* Check that C has not changed due to parallelism */

    /*
     * Everything is full rank we apply directly a GEMM
     */
    cblas_zgemm( CblasColMajor, (CBLAS_TRANSPOSE)transA, (CBLAS_TRANSPOSE)transB,
                 M, N, K,
                 CBLAS_SADDR(alpha), A->u, ldau,
                                     B->u, ldbu,
                 CBLAS_SADDR(beta),  Cptr, ldcu );
    flops = FLOPS_ZGEMM( M, N, K );

    pastix_atomic_unlock( lock );

    PASTE_CORE_ZLRMM_VOID;
    return flops;
}

/**
 *******************************************************************************
 *
 * @brief Perform the operation C = alpha * op(A) * op(B) + beta C, with A and C
 * full-rank and B low-rank.
 *
 *******************************************************************************
 *
 * @param[inout] params
 *          The LRMM structure that stores all the parameters used in the LRMM
 *          functions family.
 *          On exit, the C matrix contains the product AB aligned with its own
 *          dimensions.
 *          @sa core_zlrmm_t
 *
 *******************************************************************************
 *
 * @return The number of flops required to perform the operation.
 *
 *******************************************************************************/
pastix_fixdbl_t
core_zfrlr2fr( core_zlrmm_t *params )
{
    PASTE_CORE_ZLRMM_PARAMS( params );
    pastix_complex64_t *Cptr;
    pastix_int_t ldau, ldbu, ldbv, ldcu;
    pastix_fixdbl_t flops1 = FLOPS_ZGEMM( M, B->rk, K ) + FLOPS_ZGEMM( M, N, B->rk );
    pastix_fixdbl_t flops2 = FLOPS_ZGEMM( K, N, B->rk ) + FLOPS_ZGEMM( M, N, K     );
    pastix_fixdbl_t flops;
    int allocated = 0;

    ldau = (transA == PastixNoTrans) ? M : K;
    ldbu = (transB == PastixNoTrans) ? K : N;
    ldbv = ( B->rk == -1 ) ? -1 : B->rkmax;

    ldcu = Cm;
    Cptr = C->u;
    Cptr += ldcu * offy + offx;

    /*
     *  A(M-by-K) * B( N-by-rb x rb-by-K )^t
     */
    if ( flops1 <= flops2 ) {
        if ( (work = core_zlrmm_getws( params, M * B->rk )) == NULL ) {
            work = malloc( M * B->rk * sizeof(pastix_complex64_t) );
            allocated = 1;
        }

        /*
         *  (A * Bv) * Bu^t
         */
        cblas_zgemm( CblasColMajor, (CBLAS_TRANSPOSE)transA, (CBLAS_TRANSPOSE)transB,
                     M, B->rk, K,
                     CBLAS_SADDR(zone),  A->u, ldau,
                                         B->v, ldbv,
                     CBLAS_SADDR(zzero), work, M );

        pastix_atomic_lock( lock );
        assert( C->rk == -1 ); /* Check that C has not changed due to parallelism */
        cblas_zgemm( CblasColMajor, CblasNoTrans, (CBLAS_TRANSPOSE)transB,
                     M, N, B->rk,
                     CBLAS_SADDR(alpha), work, M,
                                         B->u, ldbu,
                     CBLAS_SADDR(beta),  Cptr, ldcu );
        flops = flops1;
        pastix_atomic_unlock( lock );
    }
    else {
        if ( (work = core_zlrmm_getws( params, K * N )) == NULL ) {
            work = malloc( K * N * sizeof(pastix_complex64_t) );
            allocated = 1;
        }

        /*
         *  A * (Bu * Bv^t)^t
         */
        cblas_zgemm( CblasColMajor, CblasNoTrans, CblasNoTrans,
                     K, N, B->rk,
                     CBLAS_SADDR(zone),  B->u, ldbu,
                                         B->v, ldbv,
                     CBLAS_SADDR(zzero), work, K );

        pastix_atomic_lock( lock );
        assert( C->rk == -1 ); /* Check that C has not changed due to parallelism */
        cblas_zgemm( CblasColMajor, (CBLAS_TRANSPOSE)transA, (CBLAS_TRANSPOSE)transB,
                     M, N, K,
                     CBLAS_SADDR(alpha), A->u, ldau,
                                         work, K,
                     CBLAS_SADDR(beta),  Cptr, ldcu );

        flops = flops2;
        pastix_atomic_unlock( lock );
    }

    if ( allocated ) {
        free( work );
    }
    PASTE_CORE_ZLRMM_VOID;
    return flops;
}

/**
 *******************************************************************************
 *
 * @brief Perform the operation C = alpha * op(A) * op(B) + beta C, with B and C
 * full-rank and A low-rank.
 *
 *******************************************************************************
 *
 * @param[inout] params
 *          The LRMM structure that stores all the parameters used in the LRMM
 *          functions family.
 *          On exit, the C matrix contains the product AB aligned with its own
 *          dimensions.
 *          @sa core_zlrmm_t
 *
 *******************************************************************************
 *
 * @return The number of flops required to perform the operation.
 *
 *******************************************************************************/
pastix_fixdbl_t
core_zlrfr2fr( core_zlrmm_t *params )
{
    PASTE_CORE_ZLRMM_PARAMS( params );
    pastix_complex64_t *Cptr;
    pastix_int_t ldau, ldav, ldbu, ldcu;
    pastix_fixdbl_t flops1 = FLOPS_ZGEMM( A->rk, N, K ) + FLOPS_ZGEMM( M, N, A->rk );
    pastix_fixdbl_t flops2 = FLOPS_ZGEMM( M, K, A->rk ) + FLOPS_ZGEMM( M, N, K     );
    pastix_fixdbl_t flops;
    int allocated = 0;

    ldau = (transA == PastixNoTrans) ? M : K;
    ldav = ( A->rk == -1 ) ? -1 : A->rkmax;
    ldbu = (transB == PastixNoTrans) ? K : N;

    ldcu = Cm;
    Cptr = C->u;
    Cptr += ldcu * offy + offx;

    /*
     *  A( M-by-ra x ra-by-K ) * B(N-by-K)^t
     */
    if ( flops1 <= flops2 ) {
        if ( (work = core_zlrmm_getws( params, A->rk * N )) == NULL ) {
            work = malloc( A->rk * N * sizeof(pastix_complex64_t) );
            allocated = 1;
        }

        /*
         *  Au * (Av^t * B^t)
         */
        cblas_zgemm( CblasColMajor, CblasNoTrans, (CBLAS_TRANSPOSE)transB,
                     A->rk, N, K,
                     CBLAS_SADDR(zone),  A->v, ldav,
                                         B->u, ldbu,
                     CBLAS_SADDR(zzero), work, A->rk );

        pastix_atomic_lock( lock );
        assert( C->rk == -1 ); /* Check that C has not changed due to parallelism */
        cblas_zgemm( CblasColMajor, CblasNoTrans, CblasNoTrans,
                     M, N, A->rk,
                     CBLAS_SADDR(alpha), A->u, ldau,
                                         work, A->rk,
                     CBLAS_SADDR(beta),  Cptr, ldcu );

        flops = flops1;
        pastix_atomic_unlock( lock );
    }
    else {
        if ( (work = core_zlrmm_getws( params, M * K )) == NULL ) {
            work = malloc( M * K * sizeof(pastix_complex64_t) );
            allocated = 1;
        }

        /*
         *  (Au * Av^t) * B^t
         */
        cblas_zgemm( CblasColMajor, CblasNoTrans, CblasNoTrans,
                     M, K, A->rk,
                     CBLAS_SADDR(zone),  A->u, ldau,
                                         A->v, ldav,
                     CBLAS_SADDR(zzero), work, M );

        pastix_atomic_lock( lock );
        assert( C->rk == -1 ); /* Check that C has not changed due to parallelism */
        cblas_zgemm( CblasColMajor, (CBLAS_TRANSPOSE)transA, (CBLAS_TRANSPOSE)transB,
                     M, N, K,
                     CBLAS_SADDR(alpha), work, M,
                                         B->u, ldbu,
                     CBLAS_SADDR(beta),  Cptr, ldcu );

        flops = flops2;
        pastix_atomic_unlock( lock );
    }

    if ( allocated ) {
        free( work );
    }
    PASTE_CORE_ZLRMM_VOID;
    return flops;
}

/**
 *******************************************************************************
 *
 * @brief Perform the operation C = alpha * op(A) * op(B) + beta C, with A and B
 * low-rank and C full-rank.
 *
 *******************************************************************************
 *
 * @param[inout] params
 *          The LRMM structure that stores all the parameters used in the LRMM
 *          functions family.
 *          On exit, the C matrix contains the product AB aligned with its own
 *          dimensions.
 *          @sa core_zlrmm_t
 *
 *******************************************************************************
 *
 * @return The number of flops required to perform the operation.
 *
 *******************************************************************************/
pastix_fixdbl_t
core_zlrlr2fr( core_zlrmm_t *params )
{
    PASTE_CORE_ZLRMM_PARAMS( params );
    pastix_complex64_t *Cptr;
    pastix_int_t        ldcu;
    pastix_lrblock_t    AB;
    pastix_trans_t      trans = PastixNoTrans;
    int                 infomask = 0;
    pastix_fixdbl_t     flops;

    ldcu = Cm;
    Cptr = C->u;
    Cptr += ldcu * offy + offx;

    flops = core_zlrlr2lr( params, &AB, &infomask );
    assert( AB.rk != -1 );
    assert( AB.rkmax != -1 );

    if ( infomask & PASTIX_LRM3_TRANSB ) {
        trans = transB;
    }

    if ( AB.rk > 0 ) {
        pastix_int_t ldabv = (trans == PastixNoTrans) ? AB.rkmax : N;

        pastix_atomic_lock( lock );
        assert( C->rk == -1 ); /* Check that C has not changed due to parallelism */

        cblas_zgemm( CblasColMajor, CblasNoTrans, (CBLAS_TRANSPOSE)trans,
                     M, N, AB.rk,
                     CBLAS_SADDR(alpha), AB.u, M,
                                         AB.v, ldabv,
                     CBLAS_SADDR(beta),  Cptr, ldcu );
        flops = FLOPS_ZGEMM( M, N, AB.rk );
        pastix_atomic_unlock( lock );
    }

    /* Free memory from zlrm3 */
    if ( infomask & PASTIX_LRM3_ALLOCU ) {
        free(AB.u);
    }
    if ( infomask & PASTIX_LRM3_ALLOCV ) {
        free(AB.v);
    }

    PASTE_CORE_ZLRMM_VOID;
    return flops;
}
