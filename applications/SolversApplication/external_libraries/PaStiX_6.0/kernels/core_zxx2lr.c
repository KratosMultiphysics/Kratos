/**
 *
 * @file core_zxx2lr.c
 *
 * PaStiX low-rank kernel routines that form the product of two matrices A and B
 * into a low-rank form for an update on a null or low-rank matrix.
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
#include "flops.h"
#include "blend/solver.h"
#include "pastix_zlrcores.h"

#ifndef DOXYGEN_SHOULD_SKIP_THIS
static pastix_complex64_t zone  =  1.0;
static pastix_complex64_t zzero =  0.0;
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/**
 *******************************************************************************
 *
 * @brief Perform the operation AB = op(A) * op(B), with A and B full-rank and AB
 * low-rank.
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
 * @param[inout] AB
 *          The low-rank structure of the AB matrix in which to store the AB product.
 *
 * @param[inout] infomask
 *          The mask of informations returned by the core_zxx2lr() functions.
 *          - If AB.u is orthogonal on exit, then PASTIX_LRM3_ORTHOU is set.
 *          - If AB.u is allocated, then PASTIX_LRM3_ALLOCU is set.
 *          - If AB.v is allocated, then PASTIX_LRM3_ALLOCV is set.
 *          - If AB.v is inistialized as one of the given pointer and op(B) is not
 *          applyed, then PASTIX_LRM3_TRANSB is set.
 *
 * @param[in] Kmax
 *          The maximum K value for which the AB product is contructed as AB.u =
 *          A, and AB.v = B
 *
 *******************************************************************************
 *
 * @return The number of flops required to perform the operation.
 *
 *******************************************************************************/
pastix_fixdbl_t
core_zfrfr2lr( core_zlrmm_t     *params,
               pastix_lrblock_t *AB,
               int              *infomask,
               pastix_int_t      Kmax )
{
    PASTE_CORE_ZLRMM_PARAMS( params );
    pastix_int_t ldau, ldbu;
    pastix_fixdbl_t flops = 0.0;

    ldau = (transA == PastixNoTrans) ? M : K;
    ldbu = (transB == PastixNoTrans) ? K : N;

    /*
     * Everything is full rank
     */
    if ( K < Kmax ) {
        /*
         * Let's build a low-rank matrix of rank K
         */
        AB->rk = K;
        AB->rkmax = K;
        AB->u = A->u;
        AB->v = B->u;
        *infomask |= PASTIX_LRM3_TRANSB;
    }
    else {
        /*
         * Let's compute the product to form a full-rank matrix of rank
         * pastix_imin( M, N )
         */
        if ( (work = core_zlrmm_getws( params, M * N )) == NULL ) {
            work = malloc( M * N * sizeof(pastix_complex64_t) );
            *infomask |= PASTIX_LRM3_ALLOCU;
        }
        AB->rk = -1;
        AB->rkmax = M;
        AB->u = work;
        AB->v = NULL;

        cblas_zgemm( CblasColMajor, CblasNoTrans, (CBLAS_TRANSPOSE)transB,
                     M, N, K,
                     CBLAS_SADDR(zone),  A->u,  ldau,
                                         B->u,  ldbu,
                     CBLAS_SADDR(zzero), AB->u, M );
        flops = FLOPS_ZGEMM( M, N, K );
    }

    PASTE_CORE_ZLRMM_VOID;
    return flops;
}

/**
 *******************************************************************************
 *
 * @brief Perform the operation AB = op(A) * op(B), with A full-rank and B and AB
 * low-rank.
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
 * @param[inout] AB
 *          The low-rank structure of the AB matrix in which to store the AB product.
 *
 * @param[inout] infomask
 *          The mask of informations returned by the core_zxx2lr() functions.
 *          - If AB.u is orthogonal on exit, then PASTIX_LRM3_ORTHOU is set.
 *          - If AB.u is allocated, then PASTIX_LRM3_ALLOCU is set.
 *          - If AB.v is allocated, then PASTIX_LRM3_ALLOCV is set.
 *          - If AB.v is inistialized as one of the given pointer and op(B) is not
 *          applyed, then PASTIX_LRM3_TRANSB is set.
 *
 * @param[in] Brkmin
 *          Threshold for which B->rk is considered as the final rank of AB
 *
 *******************************************************************************
 *
 * @return The number of flops required to perform the operation.
 *
 *******************************************************************************/
pastix_fixdbl_t
core_zfrlr2lr( core_zlrmm_t     *params,
               pastix_lrblock_t *AB,
               int              *infomask,
               pastix_int_t      Brkmin )
{
    PASTE_CORE_ZLRMM_PARAMS( params );
    pastix_int_t ldau, ldbu, ldbv;
    pastix_fixdbl_t flops;

    ldau = (transA == PastixNoTrans) ? M : K;
    ldbu = (transB == PastixNoTrans) ? K : N;
    ldbv = ( B->rk == -1 ) ? -1 : B->rkmax;

    /*
     *  A(M-by-K) * B( N-by-rb x rb-by-K )^t
     */
    if ( B->rk > Brkmin ) {
        /*
         * We are in a similar case to the _Cfr function, and we
         * choose the optimal number of flops.
         */
        pastix_fixdbl_t flops1 = FLOPS_ZGEMM( M, B->rk, K ) + FLOPS_ZGEMM( M, N, B->rk );
        pastix_fixdbl_t flops2 = FLOPS_ZGEMM( K, N, B->rk ) + FLOPS_ZGEMM( M, N, K     );
        pastix_complex64_t *tmp;

        AB->rk    = -1;
        AB->rkmax = M;
        AB->v     = NULL;

        if ( flops1 <= flops2 ) {
            if ( (work = core_zlrmm_getws( params, M * B->rk + M * N )) == NULL ) {
                work = malloc( (M * B->rk + M * N) * sizeof(pastix_complex64_t) );
                *infomask |= PASTIX_LRM3_ALLOCU;
            }

            /* AB->u will be destroyed later */
            AB->u = work;
            tmp   = work + M * N;

            /*
             *  (A * Bv) * Bu^t
             */
            cblas_zgemm( CblasColMajor, (CBLAS_TRANSPOSE)transA, (CBLAS_TRANSPOSE)transB,
                         M, B->rk, K,
                         CBLAS_SADDR(zone),  A->u, ldau,
                                             B->v, ldbv,
                         CBLAS_SADDR(zzero), tmp,  M );

            cblas_zgemm( CblasColMajor, CblasNoTrans, (CBLAS_TRANSPOSE)transB,
                         M, N, B->rk,
                         CBLAS_SADDR(zone),  tmp,   M,
                                             B->u,  ldbu,
                         CBLAS_SADDR(zzero), AB->u, M );

            flops = flops1;
        }
        else {
            if ( (work = core_zlrmm_getws( params, K * N + M * N )) == NULL ) {
                work = malloc( (K * N + M * N) * sizeof(pastix_complex64_t) );
                *infomask |= PASTIX_LRM3_ALLOCU;
            }

            /* AB->u will be destroyed later */
            AB->u = work;
            tmp   = work + M * N;

            /*
             *  A * (Bu * Bv^t)^t
             */
            cblas_zgemm( CblasColMajor, CblasNoTrans, CblasNoTrans,
                         K, N, B->rk,
                         CBLAS_SADDR(zone),  B->u, ldbu,
                                             B->v, ldbv,
                         CBLAS_SADDR(zzero), tmp,  K );

            cblas_zgemm( CblasColMajor, (CBLAS_TRANSPOSE)transA, (CBLAS_TRANSPOSE)transB,
                         M, N, K,
                         CBLAS_SADDR(zone),  A->u,  ldau,
                                             tmp,   K,
                         CBLAS_SADDR(zzero), AB->u, M );

            flops = flops2;
        }
    }
    else {
        /*
         * B->rk is the smallest rank
         */
        AB->rk    = B->rk;
        AB->rkmax = B->rkmax;
        AB->v     = B->u;
        *infomask |= PASTIX_LRM3_TRANSB;

        if ( (work = core_zlrmm_getws( params, M * B->rk )) == NULL ) {
            work = malloc( M * B->rk * sizeof(pastix_complex64_t) );
            *infomask |= PASTIX_LRM3_ALLOCU;
        }
        AB->u = work;

        cblas_zgemm( CblasColMajor, CblasNoTrans, (CBLAS_TRANSPOSE)transB,
                     M, B->rk, K,
                     CBLAS_SADDR(zone),  A->u,  ldau,
                                         B->v,  ldbv,
                     CBLAS_SADDR(zzero), AB->u, M );
        flops = FLOPS_ZGEMM( M, B->rk, K );
    }

    PASTE_CORE_ZLRMM_VOID;
    return flops;
}

/**
 *******************************************************************************
 *
 * @brief Perform the operation AB = op(A) * op(B), with B full-rank and A and AB
 * low-rank.
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
 * @param[inout] AB
 *          The low-rank structure of the AB matrix in which to store the AB product.
 *
 * @param[inout] infomask
 *          The mask of informations returned by the core_zxx2lr() functions.
 *          - If AB.u is orthogonal on exit, then PASTIX_LRM3_ORTHOU is set.
 *          - If AB.u is allocated, then PASTIX_LRM3_ALLOCU is set.
 *          - If AB.v is allocated, then PASTIX_LRM3_ALLOCV is set.
 *          - If AB.v is inistialized as one of the given pointer and op(B) is not
 *          applyed, then PASTIX_LRM3_TRANSB is set.
 *
 * @param[in] Arkmin
 *          Threshold for which A->rk is considered as the final rank of AB
 *
 *******************************************************************************
 *
 * @return The number of flops required to perform the operation.
 *
 *******************************************************************************/
pastix_fixdbl_t
core_zlrfr2lr( core_zlrmm_t     *params,
               pastix_lrblock_t *AB,
               int              *infomask,
               pastix_int_t      Arkmin )
{
    PASTE_CORE_ZLRMM_PARAMS( params );
    pastix_int_t ldau, ldav, ldbu;
    pastix_fixdbl_t flops;

    ldau = (transA == PastixNoTrans) ? M : K;
    ldav = ( A->rk == -1 ) ? -1 : A->rkmax;
    ldbu = (transB == PastixNoTrans) ? K : N;

    /*
     *  A( M-by-ra x ra-by-K ) * B(N-by-K)^t
     */
    if ( A->rk > Arkmin ) {
        /*
         * We are in a similar case to the _Cfr function, and we
         * choose the optimal number of flops.
         */
        pastix_fixdbl_t flops1 = FLOPS_ZGEMM( A->rk, N, K ) + FLOPS_ZGEMM( M, N, A->rk );
        pastix_fixdbl_t flops2 = FLOPS_ZGEMM( M, K, A->rk ) + FLOPS_ZGEMM( M, N, K     );
        pastix_complex64_t *tmp;

        AB->rk    = -1;
        AB->rkmax = M;
        AB->v     = NULL;

        if ( flops1 <= flops2 ) {
            if ( (work = core_zlrmm_getws( params, A->rk * N + M * N )) == NULL ) {
                work = malloc( (A->rk * N + M * N) * sizeof(pastix_complex64_t) );
                *infomask |= PASTIX_LRM3_ALLOCU;
            }

            /* AB->u will be destroyed later */
            AB->u = work;
            tmp   = work + M * N;

            /*
             *  Au * (Av^t * B^t)
             */
            cblas_zgemm( CblasColMajor, CblasNoTrans, (CBLAS_TRANSPOSE)transB,
                         A->rk, N, K,
                         CBLAS_SADDR(zone),  A->v, ldav,
                                             B->u, ldbu,
                         CBLAS_SADDR(zzero), tmp,  A->rk );

            cblas_zgemm( CblasColMajor, CblasNoTrans, CblasNoTrans,
                         M, N, A->rk,
                         CBLAS_SADDR(zone),  A->u,  ldau,
                                             tmp,   A->rk,
                         CBLAS_SADDR(zzero), AB->u, M );

            flops = flops1;
        }
        else {
            if ( (work = core_zlrmm_getws( params, M * K + M * N )) == NULL ) {
                work = malloc( (M * K + M * N) * sizeof(pastix_complex64_t) );
                *infomask |= PASTIX_LRM3_ALLOCU;
            }

            /* AB->u will be destroyed later */
            AB->u = work;
            tmp   = work + M * N;

            /*
             *  (Au * Av^t) * B^t
             */
            cblas_zgemm( CblasColMajor, CblasNoTrans, CblasNoTrans,
                         M, K, A->rk,
                         CBLAS_SADDR(zone),  A->u, ldau,
                                             A->v, ldav,
                         CBLAS_SADDR(zzero), tmp,  M );

            cblas_zgemm( CblasColMajor, (CBLAS_TRANSPOSE)transA, (CBLAS_TRANSPOSE)transB,
                         M, N, K,
                         CBLAS_SADDR(zone),  tmp,   M,
                                             B->u,  ldbu,
                         CBLAS_SADDR(zzero), AB->u, M );

            flops = flops2;
        }
    }
    else {
        /*
         * A->rk is the smallest rank
         */
        AB->rk    = A->rk;
        AB->rkmax = A->rk;
        AB->u     = A->u;
        *infomask |= PASTIX_LRM3_ORTHOU;

        if ( (work = core_zlrmm_getws( params, A->rk * N )) == NULL ) {
            work = malloc( A->rk * N * sizeof(pastix_complex64_t) );
            *infomask |= PASTIX_LRM3_ALLOCV;
        }
        AB->v = work;

        cblas_zgemm( CblasColMajor, CblasNoTrans, (CBLAS_TRANSPOSE)transB,
                     A->rk, N, K,
                     CBLAS_SADDR(zone),  A->v,  ldav,
                                         B->u,  ldbu,
                     CBLAS_SADDR(zzero), AB->v, AB->rkmax );

        flops = FLOPS_ZGEMM( A->rk, N, K );
    }

    PASTE_CORE_ZLRMM_VOID;
    (void)infomask;
    return flops;
}

/**
 *******************************************************************************
 *
 * @brief Perform the operation AB = op(A) * op(B), with A, B, and AB low-rank.
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
 * @param[inout] AB
 *          The low-rank structure of the AB matrix in which to store the AB product.
 *
 * @param[inout] infomask
 *          The mask of informations returned by the core_zxx2lr() functions.
 *          - If AB.u is orthogonal on exit, then PASTIX_LRM3_ORTHOU is set.
 *          - If AB.u is allocated, then PASTIX_LRM3_ALLOCU is set.
 *          - If AB.v is allocated, then PASTIX_LRM3_ALLOCV is set.
 *          - If AB.v is inistialized as one of the given pointer and op(B) is not
 *          applyed, then PASTIX_LRM3_TRANSB is set.
 *
 *******************************************************************************
 *
 * @return The number of flops required to perform the operation.
 *
 *******************************************************************************/
pastix_fixdbl_t
core_zlrlr2lr( core_zlrmm_t     *params,
               pastix_lrblock_t *AB,
               int              *infomask )
{
    PASTE_CORE_ZLRMM_PARAMS( params );
    pastix_int_t ldau, ldav, ldbu, ldbv;
    pastix_complex64_t *work2;
    pastix_lrblock_t rArB;
    pastix_fixdbl_t flops = 0.0;
    int allocated = 0;

    assert( A->rk <= A->rkmax && A->rk > 0 );
    assert( B->rk <= B->rkmax && B->rk > 0 );
    assert( transA == PastixNoTrans );
    assert( transB != PastixNoTrans );

    *infomask = 0;
    ldau = (A->rk == -1) ? A->rkmax : M;
    ldav = A->rkmax;
    ldbu = (B->rk == -1) ? B->rkmax : N;
    ldbv = B->rkmax;

    if ( (work2 = core_zlrmm_getws( params, A->rk * B->rk )) == NULL ) {
        work2 = malloc( A->rk * B->rk * sizeof(pastix_complex64_t) );
        allocated = 1;
    }

    /*
     * Let's compute A * B' = Au Av^h (Bu Bv^h)' with the smallest ws
     */
    cblas_zgemm( CblasColMajor, CblasNoTrans, (CBLAS_TRANSPOSE)transB,
                 A->rk, B->rk, K,
                 CBLAS_SADDR(zone),  A->v,  ldav,
                                     B->v,  ldbv,
                 CBLAS_SADDR(zzero), work2, A->rk );
    flops = FLOPS_ZGEMM( A->rk, B->rk, K );

    /*
     * Try to compress (Av^h Bv^h')
     */
    flops += lowrank->core_ge2lr( lowrank->tolerance, -1, A->rk, B->rk, work2, A->rk, &rArB );

    /*
     * The rank of AB is not smaller than min(rankA, rankB)
     */
    if ( rArB.rk == -1 ) {
        if ( A->rk <= B->rk ) {
            /*
             *    ABu = Au
             *    ABv = (Av^h Bv^h') * Bu'
             */
            if ( (work = core_zlrmm_getws( params, A->rk * N )) == NULL ) {
                work = malloc( A->rk * N * sizeof(pastix_complex64_t) );
                *infomask |= PASTIX_LRM3_ALLOCV;
            }

            AB->rk    = A->rk;
            AB->rkmax = A->rk;
            AB->u     = A->u;
            AB->v     = work;
            *infomask |= PASTIX_LRM3_ORTHOU;

            cblas_zgemm( CblasColMajor, CblasNoTrans, (CBLAS_TRANSPOSE)transB,
                         A->rk, N, B->rk,
                         CBLAS_SADDR(zone),  work2, A->rk,
                                             B->u,  ldbu,
                         CBLAS_SADDR(zzero), AB->v, AB->rkmax );
            flops += FLOPS_ZGEMM( A->rk, N, B->rk );
        }
        else {
            /*
             *    ABu = Au * (Av^h Bv^h')
             *    ABv = Bu'
             */
            if ( (work = core_zlrmm_getws( params, B->rk * M )) == NULL ) {
                work = malloc( B->rk * M * sizeof(pastix_complex64_t) );
                *infomask |= PASTIX_LRM3_ALLOCU;
            }

            AB->rk = B->rk;
            AB->rkmax = B->rk;
            AB->u = work;
            AB->v = B->u;

            cblas_zgemm( CblasColMajor, CblasNoTrans, CblasNoTrans,
                         M, B->rk, A->rk,
                         CBLAS_SADDR(zone),  A->u,   ldau,
                                             work2,  A->rk,
                         CBLAS_SADDR(zzero), AB->u,  M );
            flops += FLOPS_ZGEMM( M, B->rk, A->rk );

            *infomask |= PASTIX_LRM3_TRANSB;
        }
    }
    else if ( rArB.rk == 0 ) {
        AB->rk    = 0;
        AB->rkmax = 0;
        AB->u = NULL;
        AB->v = NULL;
        *infomask |= PASTIX_LRM3_ORTHOU;
    }
    /**
     * The rank of AB is smaller than min(rankA, rankB)
     */
    else {
        if ( (work = core_zlrmm_getws( params, (M + N) * rArB.rk )) == NULL ) {
            work = malloc( (M + N) * rArB.rk * sizeof(pastix_complex64_t) );
            *infomask |= PASTIX_LRM3_ALLOCU;
        }

        AB->rk    = rArB.rk;
        AB->rkmax = rArB.rk;
        AB->u = work;
        AB->v = work + M * rArB.rk;
        *infomask |= PASTIX_LRM3_ORTHOU;

        cblas_zgemm( CblasColMajor, CblasNoTrans, CblasNoTrans,
                     M, rArB.rk, A->rk,
                     CBLAS_SADDR(zone),  A->u,   ldau,
                                         rArB.u, A->rk,
                     CBLAS_SADDR(zzero), AB->u,  M );

        cblas_zgemm( CblasColMajor, CblasNoTrans, (CBLAS_TRANSPOSE)transB,
                     rArB.rk, N, B->rk,
                     CBLAS_SADDR(zone),  rArB.v, rArB.rkmax,
                                         B->u,   ldbu,
                     CBLAS_SADDR(zzero), AB->v,  rArB.rk );

        flops += FLOPS_ZGEMM( M, rArB.rk, A->rk ) + FLOPS_ZGEMM( rArB.rk, N, B->rk );
    }
    core_zlrfree(&rArB);

    if ( allocated ) {
        free( work2 );
    }
    PASTE_CORE_ZLRMM_VOID;
    return flops;
}
