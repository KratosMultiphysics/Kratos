/**
 *
 * @file core_zlrmm.c
 *
 * PaStiX low-rank kernel routines to compute a matrix matrix product in either
 * form low or full rank.
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
#include "kernels_trace.h"
#include "pastix_zlrcores.h"

/**
 *******************************************************************************
 *
 * @brief Compute the matrix matrix product when applied to a full rank matrix.
 *
 * This function considers that the C matrix is full rank, and A and B are
 * either full-rank or low-rank.  The result of the product is directly applied
 * to the C matrix.
 *
 *******************************************************************************
 *
 * @param[inout] params
 *          The LRMM structure that stores all the parameters used in the LRMM
 *          functions family.
 *          @sa core_zlrmm_t
 *
 *******************************************************************************
 *
 * @return The number of flops required to perform the operation.
 *
 *******************************************************************************/
static inline pastix_fixdbl_t
core_zlrmm_Cfr( core_zlrmm_t *params )
{
    const pastix_lrblock_t *A = params->A;
    const pastix_lrblock_t *B = params->B;
    pastix_fixdbl_t flops = 0.0;

    assert( params->transA == PastixNoTrans );
    assert( params->transB != PastixNoTrans );
    assert( A->rk <= A->rkmax && A->rk != 0 );
    assert( B->rk <= B->rkmax && B->rk != 0 );
    assert( params->C->rk == -1 );

    if ( A->rk == -1 ) {
        if ( B->rk == -1 ) {
            kernel_trace_start_lvl2( PastixKernelLvl2_LR_FRFR2FR );
            flops = core_zfrfr2fr( params );
            kernel_trace_stop_lvl2( flops );
        }
        else {
            kernel_trace_start_lvl2( PastixKernelLvl2_LR_FRLR2FR );
            flops = core_zfrlr2fr( params );
            kernel_trace_stop_lvl2( flops );
        }
    }
    else {
        if ( B->rk == -1 ) {
            kernel_trace_start_lvl2( PastixKernelLvl2_LR_LRFR2FR );
            flops = core_zlrfr2fr( params );
            kernel_trace_stop_lvl2( flops );
        }
        else {
            kernel_trace_start_lvl2( PastixKernelLvl2_LR_LRLR2FR );
            flops = core_zlrlr2fr( params );
            kernel_trace_stop_lvl2( flops );
        }
    }

    assert( params->C->rk == -1 );

    return flops;
}

/**
 *******************************************************************************
 *
 * @brief Compute the matrix matrix product when applied to a null matrix.
 *
 * This function considers that the C matrix is null, and A and B are either
 * full-rank or low-rank.  The result of the product is directly applied to the
 * C matrix.
 *
 *******************************************************************************
 *
 * @param[inout] params
 *          The LRMM structure that stores all the parameters used in the LRMM
 *          functions family.
 *          @sa core_zlrmm_t
 *
 *******************************************************************************
 *
 * @return The number of flops required to perform the operation.
 *
 *******************************************************************************/
static inline pastix_fixdbl_t
core_zlrmm_Cnull( core_zlrmm_t *params )
{
    PASTE_CORE_ZLRMM_PARAMS( params );
    pastix_lrblock_t AB;
    pastix_trans_t transV = PastixNoTrans;
    int infomask = 0;
    pastix_fixdbl_t flops = 0.0;

    assert(transA == PastixNoTrans);
    assert(transB != PastixNoTrans);
    assert( A->rk <= A->rkmax && A->rk != 0 );
    assert( B->rk <= B->rkmax && B->rk != 0 );

    if ( A->rk == -1 ) {
        if ( B->rk == -1 ) {
            kernel_trace_start_lvl2( PastixKernelLvl2_LR_FRFR2null );
            flops = core_zfrfr2lr( params, &AB, &infomask,
                                   pastix_imin( pastix_imin( M, N ),
                                                core_get_rklimit( Cm, Cn ) ) );
            kernel_trace_stop_lvl2( flops );
        }
        else {
            kernel_trace_start_lvl2( PastixKernelLvl2_LR_FRLR2null );
            flops = core_zfrlr2lr( params, &AB, &infomask,
                                   pastix_imin( M, core_get_rklimit( Cm, Cn ) ) );
            kernel_trace_stop_lvl2( flops );
        }
    }
    else {
        if ( B->rk == -1 ) {
            kernel_trace_start_lvl2( PastixKernelLvl2_LR_LRFR2null );
            flops = core_zlrfr2lr( params, &AB, &infomask,
                                   pastix_imin( N, core_get_rklimit( Cm, Cn ) ) );
            kernel_trace_stop_lvl2( flops );
        }
        else {
            kernel_trace_start_lvl2( PastixKernelLvl2_LR_LRLR2null );
            flops = core_zlrlr2lr( params, &AB, &infomask );
            kernel_trace_stop_lvl2( flops );

            assert( AB.rk != -1 );
            assert( AB.rkmax != -1 );
        }
    }

    if ( infomask & PASTIX_LRM3_TRANSB ) {
        transV = transB;
    }

    if ( AB.rk != 0 ) {
        pastix_atomic_lock( lock );
        switch ( C->rk ) {
        case -1:
            /*
             * C became full rank
             */
            flops += core_zlr2fr( params, &AB, transV );
            break;

        case 0:
            /*
             * C is still null
             */
            flops += core_zlr2null( params, &AB, transV, infomask );
            break;

        default:
            /*
             * C is low-rank of rank k
             */
            flops += core_zlr2lr( params, &AB, transV );
        }
        pastix_atomic_unlock( lock );
    }

    /* Free memory from zlrm3 */
    if ( infomask & PASTIX_LRM3_ALLOCU ) {
        free(AB.u);
    }
    if ( infomask & PASTIX_LRM3_ALLOCV ) {
        free(AB.v);
    }

    assert( C->rk <= C->rkmax);
    PASTE_CORE_ZLRMM_VOID;

    return flops;
}

/**
 *******************************************************************************
 *
 * @brief Compute the matrix matrix product when applied to a low rank matrix.
 *
 * This function considers that the C matrix is low rank, and A and B are
 * either full-rank or low-rank.  The result of the product is directly applied
 * to the C matrix.
 *
 *******************************************************************************
 *
 * @param[inout] params
 *          The LRMM structure that stores all the parameters used in the LRMM
 *          functions family.
 *          @sa core_zlrmm_t
 *
 *******************************************************************************
 *
 * @return The number of flops required to perform the operation.
 *
 *******************************************************************************/
static inline pastix_fixdbl_t
core_zlrmm_Clr( core_zlrmm_t *params )
{
    PASTE_CORE_ZLRMM_PARAMS( params );
    pastix_lrblock_t AB;
    pastix_trans_t transV = PastixNoTrans;
    int infomask = 0;
    pastix_fixdbl_t flops = 0.0;

    assert(transA == PastixNoTrans);
    assert(transB != PastixNoTrans);
    assert( A->rk <= A->rkmax && A->rk != 0 );
    assert( B->rk <= B->rkmax && B->rk != 0 );

    if ( A->rk == -1 ) {
        if ( B->rk == -1 ) {
            kernel_trace_start_lvl2( PastixKernelLvl2_LR_FRFR2LR );
            flops = core_zfrfr2lr( params, &AB, &infomask,
                                   pastix_imin( M, N ) );
            kernel_trace_stop_lvl2( flops );
        }
        else {
            kernel_trace_start_lvl2( PastixKernelLvl2_LR_FRLR2LR );
            flops = core_zfrlr2lr( params, &AB, &infomask, M );
            kernel_trace_stop_lvl2( flops );
        }
    }
    else {
        if ( B->rk == -1 ) {
            kernel_trace_start_lvl2( PastixKernelLvl2_LR_LRFR2LR );
            flops = core_zlrfr2lr( params, &AB, &infomask, N );
            kernel_trace_stop_lvl2( flops );
        }
        else {
            kernel_trace_start_lvl2( PastixKernelLvl2_LR_LRLR2LR );
            flops = core_zlrlr2lr( params, &AB, &infomask );
            kernel_trace_stop_lvl2( flops );

            assert( AB.rk != -1 );
            assert( AB.rkmax != -1 );
        }
    }

    if ( infomask & PASTIX_LRM3_TRANSB ) {
        transV = transB;
    }

    if ( AB.rk != 0 ) {
        pastix_atomic_lock( lock );
        switch ( C->rk ) {
        case -1:
            /*
             * C became full rank
             */
            flops += core_zlr2fr( params, &AB, transV );
            break;

        case 0:
            /*
             * C is still null
             */
            flops += core_zlr2null( params, &AB, transV, infomask );
            break;

        default:
            /*
             * C is low-rank of rank k
             */
            flops += core_zlr2lr( params, &AB, transV );
        }
        pastix_atomic_unlock( lock );
    }

    /* Free memory from zlrm3 */
    if ( infomask & PASTIX_LRM3_ALLOCU ) {
        free(AB.u);
    }
    if ( infomask & PASTIX_LRM3_ALLOCV ) {
        free(AB.v);
    }

    assert( C->rk <= C->rkmax);
    PASTE_CORE_ZLRMM_VOID;

    return flops;
}

/**
 *******************************************************************************
 *
 * @brief Compute the matrix matrix product when involved matrices are stored in
 * a low-rank structure.
 *
 * This function considers the generic matrix matrix product added to a third
 * matric C. All matrices are either null, low-rank or full-rank.
 *
 *******************************************************************************
 *
 * @param[inout] params
 *          The LRMM structure that stores all the parameters used in the LRMM
 *          functions family.
 *          @sa core_zlrmm_t
 *
 *******************************************************************************
 *
 * @return The number of flops required to perform the operation.
 *
 *******************************************************************************/
pastix_fixdbl_t
core_zlrmm( core_zlrmm_t *params )
{
    PASTE_CORE_ZLRMM_PARAMS( params );
    pastix_fixdbl_t flops;

    assert( transA == PastixNoTrans );
    assert( transB != PastixNoTrans );
    assert( A->rk <= A->rkmax);
    assert( B->rk <= B->rkmax);
    assert( C->rk <= C->rkmax);

    /* Quick return if multiplication by 0 */
    if ( A->rk == 0 || B->rk == 0 ) {
        return 0.0;
    }

    params->lwused = 0;

    /* TODO: this is a temporary fix */
    if ( lwork == 0 ) {
        params->work = work = NULL;
    }
    assert( ((work != NULL) && (lwork >  0)) ||
            ((work == NULL) && (lwork <= 0)) );

    if ( C->rk == 0 ) {
#if 0
        /* Debug case with NULL matrix replaced by the unit vector times 0 */
        pastix_cblk_lock( fcblk );
        if ( C->rk == 0 ) {
            core_zlrsze( 0, Cm, Cn, C, 1, 1, -1 );
            memset(C->u, 0, Cm * sizeof(pastix_complex64_t) );
            memset(C->v, 0, Cn * sizeof(pastix_complex64_t) );
            ((pastix_complex64_t*)(C->u))[0] = 1.;
        }
        pastix_cblk_unlock( fcblk );
        flops = core_zlrmm_Clr( params );
#else
        flops = core_zlrmm_Cnull( params );
#endif
    }
    else if ( C->rk == -1 ) {
        flops = core_zlrmm_Cfr( params );
    }
    else {
        flops = core_zlrmm_Clr( params );
    }

#if defined(PASTIX_DEBUG_LR)
    pastix_atomic_lock( lock );
    if ( (C->rk > 0) && (lowrank->compress_method == PastixCompressMethodRRQR) ) {
        int rc = core_zlrdbg_check_orthogonality( Cm, C->rk, (pastix_complex64_t*)C->u, Cm );
        if (rc == 1) {
            fprintf(stderr, "Failed to have u orthogonal in exit of lrmm\n" );
        }
    }
    pastix_atomic_unlock( lock );
#endif

    PASTE_CORE_ZLRMM_VOID;
    return flops;
}
