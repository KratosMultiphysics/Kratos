/**
 *
 * @file cpucblk_zcompress.c
 *
 * Precision dependent function to compress/uncompress the coefficients
 *
 * @copyright 2015-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.1
 * @author Gregoire Pichon
 * @author Mathieu Faverge
 * @date 2018-07-16
 *
 * @precisions normal z -> s d c
 *
 **/
#include "common.h"
#include "solver.h"
#include <lapacke.h>
#include "kernels_trace.h"
#include "pastix_zcores.h"
#include "pastix_zlrcores.h"

/**
 *******************************************************************************
 *
 * @brief Compress a single column block from full-rank to low-rank format
 *
 * The compression to low-rank format is parameterized by the input information
 * stored in the lowrank structure.
 *
 *******************************************************************************
 *
 * @param[in] side
 *          Define which side of the matrix must be initialized.
 *          @arg PastixLCoef if lower part only
 *          @arg PastixUCoef if upper part only
 *          @arg PastixLUCoef if both sides.
 *
 * @param[inout] cblk
 *          The column block to compress.
 *
 * @param[in] lowrank
 *          The low-rank parameters to describe the low-rank compression format.
 *
 *******************************************************************************
 *
 * @return The memory gain resulting from the compression to low-rank format in
 *         number of elements.
 *
 *******************************************************************************/
pastix_int_t
cpucblk_zcompress( pastix_coefside_t side,
                   SolverCblk       *cblk,
                   pastix_lr_t       lowrank )
{
    pastix_lrblock_t   *lrA;
    pastix_complex64_t *A;
    SolverBlok         *blok  = cblk[0].fblokptr + 1;
    SolverBlok         *lblok = cblk[1].fblokptr;
    pastix_int_t        ncols = cblk_colnbr( cblk );
    pastix_int_t        gain;
    pastix_int_t        gainL = 0;
    pastix_int_t        gainU = 0;
    pastix_fixdbl_t     flops;

    assert( cblk->cblktype & CBLK_LAYOUT_2D  );
    assert( cblk->cblktype & CBLK_COMPRESSED );
    assert( ncols > lowrank.compress_min_width );

    for (; blok<lblok; blok++)
    {
        pastix_int_t nrows = blok_rownbr( blok );

        if ( nrows > lowrank.compress_min_height )
        {
            gain = nrows * ncols;

            /* Lower part */
            if ( side != PastixUCoef ) {
                lrA = blok->LRblock;

                /* Try to compress */
                if ( lrA->rk == -1 ) {
                    A = lrA->u;

                    kernel_trace_start_lvl2( PastixKernelLvl2_LR_init_compress );
                    flops = lowrank.core_ge2lr( lowrank.tolerance, -1, nrows, ncols,
                                                A, nrows, lrA );
                    kernel_trace_stop_lvl2_rank( flops, lrA->rk );

                    free( A );
                }

                if  ( lrA->rk != -1 ) {
                    gainL += gain - ((nrows+ncols) * lrA->rk);
                }
            }

            /* Upper part */
            if ( side != PastixLCoef ) {
                lrA = blok->LRblock + 1;

                if( lrA->rk == -1 ) {
                    A = lrA->u;

                    kernel_trace_start_lvl2( PastixKernelLvl2_LR_init_compress );
                    flops = lowrank.core_ge2lr( lowrank.tolerance, -1, nrows, ncols,
                                                A, nrows, lrA );
                    kernel_trace_stop_lvl2_rank( flops, lrA->rk );

                    free( A );
                }

                if  ( lrA->rk != -1 ) {
                    gainU += gain - ((nrows+ncols) * lrA->rk);
                }
            }
        }
    }

    return gainL + gainU;
}

/**
 *******************************************************************************
 *
 * @brief Uncompress a single column block from low-rank format to full-rank
 * format.
 *
 *******************************************************************************
 *
 * @param[in] side
 *          Define which side of the matrix must be initialized.
 *          @arg PastixLCoef if lower part only
 *          @arg PastixUCoef if upper part only
 *          @arg PastixLUCoef if both sides.
 *
 * @param[inout] cblk
 *          The column block to uncompress.
 *
 *******************************************************************************/
void
cpucblk_zuncompress( pastix_coefside_t side,
                     SolverCblk       *cblk )
{
    SolverBlok  *blok, *lblok;
    pastix_int_t ncols = cblk_colnbr( cblk );
    int ret;

    if ( side != PastixUCoef ) {
        blok  = cblk[0].fblokptr;
        lblok = cblk[1].fblokptr;
        for (; blok<lblok; blok++)
        {
            pastix_lrblock_t lrtmp;
            pastix_int_t nrows = blok_rownbr( blok );

            memcpy( &lrtmp, blok->LRblock, sizeof(pastix_lrblock_t) );

            core_zlralloc( nrows, ncols, -1, blok->LRblock );
            ret = core_zlr2ge( PastixNoTrans, nrows, ncols,
                               &lrtmp,
                               blok->LRblock[0].u, nrows );
            assert( ret == 0 );
            core_zlrfree( &lrtmp );
        }
    }

    if ( side != PastixLCoef ) {
        blok  = cblk[0].fblokptr;
        lblok = cblk[1].fblokptr;
        for (; blok<lblok; blok++)
        {
            pastix_lrblock_t lrtmp;
            pastix_int_t nrows = blok_rownbr( blok );

            memcpy( &lrtmp, blok->LRblock+1, sizeof(pastix_lrblock_t) );

            core_zlralloc( nrows, ncols, -1, blok->LRblock+1 );
            ret = core_zlr2ge( PastixNoTrans, nrows, ncols,
                               &lrtmp,
                               blok->LRblock[1].u, nrows );
            assert( ret == 0 );
            core_zlrfree( &lrtmp );
        }
    }

    (void)ret;
}

/**
 *******************************************************************************
 *
 * @brief Return the memory gain of the low-rank form over the full-rank form
 * for a single column-block.
 *
 * This function returns the memory gain in number of elements for a single
 * column block when it is stored in low-rank format compared to a full rank
 * storage.
 *
 *******************************************************************************
 *
 * @param[in] side
 *          Define which side of the matrix must be initialized.
 *          @arg PastixLCoef if lower part only
 *          @arg PastixUCoef if upper part only
 *          @arg PastixLUCoef if both sides.
 *
 * @param[in] solvmtx
 *          The pointer to the solver structure.
 *
 * @param[in] cblk
 *          The column block to study.
 *
 * @param[in] gain
 *          The structure that counts gain on each type of block.
 *
 *******************************************************************************/
void
cpucblk_zmemory( pastix_coefside_t   side,
                 const SolverMatrix *solvmtx,
                 const SolverCblk   *cblk,
                 pastix_int_t       *gain )
{
    SolverBlok *blok  = cblk[0].fblokptr + 1;
    SolverBlok *lblok = cblk[1].fblokptr;

    pastix_int_t ncols = cblk_colnbr( cblk );
    pastix_int_t size;
    pastix_int_t gainblok;

    for (; blok<lblok; blok++)
    {
        pastix_int_t nrows = blok_rownbr( blok );
        size = nrows * ncols;
        gainblok = 0;

        /* Lower part */
        if ( (side != PastixUCoef) &&
             (blok->LRblock[0].rk >= 0) )
        {
            gainblok += (size - ((nrows+ncols) * blok->LRblock[0].rkmax));
        }

        if ( (side != PastixLCoef) &&
             (blok->LRblock[1].rk >= 0) )
        {
            gainblok += (size - ((nrows+ncols) * blok->LRblock[1].rkmax));
        }

        if ( (solvmtx->cblktab + blok->fcblknm)->sndeidx == cblk->sndeidx ) {
            gain[LR_InDiag] += gainblok;
        }
        else {
            gain[LR_OffDiag] += gainblok;
        }
    }

    return;
}
