/**
 *
 * @file coeftab_z.c
 *
 * Precision dependent sequential routines to apply operation of the full matrix.
 *
 * @copyright 2015-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.1
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @author Gregoire Pichon
 * @author Mathieu Faverge
 * @date 2018-07-16
 *
 * @precisions normal z -> s d c
 *
 **/
#include "common.h"
#include "solver.h"
#include "lapacke.h"
#include "sopalin/coeftab_z.h"
#include "pastix_zcores.h"

/**
 *******************************************************************************
 *
 * @brief Dump the solver matrix coefficients into a file in human readable
 * format.
 *
 * All non-zeroes coefficients are dumped in the format:
 *    i j val
 * with one value per row.
 *
 *******************************************************************************
 *
 * @param[inout] pastix_data
 *          The pastix_data instance to access the unique directory id in which
 *          output the files.
 *
 * @param[in] solvmtx
 *          The solver matrix to print.
 *
 * @param[in] filename
 *          The filename where to store the output matrix.
 *
 *******************************************************************************/
void
coeftab_zdump( pastix_data_t      *pastix_data,
               const SolverMatrix *solvmtx,
               const char         *filename )
{
    SolverCblk *cblk = solvmtx->cblktab;
    pastix_int_t itercblk;
    FILE *stream = NULL;

    stream = pastix_fopenw( &(pastix_data->dirtemp), filename, "w" );
    if ( stream == NULL ){
        return;
    }

    /*
     * TODO: there is a problem right here for now, because there are no
     * distinctions between L and U coeffcients in the final file
     */
    for (itercblk=0; itercblk<solvmtx->cblknbr; itercblk++, cblk++)
    {
        cpucblk_zdump( PastixLCoef, cblk, stream );
        if ( NULL != cblk->ucoeftab )
            cpucblk_zdump( PastixUCoef, cblk, stream );
    }

    fclose( stream );
}

/**
 *******************************************************************************
 *
 * @brief Compare two solver matrices in full-rank format.
 *
 * The second solver matrix is overwritten by the difference of the two
 * matrices.  The frobenius norm of the difference of each column block is
 * computed and the functions returns 0 if the result for all the column blocks
 * of:
 *      || B_k - A_k || / ( || A_k || * eps )
 *
 * is below 10. Otherwise, an error message is printed and 1 is returned.
 *
 *******************************************************************************
 *
 * @param[in] side
 *          Define which side of the cblk must be tested.
 *          @arg PastixLCoef if lower part only
 *          @arg PastixUCoef if upper part only
 *          @arg PastixLUCoef if both sides.
 *
 * @param[in] solvA
 *          The solver matrix A.
 *
 * @param[inout] solvB
 *          The solver matrix B.
 *          On exit, B coefficient arrays are overwritten by the result of
 *          (B-A).
 *
 *******************************************************************************
 *
 * @return 0 if the test is passed, >= 0 otherwise.
 *
 *******************************************************************************/
int
coeftab_zdiff( pastix_coefside_t   side,
               const SolverMatrix *solvA,
               SolverMatrix       *solvB )
{
    SolverCblk *cblkA = solvA->cblktab;
    SolverCblk *cblkB = solvB->cblktab;
    pastix_int_t cblknum;
    int rc       = 0;
    int saved_rc = 0;

    for(cblknum=0; cblknum<solvA->cblknbr; cblknum++, cblkA++, cblkB++) {
        rc += cpucblk_zdiff( side, cblkA, cblkB );
        if ( rc != saved_rc ){
            fprintf(stderr, "CBLK %ld was not correctly compressed\n", (long)cblknum);
            saved_rc = rc;
        }
    }

    return rc;
}

/**
 *******************************************************************************
 *
 * @brief Compress all the cblks marked as valid for low-rank format.
 *
 * All the cblk in the top levels of the elimination tree markes as candidates
 * for compression are compressed if there is a gain to compress them. The
 * compression to low-rank format is parameterized by the input information
 * stored in the lowrank structure. On exit, all the cblks marked for
 * compression are stored through the low-rank structure, even if they are kept
 * in their full-rank form.
 *
 * @remark This routine is sequential
 *
 *******************************************************************************
 *
 * @param[inout] solvmtx
 *          The solver matrix of the problem to compress.
 *
 *******************************************************************************
 *
 * @return The memory gain resulting from the compression to low-rank format in
 * Bytes.
 *
 *******************************************************************************/
pastix_int_t
coeftab_zcompress( SolverMatrix *solvmtx )
{
    SolverCblk *cblk  = solvmtx->cblktab;
    pastix_coefside_t side = (solvmtx->factotype == PastixFactLU) ? PastixLUCoef : PastixLCoef;
    pastix_int_t cblknum, gain = 0;

    for(cblknum=0; cblknum<solvmtx->cblknbr; cblknum++, cblk++) {
        if ( cblk->cblktype & CBLK_COMPRESSED ) {
            gain += cpucblk_zcompress( side, cblk, solvmtx->lowrank );
        }
    }
    return gain;
}

/**
 *******************************************************************************
 *
 * @brief Uncompress all column block in low-rank format into full-rank format
 *
 *******************************************************************************
 *
 * @param[inout] solvmtx
 *          The solver matrix of the problem.
 *
 *******************************************************************************/
void
coeftab_zuncompress( SolverMatrix *solvmtx )
{
    SolverCblk  *cblk   = solvmtx->cblktab;
    pastix_int_t cblknum;
    pastix_coefside_t side = (solvmtx->factotype == PastixFactLU) ? PastixLUCoef : PastixLCoef;

    for(cblknum=0; cblknum<solvmtx->cblknbr; cblknum++, cblk++) {
        if (cblk->cblktype & CBLK_COMPRESSED) {
            cpucblk_zuncompress( side, cblk );
        }
    }
}

/**
 *******************************************************************************
 *
 * @brief Compute the memory gain of the low-rank form over the full-rank form
 * for the entire matrix.
 *
 * This function returns the memory gain in bytes for the full matrix when
 * column blocks are stored in low-rank format compared to a full rank storage.
 *
 *******************************************************************************
 *
 * @param[in] solvmtx
 *          The solver matrix of the problem.
 *
 *******************************************************************************/
void
coeftab_zmemory( const SolverMatrix *solvmtx )
{
    pastix_coefside_t side = (solvmtx->factotype == PastixFactLU) ? PastixLUCoef : PastixLCoef;
    SolverCblk  *cblk = solvmtx->cblktab;
    SolverBlok  *blok;
    pastix_int_t i, cblknum, in_height, off_height;
    pastix_int_t gain[4] = { 0, 0, 0, 0 };
    pastix_int_t orig[4] = { 0, 0, 0, 0 };
    pastix_fixdbl_t memgain[4];
    pastix_fixdbl_t memorig[4];
    pastix_fixdbl_t totgain, totorig;

    for(cblknum=0; cblknum<solvmtx->cblknbr; cblknum++, cblk++) {

        in_height = 0;
        blok = cblk->fblokptr;
        while( (blok < cblk[1].fblokptr) &&
               ((solvmtx->cblktab + blok->fcblknm)->sndeidx == cblk->sndeidx) )
        {
            in_height += blok_rownbr( blok );
            blok++;
        }
        off_height = cblk->stride - in_height;

        if ( !(cblk->cblktype & CBLK_COMPRESSED) ) {
            orig[FR_InDiag]  += cblk_colnbr( cblk ) * in_height;
            orig[FR_OffDiag] += cblk_colnbr( cblk ) * off_height;
        }
        else {
            orig[LR_InDiag]  += cblk_colnbr( cblk ) * in_height;
            orig[LR_OffDiag] += cblk_colnbr( cblk ) * off_height;
        }

        if (cblk->cblktype & CBLK_COMPRESSED) {
            cpucblk_zmemory( side, solvmtx, cblk, gain );
        }
    }

    if ( side == PastixLUCoef ) {
        orig[FR_InDiag]  *= 2;
        orig[FR_OffDiag] *= 2;
        orig[LR_InDiag]  *= 2;
        orig[LR_OffDiag] *= 2;
    }

    totgain = 0.;
    totorig = 0.;
    for (i=0; i<4; i++) {
        memgain[i] = (orig[i] - gain[i]) * pastix_size_of( PastixComplex64 );
        memorig[i] =  orig[i]            * pastix_size_of( PastixComplex64 );
        totgain += memgain[i];
        totorig += memorig[i];
    }

    pastix_print(0, 0,
                  "    Compression:\n"
                  "      ------------------------------------------------\n"
                  "      Full-rank cblk\n"
                  "        Inside supernodes                     %8.3g %co\n"
                  "        Outside supernodes                    %8.3g %co\n"
                  "      Low-rank cblk\n"
                  "        Inside supernodes       %8.3g %co / %8.3g %co\n"
                  "        Outside supernodes      %8.3g %co / %8.3g %co\n"
                  "      ------------------------------------------------\n"
                  "      Total                     %8.3g %co / %8.3g %co\n",
                  pastix_print_value(memorig[FR_InDiag] ), pastix_print_unit(memorig[FR_InDiag] ),
                  pastix_print_value(memorig[FR_OffDiag]), pastix_print_unit(memorig[FR_OffDiag]),
                  pastix_print_value(memgain[LR_InDiag] ), pastix_print_unit(memgain[LR_InDiag] ),
                  pastix_print_value(memorig[LR_InDiag] ), pastix_print_unit(memorig[LR_InDiag] ),
                  pastix_print_value(memgain[LR_OffDiag]), pastix_print_unit(memgain[LR_OffDiag]),
                  pastix_print_value(memorig[LR_OffDiag]), pastix_print_unit(memorig[LR_OffDiag]),
                  pastix_print_value(totgain),             pastix_print_unit(totgain),
                  pastix_print_value(totorig),             pastix_print_unit(totorig) );

    return;
}

/**
 *******************************************************************************
 *
 * @brief Extract the Schur complement
 *
 * This routine is sequential and returns the full Schur complement
 * uncommpressed in Lapack format.
 *
 *******************************************************************************
 *
 * @param[in] solvmtx
 *          The solver matrix structure describing the problem.
 *
 * @param[inout] S
 *          The pointer to the allocated matrix array that will store the Schur
 *          complement.
 *
 * @param[in] lds
 *          The leading dimension of the S array.
 *
 *******************************************************************************/
void
coeftab_zgetschur( const SolverMatrix *solvmtx,
                   pastix_complex64_t *S, pastix_int_t lds )
{
    SolverCblk *cblk = solvmtx->cblktab + solvmtx->cblkschur;
    pastix_complex64_t *localS;
    pastix_int_t itercblk, fcolnum, nbcol;
    int upper_part = (solvmtx->factotype == PastixFactLU);
    fcolnum = cblk->fcolnum;

    nbcol = solvmtx->nodenbr - fcolnum;
    assert( nbcol <= lds );

    /* Initialize the array to 0 */
    LAPACKE_zlaset_work( LAPACK_COL_MAJOR, 'A', nbcol, nbcol, 0., 0., S, lds );

    for (itercblk=solvmtx->cblkschur; itercblk<solvmtx->cblknbr; itercblk++, cblk++)
    {
        assert( cblk->cblktype & CBLK_IN_SCHUR );
        assert( lds >= cblk->stride );

        localS = S + (cblk->fcolnum - fcolnum) * lds + (cblk->fcolnum - fcolnum);

        cpucblk_zgetschur( cblk, upper_part, localS, lds );
    }
}

/**
 *******************************************************************************
 *
 * @brief Extract the diagonal
 *
 * This routine is sequential and returns the full diagonal in the vector D,
 * such that:
 *     D[incD*i]= A(i, i)
 *
 *******************************************************************************
 *
 * @param[in] solvmtx
 *          The solver matrix structure describing the problem.
 *
 * @param[inout] D
 *          The pointer to the allocated vector array that will store the diagonal.
 *          D must be of size solvmtx->nodenbr * incD.
 *
 * @param[in] incD
 *          The increment bewteen two elements of D. incD > 0.
 *
 *******************************************************************************/
void
coeftab_zgetdiag( const SolverMatrix *solvmtx,
                  pastix_complex64_t *D, pastix_int_t incD )
{
    SolverCblk *cblk = solvmtx->cblktab;
    pastix_complex64_t *A;
    pastix_int_t lda, itercblk, nbcol, i;

    for (itercblk=0; itercblk<solvmtx->cblknbr; itercblk++, cblk++)
    {
        nbcol = cblk_colnbr( cblk );
        if ( cblk->cblktype & CBLK_COMPRESSED ) {
            assert( cblk->fblokptr->LRblock[0].rk == -1 );
            A   = cblk->fblokptr->LRblock[0].u;
            lda = cblk_colnbr( cblk ) + 1;
        }
        else {
            A = cblk->lcoeftab;

            if ( cblk->cblktype & CBLK_LAYOUT_2D ) {
                lda = cblk_colnbr( cblk ) + 1;
            }
            else {
                lda = cblk->stride + 1;
            }
        }

        for (i=0; i<nbcol; i++, D += incD, A += lda ) {
            *D = *A;
        }
    }
}
