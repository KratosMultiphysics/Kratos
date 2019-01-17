/**
 *
 * @file coeftab_zinit.c
 *
 * Precision dependent coeficient array initialization routines.
 *
 * @copyright 2015-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.1
 * @author Xavier Lacoste
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @date 2018-07-16
 *
 * @precisions normal z -> s d c
 *
 **/
#define _GNU_SOURCE 1
#include "common.h"
#include "solver.h"
#include "bcsc.h"
#include "sopalin/coeftab_z.h"
#include "pastix_zcores.h"

/**
 *******************************************************************************
 *
 * @brief Dump a single column block into a FILE in a human readale format.
 *
 * All non-zeroes coefficients are dumped in the format:
 *    i j val
 * with one value per row.
 *
 *******************************************************************************
 *
 * @param[in] side
 *          Define which side of the cblk must be printed.
 *          @arg PastixLCoef if lower part only
 *          @arg PastixUCoef if upper part only
 *
 * @param[in] cblk
 *          The column block to dump into the file.
 *
 * @param[inout] stream
 *          The FILE structure opened in write mode.
 *
 *******************************************************************************/
void
cpucblk_zdump( pastix_coefside_t side,
               const SolverCblk *cblk,
               FILE             *stream )
{
    const pastix_complex64_t *coeftab = side == PastixUCoef ? cblk->ucoeftab : cblk->lcoeftab;
    SolverBlok  *blok;
    pastix_int_t itercol;
    pastix_int_t iterrow;
    pastix_int_t coefindx;

    /* We don't know how to dump the compressed block for now */
    if ( cblk->cblktype & CBLK_COMPRESSED ) {
        fprintf(stderr, "coeftab_zcblkdump: Can't dump a compressed cblk\n");
        return;
    }

    for (itercol  = cblk->fcolnum;
         itercol <= cblk->lcolnum;
         itercol++)
    {
        /* Diagonal Block */
        blok     = cblk->fblokptr;
        coefindx = blok->coefind;
        if (cblk->cblktype & CBLK_LAYOUT_2D) {
            coefindx += (itercol - cblk->fcolnum) * blok_rownbr( blok );
        }
        else {
            coefindx += (itercol - cblk->fcolnum) * cblk->stride;
        }

        for (iterrow  = blok->frownum;
             iterrow <= blok->lrownum;
             iterrow++, coefindx++)
        {
            if ((cabs( coeftab[coefindx] ) > 0.) &&
                (itercol <= iterrow))
            {
                if ( side == PastixUCoef ) {
#if defined(PRECISION_z) || defined(PRECISION_c)
                    fprintf(stream, "%ld %ld (%13e,%13e) [U]\n",
                            (long)itercol, (long)iterrow,
                            creal(coeftab[coefindx]), cimag(coeftab[coefindx]));
#else
                    fprintf(stream, "%ld %ld %13e [U]\n",
                            (long)itercol, (long)iterrow,
                            coeftab[coefindx]);
#endif
                }
                else {
#if defined(PRECISION_z) || defined(PRECISION_c)
                    fprintf(stream, "%ld %ld (%13e,%13e) [L]\n",
                            (long)iterrow, (long)itercol,
                            creal(coeftab[coefindx]), cimag(coeftab[coefindx]));
#else
                    fprintf(stream, "%ld %ld %13e [L]\n",
                            (long)iterrow, (long)itercol,
                            coeftab[coefindx]);
#endif
                }
            }
        }

        /* Off diagonal blocks */
        blok++;
        while( blok < (cblk+1)->fblokptr )
        {
            coefindx  = blok->coefind;
            if (cblk->cblktype & CBLK_LAYOUT_2D) {
                coefindx += (itercol - cblk->fcolnum) * blok_rownbr( blok );
            }
            else {
                coefindx += (itercol - cblk->fcolnum) * cblk->stride;
            }

            for (iterrow  = blok->frownum;
                 iterrow <= blok->lrownum;
                 iterrow++, coefindx++)
            {
                if (cabs( coeftab[coefindx]) > 0.)
                {
                    if ( side == PastixUCoef ) {
#if defined(PRECISION_z) || defined(PRECISION_c)
                        fprintf(stream, "%ld %ld (%13e,%13e) [U]\n",
                                (long)itercol, (long)iterrow,
                                creal(coeftab[coefindx]), cimag(coeftab[coefindx]));
#else
                        fprintf(stream, "%ld %ld %13e [U]\n",
                                (long)itercol, (long)iterrow,
                                coeftab[coefindx]);
#endif
                    }
                    else {
#if defined(PRECISION_z) || defined(PRECISION_c)
                        fprintf(stream, "%ld %ld (%13e,%13e) [L]\n",
                                (long)iterrow, (long)itercol,
                                creal(coeftab[coefindx]), cimag(coeftab[coefindx]));
#else
                        fprintf(stream, "%ld %ld %13e [L]\n",
                                (long)iterrow, (long)itercol,
                                coeftab[coefindx]);
#endif
                    }
                }
            }
            blok++;
        }
    }
}

/**
 *******************************************************************************
 *
 * @brief Fully initialize a single cblk.
 *
 * The cblk is allocated, intialized from the bcsc, and compressed if necessary.
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
 *          The solver matrix data structure.
 *
 * @param[in] bcsc
 *          The internal block CSC structure to fill-in the matrix.
 *
 * @param[in] itercblk
 *          The index of the cblk to initialize.
 *
 * @param[inout] directory
 *          The pointer to the temporary directory where to store the output
 *          files.  Used only if PASTIX_DEBUG_DUMP_COEFTAB is defined.
 *
 *******************************************************************************/
void
cpucblk_zinit( pastix_coefside_t    side,
               const SolverMatrix  *solvmtx,
               const pastix_bcsc_t *bcsc,
               pastix_int_t         itercblk,
               char               **directory )
{
    pastix_int_t compress_when = solvmtx->lowrank.compress_when;
    SolverCblk  *cblk = solvmtx->cblktab + itercblk;

    cpucblk_zalloc( side, cblk );
    cpucblk_zfillin( side, solvmtx, bcsc, itercblk );

#if defined(PASTIX_DEBUG_DUMP_COEFTAB)
    /*
     * Rk: This function is not in the kernel directory to avoid the double
     * dependency with the pastix library due to pastix_fopenw()
     */
    {
        FILE *f = NULL;
        char *filename;
        int rc;

        /* Lower part */
        if ( side != PastixUCoef )
        {
            rc = asprintf( &filename, "Lcblk%05ld_init.txt", itercblk );
            f  = pastix_fopenw( directory, filename, "w" );
            if ( f != NULL ) {
                cpucblk_zdump( PastixLCoef, cblk, f );
                fclose( f );
            }
            free( filename );
        }

        /* Upper part */
        if ( side != PastixLCoef )
        {
            rc = asprintf( &filename, "Ucblk%05ld_init.txt", itercblk );
            f  = pastix_fopenw( directory, filename, "w" );
            if ( f != NULL ) {
                cpucblk_zdump( PastixUCoef, cblk, f );
                fclose( f );
            }
            free( filename );
        }
        (void)rc;
    }
#endif /* defined(PASTIX_DEBUG_DUMP_COEFTAB) */

    /**
     * Try to compress the cblk if needs to be compressed
     * TODO: change the criteria based on the level in the tree
     */
    if ( (cblk->cblktype & CBLK_COMPRESSED)                          &&
         (compress_when == PastixCompressWhenBegin)                  &&
         (cblk_colnbr( cblk ) > solvmtx->lowrank.compress_min_width) )
    {
        cpucblk_zcompress( side, cblk, solvmtx->lowrank );
    }

    (void)directory;
}
