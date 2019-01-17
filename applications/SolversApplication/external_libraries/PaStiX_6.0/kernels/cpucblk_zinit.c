/**
 *
 * @file cpucblk_zinit.c
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
#include "common.h"
#include "solver.h"
#include "bcsc.h"
#include "pastix_zcores.h"
#include "pastix_zlrcores.h"

/**
 *******************************************************************************
 *
 * @brief Allocate the cblk structure to store the coefficient
 *
 * When stored in low-rank format, the data pointer in the low-rank structure of
 * each block must be initialized.
 * This routines performs only the allocation and is thread-safe if called in
 * parallel on the Lower and upper part.
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
 *          The column block to allocate.
 *
 *******************************************************************************/
void
cpucblk_zalloc( pastix_coefside_t  side,
                SolverCblk        *cblk )
{
    pastix_int_t ncols = cblk_colnbr( cblk );

    if ( cblk->cblktype & CBLK_COMPRESSED ) {
        pastix_lrblock_t   *LRblocks = NULL;
        SolverBlok         *blok     = cblk[0].fblokptr;
        SolverBlok         *lblok    = cblk[1].fblokptr;

        /* H then split */
        assert( cblk->cblktype & CBLK_LAYOUT_2D );

        LRblocks = blok->LRblock;

        if ( LRblocks == NULL ) {
            /* One allocation per cblk */
            LRblocks = malloc( 2 * (lblok - blok) * sizeof(pastix_lrblock_t) );
            memset( LRblocks, 0, 2 * (lblok - blok) * sizeof(pastix_lrblock_t) );
            if (!pastix_atomic_cas_xxb( &(blok->LRblock), (uint64_t)NULL, (uint64_t)LRblocks, sizeof(void*) )) {
                free( LRblocks );
                LRblocks = blok->LRblock;
            }
        }
        assert( LRblocks != NULL );

        for (; blok<lblok; blok++)
        {
            pastix_int_t nrows = blok_rownbr( blok );
            blok->LRblock = LRblocks;

            if ( side != PastixUCoef ) {
                core_zlralloc( nrows, ncols, -1, LRblocks );
            }
            LRblocks++;

            if ( side != PastixLCoef ) {
                core_zlralloc( nrows, ncols, -1, LRblocks );
            }
            LRblocks++;
        }
    }
    else {
        size_t coefnbr = cblk->stride * ncols;

        if ( side != PastixUCoef ) {
            assert( cblk->lcoeftab == NULL );

            MALLOC_INTERN( cblk->lcoeftab, coefnbr, pastix_complex64_t );
            memset( cblk->lcoeftab, 0, coefnbr * sizeof(pastix_complex64_t) );
        }

        if ( side != PastixLCoef ) {
            assert( cblk->ucoeftab == NULL );

            MALLOC_INTERN( cblk->ucoeftab, coefnbr, pastix_complex64_t );
            memset( cblk->ucoeftab, 0, coefnbr * sizeof(pastix_complex64_t) );
        }
    }
}

/**
 *******************************************************************************
 *
 * @brief Initialize the full-rank coeftab structure from the internat bcsc.
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
 *          PaStiX structure to store numerical data and flags
 *
 * @param[in] bcsc
 *          The internal bcsc structure that hold the graph with permutation
 *          stored by cblk.
 *
 * @param[in] itercblk
 *          The index of the cblk to fill in both bcsc and solvmtx structures.
 *
 *******************************************************************************/
static inline void
cpucblk_zfillin_fr( pastix_coefside_t    side,
                    const SolverMatrix  *solvmtx,
                    const pastix_bcsc_t *bcsc,
                    pastix_int_t         itercblk )
{
    const bcsc_cblk_t *csccblk = bcsc->cscftab + itercblk;
    SolverCblk *solvcblk = solvmtx->cblktab + itercblk;
    SolverBlok *solvblok;
    SolverBlok *lsolvblok = (solvcblk+1)->fblokptr;
    pastix_complex64_t *lcoeftab = solvcblk->lcoeftab;
    pastix_complex64_t *ucoeftab = solvcblk->ucoeftab;
    pastix_complex64_t *Lvalues  = bcsc->Lvalues;
    pastix_complex64_t *Uvalues  = bcsc->Uvalues;
    pastix_int_t ldd = solvcblk->stride;
    pastix_int_t itercoltab, iterval, coefindx;
    int is2d = solvcblk->cblktype & CBLK_LAYOUT_2D;

    assert( (side != PastixUCoef) || (ucoeftab != NULL) );

    for (itercoltab=0; itercoltab<csccblk->colnbr; itercoltab++)
    {
        pastix_int_t frow = csccblk->coltab[itercoltab];
        pastix_int_t lrow = csccblk->coltab[itercoltab+1];
        solvblok = solvcblk->fblokptr;
        if ( is2d ) {
            ldd = blok_rownbr( solvblok );
        }

        for (iterval=frow; iterval<lrow; iterval++)
        {
            pastix_int_t rownum = bcsc->rowtab[iterval];

            /* If values in the lower part of the matrix */
            if (rownum >= (solvcblk->fcolnum+itercoltab))
            {
                while ((solvblok < lsolvblok) &&
                       ((solvblok->lrownum < rownum) ||
                        (solvblok->frownum > rownum)))
                {
                    solvblok++;
                    if ( is2d ) {
                        ldd = blok_rownbr( solvblok );
                    }
                }

                if ( solvblok < lsolvblok )
                {
                    coefindx  = solvblok->coefind;
                    coefindx += rownum - solvblok->frownum; /* Row shift    */
                    coefindx += itercoltab * ldd;           /* Column shift */

                    if ( side != PastixUCoef ) {
                        lcoeftab[coefindx] = Lvalues[iterval];
                    }

                    if ( (side != PastixLCoef) &&
                         (rownum > (solvcblk->fcolnum + itercoltab)) )
                    {
#if defined(PRECISION_z) || defined(PRECISION_c)
                        if (bcsc->mtxtype == PastixHermitian)
                            ucoeftab[coefindx] = conj(Uvalues[iterval]);
                        else
#endif
                            ucoeftab[coefindx] = Uvalues[iterval];
                    }
                }
                else {
#if defined(PASTIX_DEBUG_COEFTAB)
                    fprintf(stderr, "cpucblk_zfillin: drop coeff from CSC c=%ld(%ld) l=%ld(%ld) cblk=%ld fcol=%ld lcol=%ld\n",
                            (long)solvcblk->fcolnum + itercoltab, (long)itercoltab,
                            (long)rownum, (long)iterval, (long)itercblk,
                            (long)solvcblk->fcolnum, (long)solvcblk->lcolnum );
#endif
                }
            }
        }
    }
}

/**
 *******************************************************************************
 *
 * @brief Initialize the low-rank coeftab structure from the internal bcsc.
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
 *          PaStiX structure to store numerical data and flags
 *
 * @param[in] bcsc
 *          The internal bcsc structure that hold the graph with permutation
 *          stored by cblk.
 *
 * @param[in] itercblk
 *          The index of the cblk to fill in both bcsc and solvmtx structures.
 *
 *******************************************************************************/
static inline void
cpucblk_zfillin_lr( pastix_coefside_t    side,
                    const SolverMatrix  *solvmtx,
                    const pastix_bcsc_t *bcsc,
                    pastix_int_t         itercblk )
{
    const bcsc_cblk_t *csccblk = bcsc->cscftab + itercblk;
    SolverCblk *solvcblk = solvmtx->cblktab + itercblk;
    SolverBlok *solvblok;
    SolverBlok *lsolvblok = (solvcblk+1)->fblokptr;
    pastix_complex64_t *lcoeftab, *ucoeftab;
    pastix_complex64_t *Lvalues  = bcsc->Lvalues;
    pastix_complex64_t *Uvalues  = bcsc->Uvalues;
    pastix_int_t itercoltab, iterval, coefindx, ldd;

    assert( solvcblk->cblktype & CBLK_LAYOUT_2D );

    for (itercoltab=0; itercoltab<csccblk->colnbr; itercoltab++)
    {
        pastix_int_t frow = csccblk->coltab[itercoltab];
        pastix_int_t lrow = csccblk->coltab[itercoltab+1];

        solvblok = solvcblk->fblokptr;
        ldd = blok_rownbr( solvblok );
        lcoeftab = (pastix_complex64_t*)(solvblok->LRblock[0].u);
        ucoeftab = (pastix_complex64_t*)(solvblok->LRblock[1].u);

        for (iterval=frow; iterval<lrow; iterval++)
        {
            pastix_int_t rownum = bcsc->rowtab[iterval];

            /* If values in the lower part of the matrix */
            if (rownum >= (solvcblk->fcolnum+itercoltab))
            {
                while ((solvblok < lsolvblok) &&
                       ((solvblok->lrownum < rownum) ||
                        (solvblok->frownum > rownum)))
                {
                    solvblok++;
                    ldd = blok_rownbr( solvblok );
                    lcoeftab = (pastix_complex64_t*)(solvblok->LRblock[0].u);
                    ucoeftab = (pastix_complex64_t*)(solvblok->LRblock[1].u);
                }

                if ( solvblok < lsolvblok )
                {
                    coefindx  = rownum - solvblok->frownum; /* Row shift    */
                    coefindx += itercoltab * ldd;           /* Column shift */

                    if ( side != PastixUCoef ) {
                        lcoeftab[coefindx] = Lvalues[iterval];
                    }

                    if ( (side != PastixLCoef) &&
                         (rownum > (solvcblk->fcolnum + itercoltab)) )
                    {
#if defined(PRECISION_z) || defined(PRECISION_c)
                        if (bcsc->mtxtype == PastixHermitian)
                            ucoeftab[coefindx] = conj(Uvalues[iterval]);
                        else
#endif
                            ucoeftab[coefindx] = Uvalues[iterval];
                    }
                }
                else {
#if defined(PASTIX_DEBUG_COEFTAB)
                    fprintf(stderr, "cpucblk_zfillin: drop coeff from CSC c=%ld(%ld) l=%ld(%ld) cblk=%ld fcol=%ld lcol=%ld\n",
                            (long)solvcblk->fcolnum + itercoltab, (long)itercoltab,
                            (long)rownum, (long)iterval, (long)itercblk,
                            (long)solvcblk->fcolnum, (long)solvcblk->lcolnum );
#endif
                }
            }
        }
    }
}

/**
 *******************************************************************************
 *
 * @brief Initialize the coeftab structure from the internal bcsc.
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
 *          PaStiX structure to store numerical data and flags
 *
 * @param[in] bcsc
 *          The internal bcsc structure that hold the graph with permutation
 *          stored by cblk.
 *
 * @param[in] itercblk
 *          The index of the cblk to fill in both bcsc and solvmtx structures.
 *
 *******************************************************************************/
void
cpucblk_zfillin( pastix_coefside_t    side,
                 const SolverMatrix  *solvmtx,
                 const pastix_bcsc_t *bcsc,
                 pastix_int_t         itercblk )
{
    if ( (solvmtx->cblktab + itercblk)->cblktype & CBLK_COMPRESSED ) {
        cpucblk_zfillin_lr( side, solvmtx, bcsc, itercblk );
    }
    else {
        cpucblk_zfillin_fr( side, solvmtx, bcsc, itercblk );
    }
}

