/**
 *
 * @file symbol_io.c
 *
 * PaStiX symbol structure IO functions.
 *
 * @copyright 2004-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.1
 * @author David Goudin
 * @author Pascal Henon
 * @author Francois Pellegrini
 * @author Pierre Ramet
 * @date 2018-07-16
 *
 *  Dates:
 *   Version 0.0 - from 23 aug 1998 to 07 oct 1998
 *   Version 0.1 - from 21 mar 2002 to 21 mar 2002
 *   Version 1.0 - from 03 jun 2002 to 08 sep 2003
 *   Version 3.0 - from 29 feb 2004 to 29 feb 2004
 *
 * @addtogroup pastix_symbol
 * @{
 *
 */
#include "common.h"
#include "symbol.h"

/**
 *******************************************************************************
 *
 * @brief Load the given block matrix structure from the given stream.
 *
 *******************************************************************************
 *
 * @param[inout] symbptr
 *          The symbolic matrix structure to fill in.
 *
 * @param[inout] stream
 *          The stream from which to read the structure.
 *
 *******************************************************************************
 *
 * @retval 0  on success.
 * @retval !0 on failure.
 *
 *******************************************************************************/
int
pastixSymbolLoad ( symbol_matrix_t * const symbptr,
                   FILE            * const stream )
{
    pastix_int_t                 versval;
    pastix_int_t                 baseval;
    pastix_int_t                 nodenbr;
    pastix_int_t                 cblknbr;
    pastix_int_t                 cblknum;
    pastix_int_t                 bloknbr;
    pastix_int_t                 bloknum;

    if ((intLoad (stream, &versval) +               /* Read header */
         intLoad (stream, &cblknbr) +
         intLoad (stream, &bloknbr) +
         intLoad (stream, &nodenbr) +
         intLoad (stream, &baseval) != 5) ||
        (versval < 0)                     ||        /* Version should be 0 or 1 */
        (versval > 1)                     ||
        (bloknbr < cblknbr)               ||
        (nodenbr < cblknbr)) {
        errorPrint ("symbolLoad: bad input (1)");
        return     (1);
    }

    if (((symbptr->cblktab = (symbol_cblk_t *) memAlloc ((cblknbr + 1) * sizeof (symbol_cblk_t))) == NULL) ||
        ((symbptr->bloktab = (symbol_blok_t *) memAlloc ( bloknbr      * sizeof (symbol_blok_t))) == NULL)) {
        errorPrint ("symbolLoad: out of memory");
        pastixSymbolExit (symbptr);
        return     (1);
    }
    symbptr->baseval = baseval;
    symbptr->cblknbr = cblknbr;
    symbptr->bloknbr = bloknbr;
    symbptr->nodenbr = nodenbr;

    for (cblknum = 0; cblknum < cblknbr; cblknum ++) {
        if ((intLoad (stream, &symbptr->cblktab[cblknum].fcolnum) + /* Read column blocks */
             intLoad (stream, &symbptr->cblktab[cblknum].lcolnum) +
             intLoad (stream, &symbptr->cblktab[cblknum].bloknum) != 3) ||
            (symbptr->cblktab[cblknum].fcolnum > symbptr->cblktab[cblknum].lcolnum)) {
            errorPrint ("symbolLoad: bad input (2)");
            pastixSymbolExit (symbptr);
            return     (1);
        }
    }
    symbptr->cblktab[cblknbr].fcolnum =             /* Set last column block */
        symbptr->cblktab[cblknbr].lcolnum = nodenbr + baseval;
    symbptr->cblktab[cblknbr].bloknum = bloknbr + baseval;

    for (bloknum = 0; bloknum < bloknbr; bloknum ++) {
        if ((intLoad (stream, &symbptr->bloktab[bloknum].frownum) + /* Read column blocks */
             intLoad (stream, &symbptr->bloktab[bloknum].lrownum) +
             intLoad (stream, &symbptr->bloktab[bloknum].fcblknm) != 3) ||
            (symbptr->bloktab[bloknum].frownum > symbptr->bloktab[bloknum].lrownum)) {
            errorPrint ("symbolLoad: bad input (3)");
            pastixSymbolExit (symbptr);
            return     (1);
        }

        /* This block is required for old file compatibility: will be removed */
        if (0)
        {
            pastix_int_t tmp;
            if ((versval > 0) && (intLoad (stream, &tmp) != 1)) {
                errorPrint ("symbolLoad: bad input (4)");
                pastixSymbolExit (symbptr);
                return     (1);
            }
        }
    }

    // TODO: store the value and read it
    symbptr->dof = 1;
    return (0);
}

/**
 *******************************************************************************
 *
 * @brief Save the given block matrix structure to the given stream.
 *
 *******************************************************************************
 *
 * @param[inout] symbptr
 *          The symbolic matrix structure to write.
 *
 * @param[inout] stream
 *          The stream to which to write the structure.
 *
 *******************************************************************************
 *
 * @retval 0  on success.
 * @retval !0 on failure.
 *
 *******************************************************************************/
int
pastixSymbolSave( const symbol_matrix_t * const symbptr,
                  FILE                  * const stream )
{
    const symbol_cblk_t *cblktnd;
    const symbol_cblk_t *cblkptr;
    const symbol_blok_t *bloktnd;
    const symbol_blok_t *blokptr;
    int o;

    o = (fprintf (stream, "1\n%ld\t%ld\t%ld\t%ld\n", /* Write file header */
                  (long) symbptr->cblknbr,
                  (long) symbptr->bloknbr,
                  (long) symbptr->nodenbr,
                  (long) symbptr->baseval) == EOF);
    for (cblkptr = symbptr->cblktab, cblktnd = cblkptr + symbptr->cblknbr;
         (cblkptr < cblktnd) && (o == 0); cblkptr ++) {
        o = (fprintf (stream, "%ld\t%ld\t%ld\n",
                      (long) cblkptr->fcolnum,
                      (long) cblkptr->lcolnum,
                      (long) cblkptr->bloknum) == EOF);
    }
    for (blokptr = symbptr->bloktab, bloktnd = blokptr + symbptr->bloknbr;
         (blokptr < bloktnd) && (o == 0); blokptr ++) {
        o = (fprintf (stream, "%ld\t%ld\t%ld\n",
                      (long) blokptr->frownum,
                      (long) blokptr->lrownum,
                      (long) blokptr->fcblknm ) == EOF);
    }

    return (o);
}

/**
 *******************************************************************************
 *
 * @brief Print the given block matrix structure in human readable format.
 *
 *******************************************************************************
 *
 * @param[inout] symbptr
 *          The symbolic matrix structure to write.
 *
 * @param[inout] file
 *          The file stream to which to write the structure.
 *
 *******************************************************************************/
void
pastixSymbolPrint( const symbol_matrix_t *symbptr,
                   FILE *file )
{
    pastix_int_t i, j;
    for(i=0;i<symbptr->cblknbr;i++)
    {
        fprintf(file, "CBLK %ld [%ld : %ld ] \n",
                (long)i,
                (long)symbptr->cblktab[i].fcolnum,
                (long)symbptr->cblktab[i].lcolnum );
        for(j=symbptr->cblktab[i].bloknum;
            j<symbptr->cblktab[i+1].bloknum;j++) {
            fprintf(file, "--BLOK %ld [%ld : %ld ]\n",
                    (long)j,
                    (long)symbptr->bloktab[j].frownum,
                    (long)symbptr->bloktab[j].lrownum );
        }
        fprintf(file, "\n");
    }
}

/**
 * @}
 */
