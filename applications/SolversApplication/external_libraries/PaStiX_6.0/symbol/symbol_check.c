/**
 *
 * @file symbol_check.c
 *
 * PaStiX routines to check the correctness of the symbol structure.
 *
 * @copyright 1999-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.1
 * @author Francois Pellegrini
 * @author David Goudin
 * @author Pascal Henon
 * @author Francois Pellegrini
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @date 2018-07-16
 *
 **/
#include "common.h"
#include "symbol.h"

/**
 *******************************************************************************
 *
 * @ingroup pastix_symbol
 *
 * @brief Checks the consistency of the given symbolic block matrix.
 *
 * Because of incomplete factorization, from version 1.0, no check is performed
 * regarding the existence of facing blocks in facing columns.
 *
 * @todo Complete test set to check the brow information
 *
 *******************************************************************************
 *
 * @param[in] symbptr
 *          The symbol structure to check.
 *
 *******************************************************************************
 *
 * @retval 0 if the symbol matrix is correct
 * @retval 1 if incorrect
 *
 *******************************************************************************/
int
pastixSymbolCheck(const symbol_matrix_t * const  symbptr)
{
    pastix_int_t         baseval; /* Base value                           */
    const symbol_cblk_t *cblktax; /* Based access to cblktab              */
    pastix_int_t         cblkmax; /* Maximum column block index           */
    pastix_int_t         cblknum; /* Based number of current column block */
    const symbol_blok_t *bloktax; /* Based access to bloktab              */
    pastix_int_t         blokmax; /* Maximum block index                  */
    pastix_int_t         bloknum; /* Based number of current block        */
    pastix_int_t         nodemax; /* Maximum node index                   */

    baseval = symbptr->baseval;
    cblktax = symbptr->cblktab - baseval;
    cblkmax = symbptr->cblknbr + (baseval - 1);
    bloktax = symbptr->bloktab - baseval;
    blokmax = symbptr->bloknbr + baseval;
    nodemax = symbptr->nodenbr;

    for (cblknum = bloknum = baseval;
         cblknum <= cblkmax; cblknum ++) {
        if ((cblktax[cblknum].fcolnum     <  baseval)                  ||
            (cblktax[cblknum].lcolnum     >  nodemax)                  ||
            (cblktax[cblknum].bloknum     >  blokmax)                  ||
            (cblktax[cblknum].fcolnum     >  cblktax[cblknum].lcolnum) ||
            ((cblktax[cblknum].brownum != -1) &&
             (cblktax[cblknum + 1].brownum  <  cblktax[cblknum].brownum))||
            (cblktax[cblknum + 1].fcolnum-1 != cblktax[cblknum].lcolnum) ||
            (cblktax[cblknum + 1].bloknum   <= cblktax[cblknum].bloknum))
        {
            errorPrint ("symbolCheck: invalid column block array");
            assert(0);
            return     (1);
        }

        if ((bloktax[bloknum].frownum != cblktax[cblknum].fcolnum) ||
            (bloktax[bloknum].lrownum != cblktax[cblknum].lcolnum) ||
            (bloktax[bloknum].fcblknm != cblknum)) {
            errorPrint ("symbolCheck: invalid diagonal block");
            assert(0);
            return     (1);
        }

        for (bloknum ++; bloknum < cblktax[cblknum + 1].bloknum; bloknum ++) {
            if ((bloktax[bloknum].lcblknm != cblknum)                      ||
                (bloktax[bloknum].fcblknm <  baseval)                      ||
                (bloktax[bloknum].fcblknm >  cblkmax)                      ||
                (bloktax[bloknum].frownum <= bloktax[bloknum - 1].lrownum) ||
                (bloktax[bloknum].fcblknm <  bloktax[bloknum - 1].fcblknm)) {
                errorPrint ("symbolCheck: invalid block array");
                assert(0);
                return     (1);
            }
        }
    }

    assert( (cblktax[cblknum].brownum == -1) ||
            (cblktax[cblknum].brownum == (symbptr->bloknbr - symbptr->cblknbr)) );

    /* Check that schur starts within the symbol */
    if ( (symbptr->schurfcol < 0) || (symbptr->schurfcol > (symbptr->nodenbr + 1)) ) {
        errorPrint ("symbolCheck: invalid schurfcol");
        assert(0);
        return     (1);
    }
    return (0);
}
