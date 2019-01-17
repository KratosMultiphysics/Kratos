/**
 *
 * @file symbol_base.c
 *
 * PaStiX symbol structure base function
 *
 * @copyright 1999-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.1
 * @author Francois Pellegrini
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
 * @brief Sets the base of the given symbol matrix structure to the given
 * base value.
 *
 *******************************************************************************
 *
 * @param[inout] symbptr
 *          The symbol structure to initialize.
 *
 * @param[in] baseval
 *          The base value.
 *
 *******************************************************************************/
void
pastixSymbolBase ( symbol_matrix_t * const symbptr,
                   const pastix_int_t      baseval )
{
    symbol_cblk_t *cblk;
    symbol_blok_t *blok;
    pastix_int_t   baseadj; /* Base adjust */
    pastix_int_t   cblknum;
    pastix_int_t   bloknum;

    baseadj = baseval - symbptr->baseval; /* Set base adjust     */
    if (baseadj == 0)                     /* If base already set */
        return;

    symbptr->baseval = baseval;           /* Set graph base */

    symbptr->schurfcol += baseadj;

    cblk = symbptr->cblktab;
    for (cblknum = 0; cblknum <= symbptr->cblknbr; cblknum ++, cblk++) {
        cblk->fcolnum += baseadj;
        cblk->lcolnum += baseadj;
        cblk->bloknum += baseadj;
    }

    blok = symbptr->bloktab;
    for (bloknum = 0; bloknum < symbptr->bloknbr; bloknum ++) {
        blok->frownum += baseadj;
        blok->lrownum += baseadj;
        blok->lcblknm += baseadj;
        blok->fcblknm += baseadj;
    }

    if ( symbptr->dof < 1 ) {
        pastix_int_t  i;
        pastix_int_t *dofs = symbptr->dofs;

        assert( dofs != NULL );

        for (i = 0; i <= symbptr->nodenbr; i ++, dofs++) {
            *dofs += baseadj;
        }
    }
}
