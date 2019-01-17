/**
 *
 * @file symbol_expand.c
 *
 * @copyright 2012-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.1
 * @author Vincent Bridonneau
 * @date 2018-07-16
 *
 * @precisions normal z -> s d c p
 *
 **/
#include "common.h"
#include "symbol.h"
#include "symbol_reorder.h"

/**
 *******************************************************************************
 *
 * @ingroup symbol_dev
 *
 * @brief Compute the index of the first cblk belonging to the Schur complement
 *
 *******************************************************************************
 *
 * @param[inout] symbptr
 *          The symbol structure pointer in which to find the schur.
 *
 *******************************************************************************
 *
 * @return The index of the first cblk belonging to the Schur complement.
 *
 *******************************************************************************/
static inline pastix_int_t
symbol_expand_find_schurcblk( const symbol_matrix_t *symbptr )
{
    const symbol_cblk_t *cblk;

    cblk = symbptr->cblktab + symbptr->cblknbr;

    for (; cblk >= symbptr->cblktab; cblk--)
    {
        if ( cblk->fcolnum == symbptr->schurfcol ) {
            return cblk - symbptr->cblktab;
        }

        assert( cblk->fcolnum > symbptr->schurfcol );
    }

    return 0;
}

/**
 *******************************************************************************
 *
 * @ingroup symbol_dev
 *
 * @brief Expand the symbol matrix structure when the dof are variadic
 *
 *******************************************************************************
 *
 * @param[inout] symbptr
 *          The symbol structure to expand. On entry, the information is based
 *          on the compressed graph. On exit, the information are expanded to
 *          the size of the expanded matrix.
 *
 *******************************************************************************/
static inline void
symbol_expand_var( symbol_matrix_t *symbptr )
{
    pastix_int_t   i, baseval, cblknbr, bloknbr;
    symbol_cblk_t *cblk;
    symbol_blok_t *blok;
    const pastix_int_t *dofs;
    pastix_int_t   schuridx;

    dofs    = symbptr->dofs;
    baseval = symbptr->baseval;
    cblknbr = symbptr->cblknbr;
    bloknbr = symbptr->bloknbr;

    /*
     * Get the cblk index corresponding to schurfcol
     */
    schuridx = symbol_expand_find_schurcblk( symbptr );

    /*
     * Update cblks
     */
    cblk = symbptr->cblktab;
    for (i=0; i<cblknbr; i++, cblk++) {
        cblk->fcolnum = dofs[ cblk->fcolnum     - baseval ];
        cblk->lcolnum = dofs[ cblk->lcolnum + 1 - baseval ] - 1;
    }
    cblk->fcolnum = cblk[-1].lcolnum + 1;
    cblk->lcolnum = cblk[-1].lcolnum + 1;

    /*
     * Update bloks
     */
    blok = symbptr->bloktab;
    for (i=0; i<bloknbr; i++, blok++) {
        blok->frownum = dofs[ blok->frownum     - baseval ];
        blok->lrownum = dofs[ blok->lrownum + 1 - baseval ] - 1;
    }

    symbptr->nodenbr   = symbptr->cblktab[cblknbr ].lcolnum - baseval;
    symbptr->schurfcol = symbptr->cblktab[schuridx].fcolnum;
}

/**
 *******************************************************************************
 *
 * @ingroup symbol_dev
 *
 * @brief Expand the symbol matrix structure when the dof are constant
 *
 *******************************************************************************
 *
 * @param[inout] symbptr
 *          The symbol structure to expand. On entry, the information is based
 *          on the compressed graph. On exit, the information are expanded to
 *          the size of the expanded matrix.
 *
 *******************************************************************************/
static inline void
symbol_expand_fix( symbol_matrix_t *symbptr )
{
    pastix_int_t   col, cblknbr, bloknbr, row, dof;
    symbol_cblk_t *cblk;
    symbol_blok_t *blok;

    dof = symbptr->dof;

    /* Update cblks first */
    cblk    = symbptr->cblktab;
    cblknbr = symbptr->cblknbr;

    for (col = 0; col < cblknbr; col++, cblk++)
    {
        cblk->fcolnum =  cblk->fcolnum    * dof;
        cblk->lcolnum = (cblk->lcolnum+1) * dof - 1;
    }
    cblk->fcolnum = cblk[-1].lcolnum + 1;
    cblk->lcolnum = cblk[-1].lcolnum + 1;

    /* Update block row and column indexes */
    bloknbr = symbptr->bloknbr;
    blok    = symbptr->bloktab;

    for (row = 0; row < bloknbr; row++, blok++)
    {
        blok->frownum =  blok->frownum    * dof;
        blok->lrownum = (blok->lrownum+1) * dof - 1;
    }

    symbptr->nodenbr   *= dof;
    symbptr->schurfcol *= dof;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_symbol
 *
 * @brief Expand the symbol matrix structure based on the dof information
 * (compressed -> expanded)
 *
 *******************************************************************************
 *
 * @param[inout] symbptr
 *          The symbol structure to expand. On entry, the information is based
 *          on the compressed graph. On exit, the information are expanded to
 *          the size of the expanded matrix.
 *
 *******************************************************************************/
void
pastixSymbolExpand( symbol_matrix_t *symbptr )
{
    if ( symbptr == NULL ) {
        pastix_print_error( "pastixSymbolExpand: The symbol matrix is not initialized\n" );
        return;
    }

    pastixSymbolBase( symbptr, 0 );

    if ( symbptr->dof > 1 ) {
        symbol_expand_fix( symbptr );
    }
    else if ( symbptr->dof < 1 ) {
        symbol_expand_var( symbptr );
    }

    symbptr->dof = 1;
    memFree_null( symbptr->dofs );

    pastixSymbolCheck( symbptr );
}
