/**
 *
 * @file symbol_fax_graph.c
 *
 * PaStiX fax symbolic factorization routines fro Scotch esmumps library Part of
 * a parallel direct block solver. This is the block symbolic factorization
 * routine for graphs.
 *
 * @copyright 2004-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.1
 * @author Francois Pellegrini
 * @date 2018-07-16
 *
 *
 *  Dates:
 *   Version 0.0 - from 22 jul 1998 to 29 sep 1998
 *   Version 0.2 - from 08 may 2000 to 09 may 2000
 *   Version 1.0 - from 01 jun 2002 to 03 jun 2002
 *   Version 1.1 - from 26 jun 2002 to 26 jun 2002
 *   Version 2.0 - from 21 mar 2003 to 21 mar 2003
 *   Version 3.0 - from 02 mar 2004 to 02 mar 2004
 *
 **/
#include "common.h"
#include "symbol.h"
#include "pastix/order.h"
#include "symbol_fax.h"

/**
 *******************************************************************************
 *
 * @ingroup pastix_symbol
 *
 * @brief Compute the block symbolic factorization of
 * the given matrix graph according to the given vertex ordering.
 *
 * pastixSymbolFaxGraph() could have called pastixSymbolFax() in the regular way, as do all
 * of the grid-like factorization routines. However, for efficiency reasons, we
 * have decided to inline pastixSymbolFax(), to avoid a function call for every arc.
 *
 *******************************************************************************
 *
 * @param[inout] symbptr
 *          The symbolic matrix structure to fill in.
 *
 * @param[in] vertnbr
 *          The number of vertices.
 *
 * @param[in] verttab
 *          Array of size vertnbr+1
 *          The array of indirection to the edgetab for each vertex.
 *          edgetab[ verttab[i] ] to edgetab[ verttab[i+1] are the edges of the
 *          ith vertex.
 *
 * @param[in] edgetab
 *          Array of size edgenbr = verttab[vertnbr] - verttab[0]. The array of edges.
 *          edgetab[ verttab[i]   - verttab[0] ] to
 *          edgetab[ verttab[i+1] - verttab[0] ] are the edges of the ith vertex.
 *
 * @param[in] ordeptr
 *          The ordering structure that contains the permutation and inverse
 *          permutation vector, as well as the list of supernodes and the
 *          associated tree.
 *
 *******************************************************************************
 *
 * @retval 0  on success.
 * @retval !0 on failure.
 *
 *******************************************************************************/
int
pastixSymbolFaxGraph(       symbol_matrix_t *symbptr,
                            pastix_int_t     vertnbr,
                      const pastix_int_t    *verttab,
                      const pastix_int_t    *edgetab,
                      const pastix_order_t  *ordeptr )
{
#ifndef DOXYGEN_SHOULD_SKIP_THIS
    pastix_int_t baseval = verttab[0];
    pastix_int_t edgenbr = verttab[vertnbr] - baseval;
    const pastix_int_t * verttax;
    pastix_int_t         edgenum;
    const pastix_int_t * edgetax;

    verttax = verttab - baseval;
    edgetax = edgetab - baseval;

#define SYMBOL_FAX_ITERATOR(ngbdptr, vertnum, vertend)	\
    for (edgenum = verttax[vertnum];			\
         edgenum < verttax[vertnum + 1];                \
         edgenum ++) {					\
        vertend = edgetax[edgenum];

#define SYMBOL_FAX_VERTEX_DEGREE(ngbdptr, vertnum)	\
    (verttax[(vertnum) + 1] - verttax[(vertnum)])

    {
#define SYMBOL_FAX_INCLUDED
#include "symbol_fax.c"
    }
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
}
