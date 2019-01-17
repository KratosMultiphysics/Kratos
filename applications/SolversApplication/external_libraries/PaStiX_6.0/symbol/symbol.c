/**
 *
 * @file symbol.c
 *
 * PaStiX symbol structure routines
 *
 * @copyright 2004-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.1
 * @author David Goudin
 * @author Francois Pellegrini
 * @author Mathieu Faverge
 * @author Pascal Henon
 * @author Pierre Ramet
 * @date 2018-07-16
 *
 **/
#include "common.h"
#include "graph.h"
#include "pastix/order.h"
#include "symbol.h"

/**
 *******************************************************************************
 *
 * @ingroup symbol_dev
 *
 * @brief Add a dof array to the symbol matrix if any.
 *
 * If the dof parameter is variadic, then a permuted version of the initial dof
 * array is constructed to match the symbol matrix that is working the permuted
 * matrix A.
 *
 *******************************************************************************
 *
 * @param[in] graph
 *          The original graph of the matrix
 *
 * @param[in] order
 *          The ordering structure describing the permutation of the unknowns in
 *          the compressed form.
 *
 * @param[inout] symbptr
 *          The symbol structure to which dof array must be added.
 *
 *******************************************************************************/
static inline void
symbol_init_adddofs( const pastix_graph_t *graph,
                     const pastix_order_t *order,
                     symbol_matrix_t      *symbptr )
{
    symbptr->dof  = graph->dof;
    symbptr->dofs = NULL;

    if ( symbptr->dof < 1 ) {
        pastix_int_t symbbase = symbptr->baseval;
        pastix_int_t ordebase = order->baseval;
        pastix_int_t i, ip, n, d, *dofs;

        n = graph->gN;

        MALLOC_INTERN( symbptr->dofs, n+1, pastix_int_t );

        dofs = symbptr->dofs;
        dofs[0] = symbbase;

        for(ip=0; ip<n; ip++, dofs++) {
            i = order->peritab[ip] - ordebase;

            assert( i < n );
            d = graph->dofs[i+1] - graph->dofs[i];

            dofs[1] = dofs[0] + d;
        }
        assert( (symbptr->dofs[n] - symbbase) == (graph->dofs[n] - graph->dofs[0]) );
    }

    return;
}

/**
 * @addtogroup pastix_symbol
 * @{
 *
 *******************************************************************************
 *
 * @brief Initialize the symbol structure.
 *
 * Initialized the permuted dof array if graph and order are provided. The
 * symbol is considered as dof = 1 otherwise.
 *
 *******************************************************************************
 *
 * @param[in] graph
 *          The original graph of the matrix
 *
 * @param[in] order
 *          The ordering structure describing the permutation of the unknowns in
 *          the compressed form.
 *
 * @param[inout] symbptr
 *          The symbol structure to initialize.
 *
 *******************************************************************************/
void
pastixSymbolInit ( const pastix_graph_t *graph,
                   const pastix_order_t *order,
                        symbol_matrix_t *symbptr )
{
    memset (symbptr, 0, sizeof (symbol_matrix_t));
    symbptr->dof = 1;
    symbptr->schurfcol = -1;

    if ( (graph != NULL) && (order != NULL) ) {
        symbol_init_adddofs( graph, order, symbptr );
    }

    return;
}

/**
 *******************************************************************************
 *
 * @brief Free the content of symbolic matrix.
 *
 * All the arrays from the structure are freed and the structure is memset to 0
 * at exit, but the symbol itself is not freed. It will require a new call to
 * pastixSymbolInit() if the memory space area needs to be reused for a new
 * symbol matrix.
 *
 *******************************************************************************
 *
 * @param[inout] symbptr
 *          The pointer to the structure to free.
 *
 *******************************************************************************/
void
pastixSymbolExit( symbol_matrix_t *symbptr )
{
    if (symbptr->dofs != NULL) {
        memFree_null( symbptr->dofs );
    }
    if (symbptr->cblktab != NULL) {
        memFree_null( symbptr->cblktab );
    }
    if (symbptr->bloktab != NULL) {
        memFree_null( symbptr->bloktab );
    }
    if (symbptr->browtab != NULL) {
        memFree_null( symbptr->browtab );
    }
    pastixSymbolInit( NULL, NULL, symbptr );
}

/**
 *******************************************************************************
 *
 * @brief Reallocate the data structure to optimize the memory alignment.
 *
 * This function is used when the symbol need to be shrinked to a smaller or
 * larger set of blocks and column blocks. The original data is copied in the
 * new arrays.
 *
 *******************************************************************************
 *
 * @param[inout] symbptr
 *          The pointer to the structure that needs to be reallocated.
 *
 *******************************************************************************/
void
pastixSymbolRealloc( symbol_matrix_t *symbptr )
{
    symbol_cblk_t *cblktab = NULL;
    symbol_blok_t *bloktab = NULL;

    /* Move column block array */
    MALLOC_INTERN( cblktab, symbptr->cblknbr+1, symbol_cblk_t );
    memcpy(cblktab, symbptr->cblktab, (symbptr->cblknbr + 1) * sizeof (symbol_cblk_t));
    memFree(symbptr->cblktab);
    symbptr->cblktab = cblktab;

    /* Move block array */
    MALLOC_INTERN( bloktab, symbptr->bloknbr, symbol_blok_t );
    memcpy(bloktab, symbptr->bloktab, (symbptr->bloknbr) * sizeof (symbol_blok_t));
    memFree(symbptr->bloktab);
    symbptr->bloktab = bloktab;
}

/**
 *******************************************************************************
 *
 * @brief Search the targeted block C for a couple of blocks A and B.
 *
 * When executing the simulation run to map the data on the cores, it requires
 * to compute dependencies of each block. In that case for each couple of blocks
 * A (defined by bloknum), and B (defined by bloksrc), we need to find the block
 * that will receive the contribution for the matrix-matrix product. To speedup
 * the search, the startsearch parameter can be given to specify that the index
 * of the block searched is larger than this parameter. It returns the index of
 * the C block when found, or -1 if no block is facing the update. This happens,
 * only if ILU(k) factorization is applied and ricar is set to true. Otherwise,
 * returning -1 is an error.
 *
 *******************************************************************************
 *
 * @param[in] symbptr
 *          The pointer to the structure that needs to be reallocated.
 *
 * @param[in] bloksrc
 *          Index of the block used as A in teh matrix-matrix product.
 *
 * @param[in] bloknum
 *          Index of the block used as B in teh matrix-matrix product.
 *
 * @param[in] startsearch
 *          Index of a block belonging to the facing cblk of B used as a
 *          starting point to find the index of C.
 *
 * @param[in] ricar
 *          Booleen to enable ILU(k) factorization or not
 *
 *******************************************************************************
 *
 * @retval -1 No block where found. This is an error if ricar is disabled,
 *            otherwise it means that the path is longer than the k parameter in
 *            the ILU(k) factorization.
 *
 * @retval i The index of the block C that will receive the contribution from
 *           A * B^t
 *
 *******************************************************************************/
pastix_int_t
pastixSymbolGetFacingBloknum( const symbol_matrix_t *symbptr,
                              pastix_int_t bloksrc,
                              pastix_int_t bloknum,
                              pastix_int_t startsearch,
                              int ricar )
{
    symbol_blok_t *bsrc;
    symbol_blok_t *bdst;
    pastix_int_t i, fcblknum, fbloknum, lbloknum;

    fcblknum = symbptr->bloktab[bloksrc].fcblknm;
    fbloknum = symbptr->cblktab[fcblknum].bloknum;
    lbloknum = symbptr->cblktab[fcblknum+1].bloknum;

    if(startsearch < fbloknum )
        startsearch = fbloknum;

    assert( startsearch < lbloknum );

    /* Block in original column block */
    bsrc = (symbptr->bloktab) + bloknum;

    /* Search for the facing block in the facing column block */
    bdst = (symbptr->bloktab) + startsearch;

    if(ricar == 0)
    {
        for(i=startsearch; i<lbloknum; i++, bdst++ )
            if( bdst->lrownum >= bsrc->frownum)
                break;

        /* We should always exit the loop in non ilu(k) mode */
        assert( (bdst->frownum <= bsrc->frownum) &&
                (bdst->lrownum >= bsrc->lrownum) );

        return i;
    }
    else
    {
        for(i=startsearch; i<lbloknum; i++, bdst++)
        {
            if( ((bsrc->frownum >= bdst->frownum) && (bsrc->frownum <= bdst->lrownum)) ||
                ((bsrc->lrownum >= bdst->frownum) && (bsrc->lrownum <= bdst->lrownum)) ||
                ((bsrc->frownum <= bdst->frownum) && (bsrc->lrownum >= bdst->lrownum)) )
                return i;  /** We found the first block that matches **/

            if(bsrc->lrownum < bdst->frownum)
            {
                return -1;
            }
        }
    }
    return -1;
}

/**
 *******************************************************************************
 *
 * @brief Construct the browtab array that stores the blocks in a CSR way.
 *
 * The browtab is an equivalent of the columns array in a CSR for the symbolic
 * structure in terms of block indexes.
 *
 *******************************************************************************
 *
 * @param[inout] symbptr
 *          The pointer to the symbolic structure to update.
 *
 *******************************************************************************/
void
pastixSymbolBuildRowtab( symbol_matrix_t *symbptr )
{
    symbol_cblk_t *cblk;
    symbol_blok_t *blok;
    pastix_int_t  *innbr, *intmp, *browtab;
    pastix_int_t   itercblk;
    pastix_int_t   cblknbr;
    pastix_int_t   edgenbr = symbptr->bloknbr - symbptr->cblknbr;

    cblknbr = symbptr->cblknbr;

    MALLOC_INTERN(innbr, cblknbr, pastix_int_t );
    memset( innbr, 0, cblknbr * sizeof(pastix_int_t) );

    /* Count the number of input edge per cblk */
    cblk = symbptr->cblktab;
    blok = symbptr->bloktab;
    for(itercblk=0; itercblk<cblknbr; itercblk++, cblk++)
    {
        pastix_int_t iterblok = cblk[0].bloknum + 1;
        pastix_int_t lbloknum = cblk[1].bloknum;

        /* Skip diagonal block */
        blok++;

        /* Off-diagonal blocks */
        for( ; iterblok < lbloknum; iterblok++, blok++)
        {
            innbr[ blok->fcblknm ]++;
        }
    }

    /* Initialize the brownum fields */
    cblk = symbptr->cblktab;
    cblk->brownum = 0;
    for(itercblk=0; itercblk<cblknbr; itercblk++, cblk++)
    {
        cblk[1].brownum = cblk[0].brownum + innbr[ itercblk ];
        innbr[itercblk] = cblk[0].brownum;
    }
    assert( cblk[0].brownum == edgenbr );

    /* Initialize the browtab */
    MALLOC_INTERN(browtab, edgenbr, pastix_int_t );

    cblk = symbptr->cblktab;
    blok = symbptr->bloktab;
    for(itercblk=0; itercblk<cblknbr; itercblk++, cblk++)
    {
        pastix_int_t iterblok = cblk[0].bloknum + 1;
        pastix_int_t lbloknum = cblk[1].bloknum;

        /* Skip diagonal block */
        blok++;

        /* Off-diagonal blocks */
        for( ; iterblok < lbloknum; iterblok++, blok++)
        {
            intmp = innbr + blok->fcblknm;
            browtab[ *intmp ] = iterblok;
            (*intmp)++;
        }
    }

    if (symbptr->browtab == NULL) {
        memFree(symbptr->browtab);
    }
    symbptr->browtab = browtab;

    memFree( innbr );
    return;
}

/**
 *******************************************************************************
 *
 * @brief Print statistical information about the symbolic matrix structure
 *
 *******************************************************************************
 *
 * @param[in] symbptr
 *          The pointer to the symbolic structure.
 *
 *******************************************************************************/
void
pastixSymbolPrintStats( const symbol_matrix_t *symbptr )
{
    symbol_cblk_t *cblk;
    symbol_blok_t *blok;
    pastix_int_t itercblk, dof;
    pastix_int_t cblknbr, bloknbr;
    pastix_int_t cblkmin, cblkmax;
    pastix_int_t blokmin, blokmax;
    double cblkavg1, blokavg1;
    double cblkavg2, blokavg2;
    size_t mem = 0;

    cblknbr  = symbptr->cblknbr;
    bloknbr  = symbptr->bloknbr - cblknbr;
    cblkmin  = PASTIX_INT_MAX;
    cblkmax  = 0;
    cblkavg1 = 0;
    cblkavg2 = 0;
    blokmin  = PASTIX_INT_MAX;
    blokmax  = 0;
    blokavg1 = 0;
    blokavg2 = 0;

    cblk = symbptr->cblktab;
    blok = symbptr->bloktab;

    for(itercblk=0; itercblk<cblknbr; itercblk++, cblk++)
    {
        pastix_int_t iterblok = cblk[0].bloknum + 1;
        pastix_int_t lbloknum = cblk[1].bloknum;

        pastix_int_t colnbr = cblk->lcolnum - cblk->fcolnum + 1;

        cblkmin = pastix_imin( cblkmin, colnbr );
        cblkmax = pastix_imax( cblkmax, colnbr );
        cblkavg1 += colnbr;
        cblkavg2 += colnbr * colnbr;
        blok++;

        /* Only extra diagonal */
        for( ; iterblok < lbloknum; iterblok++, blok++)
        {
            pastix_int_t rownbr = blok->lrownum - blok->frownum + 1;

            blokmin = pastix_imin( blokmin, rownbr );
            blokmax = pastix_imax( blokmax, rownbr );
            blokavg1 += rownbr;
            blokavg2 += rownbr * rownbr;
        }
    }

    dof = symbptr->dof;
    blokmin *= dof;
    blokmax *= dof;
    cblkmin *= dof;
    cblkmax *= dof;
    cblkavg1 = (cblkavg1 * (double)dof ) / (double)cblknbr;
    blokavg1 = (blokavg1 * (double)dof ) / (double)bloknbr;
    cblkavg2 = sqrt( ((cblkavg2 * (double)dof * (double)dof) / (double)cblknbr) - cblkavg1 * cblkavg1 );
    blokavg2 = sqrt( ((blokavg2 * (double)dof * (double)dof) / (double)bloknbr) - blokavg1 * blokavg1 );

    /* Compute symbol matrix space */
    mem = sizeof( symbol_matrix_t );
    mem += sizeof( symbol_cblk_t )   * (cblknbr + 1);
    mem += sizeof( symbol_blok_t )   * symbptr->bloknbr;
    mem += sizeof( pastix_int_t ) * bloknbr;

    fprintf(stdout,
            "    Symbol Matrix statistics:\n"
            "      Number of cblk                    %10ld\n"
            "      Number of blok                    %10ld\n"
            "      Cblk width min                    %10ld\n"
            "      Cblk width max                    %10ld\n"
            "      Cblk width avg                   %11.2lf\n"
            "      Cblk width stdev                 %11.2lf\n"
            "      Blok height min                   %10ld\n"
            "      Blok height max                   %10ld\n"
            "      Blok height avg                  %11.2lf\n"
            "      Blok height stdev                %11.2lf\n"
            "      Memory space                     %11.2lf %co\n",
            (long)cblknbr, (long)bloknbr,
            (long)cblkmin, (long)cblkmax, cblkavg1, cblkavg2,
            (long)blokmin, (long)blokmax, blokavg1, blokavg2,
            pastix_print_value( mem ),
            pastix_print_unit( mem ) );
}

/**
 * @}
 */
