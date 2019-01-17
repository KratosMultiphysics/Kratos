/**
 *
 * @file symbol_cost.c
 *
 * PaStiX symbol structure cost functions
 *
 * @copyright 1999-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.1
 * @author Pascal Henon
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @date 2018-07-16
 *
 **/
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <assert.h>

#include "common.h"
#include "symbol.h"
#include "symbol_cost.h"

/**
 *******************************************************************************
 *
 * @ingroup symbol_dev_cost
 *
 * @brief Template function to compute cost on a column-block based approach
 *        with a single update per column block
 *
 *******************************************************************************
 *
 * @param[in] fptr
 *          The structure that contains the cost functions (diag, trsm and
 *          update are used)
 *
 * @param[in] symbmtx
 *          The symbolic matrix structure on which to compute the costs.
 *
 * @param[in] cblknum
 *          The index of the column-block for which the cost will be computed
 *
 *******************************************************************************
 *
 * @return The cost associated to the cblk of index cblknum and evaluated with
 * the set of given functions
 *
 *******************************************************************************/
static double
sum1d( const symbol_function_t *fptr,
       const symbol_matrix_t   *symbmtx,
             pastix_int_t       cblknum )
{
    symbol_cblk_t *cblk = symbmtx->cblktab + cblknum;
    pastix_int_t M, N, k;
    double nbops = 0.;
    double dof = (double)(symbmtx->dof);

    /*
     * Size of the factorization kernel (square)
     */
    N = (cblk->lcolnum - cblk->fcolnum + 1);

    /*
     * Height of the TRSM to which apply the TRSM
     */
    M = 0;
    for(k = cblk[0].bloknum+1; k < cblk[1].bloknum; k++)
    {
        M += (symbmtx->bloktab[k].lrownum - symbmtx->bloktab[k].frownum + 1);
    }

    if ( dof > 0. ) {
        N *= dof;
        M *= dof;
    }

    nbops = fptr->diag( N );
    if( M > 0 ) {
        nbops += fptr->trsm( M, N );
        nbops += fptr->update( N, M );
    }

    return nbops;
}

/**
 *******************************************************************************
 *
 * @ingroup symbol_dev_cost
 *
 * @brief Template function to compute cost on block based approach.
 *
 * As opposed to sum1d(), the updates are split in one per off-diagonal block
 * making it more precise to evaluate the performance cost of the GEMMs, for
 * example, as it exactly follow the 1D scheme used in the static scheduler of
 * PaStiX.
 *
 *******************************************************************************
 *
 * @param[in] fptr
 *          The structure that contains the cost functions (diag, trsm and
 *          blkupdate are used)
 *
 * @param[in] symbmtx
 *          The symbolic matrix structure on which to compute the costs.
 *
 * @param[in] cblknum
 *          The index of the column-block for which the cost will be computed
 *
 *******************************************************************************
 *
 * @return The cost associated to the cblk of index cblknum and evaluated with
 * the set of given functions
 *
 *******************************************************************************/
static double
sum2d( const symbol_function_t *fptr,
       const symbol_matrix_t      *symbmtx,
             pastix_int_t       cblknum )
{
    symbol_cblk_t *cblk = symbmtx->cblktab + cblknum;
    pastix_int_t M, N, K, l;
    double nbops = 0.;
    double dof = (double)(symbmtx->dof);

    /*
     * Size of the factorization kernel (square)
     */
    N = (cblk->lcolnum - cblk->fcolnum + 1);

    /*
     * Height of the TRSM to which apply the TRSM
     */
    M = 0;
    for(l = cblk[0].bloknum+1; l < cblk[1].bloknum; l++)
    {
        M += (symbmtx->bloktab[l].lrownum - symbmtx->bloktab[l].frownum + 1);
    }

    N *= dof;
    M *= dof;

    nbops = fptr->diag( N );
    if( M > 0 ) {
        nbops += fptr->trsm( M, N );
    }

    /*
     * Compute the cost of each GEMM
     */
    K = N;
    for(l = cblk[0].bloknum+1; l < cblk[1].bloknum; l++)
    {
        N = (symbmtx->bloktab[l].lrownum - symbmtx->bloktab[l].frownum + 1);
        N *= dof;

        nbops += fptr->blkupdate( K, M, N );

        M -= N;
    }
    return nbops;
}

/**
 *******************************************************************************
 *
 * @ingroup symbol_dev_cost
 *
 * @brief Template function to compute cost on block based approach which keeps
 * the cost per block.
 *
 * As opposed to sum2d(), the cost of each update per block is stored in the
 * blokcost array. Despite this storage, the function is completely identical.
 *
 *******************************************************************************
 *
 * @param[in] fptr
 *          The structure that contains the cost functions (diag, trsm and
 *          blkupdate are used)
 *
 * @param[in] symbmtx
 *          The symbolic matrix structure on which to compute the costs.
 *
 * @param[in] cblknum
 *          The index of the column-block for which the cost will be computed
 *
 * @param[inout] blokcost
 *          An array of size the number of blocks in the cblknum column-block in
 *          which to store the cost per block of each updates.
 *
 *******************************************************************************
 *
 * @return The cost associated to the cblk of index cblknum and evaluated with
 * the set of given functions
 *
 *******************************************************************************/
static double
sum2dext( const symbol_function_t *fptr,
          const symbol_matrix_t   *symbmtx,
                pastix_int_t       cblknum,
                double            *blokcost )
{
    symbol_cblk_t *cblk = symbmtx->cblktab + cblknum;
    pastix_int_t M, N, K, l;
    double nbops = 0.;
    double dof = (double)(symbmtx->dof);

    /*
     * Size of the factorization kernel (square)
     */
    N = (cblk->lcolnum - cblk->fcolnum + 1);

    /*
     * Height of the TRSM to which apply the TRSM
     */
    M = 0;
    for(l = cblk[0].bloknum+1; l < cblk[1].bloknum; l++)
    {
        M += (symbmtx->bloktab[l].lrownum - symbmtx->bloktab[l].frownum + 1);
    }

    N *= dof;
    M *= dof;

    nbops = fptr->diag( N );
    if( M > 0 ) {
        nbops += fptr->trsm( M, N );
    }
    *blokcost = nbops;
    blokcost++;

    /*
     * Compute the cost of each GEMM
     */
    K = N;
    for(l = cblk[0].bloknum+1; l < cblk[1].bloknum; l++, blokcost++)
    {
        N = (symbmtx->bloktab[l].lrownum - symbmtx->bloktab[l].frownum + 1);
        N *= dof;

        *blokcost = fptr->blkupdate( K, M, N );
        nbops += *blokcost;

        M -= N;
    }
    return nbops;
}

/**
 *******************************************************************************
 *
 * @ingroup symbol_dev_cost
 *
 * @brief Recursive function to compute the cost of the full symbolic structure
 * with either sum1d(), sum2d(), or sum2dext().
 *
 *******************************************************************************
 *
 * @param[in] a
 *          The first column-block index of the range to address
 *
 * @param[in] b
 *          The last column-block index of the range to address (inclusive)
 *
 * @param[in] fval
 *          The function to use to compute the cost. This can be sum1d(),
 *          sum2d(), or sum2dext()
 *
 * @param[in] fptr
 *          The set of functions that will be applied. It can be size function,
 *          floating point operation, performance models, ...
 *
 * @param[in] symbmtx
 *          The symbol matrix on which to compute the wanted information.
 *
 *******************************************************************************
 *
 * @return The cost associated to the integral of the symbol matrix structure
 *         and evaluated with the set of given functions
 *
 *******************************************************************************/
static double
recursive_sum( pastix_int_t a, pastix_int_t b,
               double (*fval)(const symbol_function_t *, const symbol_matrix_t *, pastix_int_t),
               const symbol_function_t *fptr,
               const symbol_matrix_t      *symbmtx )
{
    if(a != b)
        return recursive_sum(        a, (a+b)/2, fval, fptr, symbmtx)
            +  recursive_sum((a+b)/2+1,       b, fval, fptr, symbmtx);

    return fval(fptr, symbmtx, a);
}
/**
 * @}
 *
 * @addtogroup pastix_symbol
 * @{
 */


/**
 *******************************************************************************
 *
 * @brief Computes the number of non-zero elements in L.
 *
 * This computes the number of non-zero elements stored in the symbol matrix in
 * order to compute the fill-in. This routines returns the number of non-zero of
 * the strictly lower part of the matrix.
 *
 *******************************************************************************
 *
 * @param[in] symbptr
 *          The symbol structure to study.
 *
 *******************************************************************************
 *
 * @return  The number of non zero elements in the strictly lower part of the
 *          full symbol matrix.
 *
 *******************************************************************************/
pastix_int_t
pastixSymbolGetNNZ( const symbol_matrix_t *symbptr )
{
    symbol_cblk_t *cblk;
    symbol_blok_t *blok;
    pastix_int_t itercblk;
    pastix_int_t cblknbr;
    pastix_int_t nnz = 0;
    pastix_int_t dof = symbptr->dof;

    cblknbr = symbptr->cblknbr;
    cblk    = symbptr->cblktab;
    blok    = symbptr->bloktab;

    for(itercblk=0; itercblk<cblknbr; itercblk++, cblk++)
    {
        pastix_int_t iterblok = cblk[0].bloknum + 1;
        pastix_int_t lbloknum = cblk[1].bloknum;

        pastix_int_t colnbr = dof * (cblk->lcolnum - cblk->fcolnum + 1);

        /* Diagonal block */
        blok++;
        nnz += ( colnbr * (colnbr+1) ) / 2 - colnbr;

        /* Off-diagonal blocks */
        for( ; iterblok < lbloknum; iterblok++, blok++)
        {
            pastix_int_t rownbr = (blok->lrownum - blok->frownum + 1) * dof;

            nnz += rownbr * colnbr;
        }
    }

    return nnz;
}

/**
 *******************************************************************************
 *
 * @brief Computes the number of theoretical and real flops.
 *
 * Given a symbolic factorization structure, the type of factorization: Cholesky
 * or LU, and the arithmetic, this function will return the total number of
 * floating point operation that will be performed during the numerical
 * factorization.
 *
 *******************************************************************************
 *
 * @param[in] symbmtx
 *          The symbol structure to study.
 *
 * @param[in] flttype
 *          The floating type of the elements in the matrix.
 *          PastixPattern, PastixFloat, PastixDouble, PastixComplex32 or
 *          PastixComplex64. In case of PastixPattern, values for PastixDouble
 *          are returned.
 *
 * @param[in] factotype
 *          The factorization algorithm to perform: PastixFactLLT,
 *          PastixFactLDLT, PastixFactLLH, PastixFactLDLH or PastixFactLU.
 *
 * @param[out] thflops
 *          Returns the number of theoretical flops to perform.
 *          NULL if not asked.
 *
 * @param[out] rlflops
 *          Returns the number of real flops to perform, taking into account
 *          copies and scatter operations.
 *          NULL if not asked.
 *
 *******************************************************************************/
void
pastixSymbolGetFlops( const symbol_matrix_t *symbmtx,
                pastix_coeftype_t   flttype,
                pastix_factotype_t  factotype,
                double             *thflops,
                double             *rlflops )
{
    int iscomplex = ((flttype == PastixComplex32) || (flttype == PastixComplex64)) ? 1 : 0;

    /* Compute theoretical flops */
    if ( thflops != NULL ) {
        *thflops = recursive_sum(0, symbmtx->cblknbr-1, sum1d,
                                 &(flopstable[iscomplex][factotype]),
                                 symbmtx);
    }

    /* Compute performed flops */
    if ( rlflops != NULL ) {
        *rlflops = recursive_sum(0, symbmtx->cblknbr-1, sum2d,
                                 &(flopstable[iscomplex][factotype]),
                                 symbmtx);
    }
}

/**
 *******************************************************************************
 *
 * @brief Computes the cost of structure for the costMatrixBuild() function.
 *
 * This function iterates on the column-blocks and blocks to compute the cost of
 * the operation performed on each of those elements for the costMatrixBuild()
 * function that is used in the simulation for the data mapping. It returns an
 * array of cost for each type of element.
 *
 *******************************************************************************
 *
 * @param[in] symbmtx
 *          The symbol structure to study.
 *
 * @param[in] flttype
 *          The floating type of the elements in the matrix.
 *          PastixPattern, PastixFloat, PastixDouble, PastixComplex32 or
 *          PastixComplex64. In case of PastixPattern, values for PastixDouble
 *          are returned.
 *
 * @param[in] factotype
 *          The factorization algorithm to perform: PastixFactLLT,
 *          PastixFactLDLT, PastixFactLLH, PastixFactLDLH or PastixFactLU.
 *
 * @param[inout] cblkcost
 *          An allocated array of size cblknbr that will holds the cost per cblk
 *          on exit.
 *
 * @param[inout] blokcost
 *          An allocated array of size bloknbr that will holds the cost per blok
 *          on exit.
 *
 *******************************************************************************/
void
pastixSymbolGetTimes( const symbol_matrix_t *symbmtx,
                      pastix_coeftype_t      flttype,
                      pastix_factotype_t     factotype,
                      double                *cblkcost,
                      double                *blokcost )
{
    symbol_function_t *f;
    double *cblkptr, *blokptr;
    pastix_int_t i;
    int iscomplex = ((flttype == PastixComplex32) || (flttype == PastixComplex64)) ? 1 : 0;
    f = &(perfstable[iscomplex][factotype]);

    /* Initialize costs */
    cblkptr = cblkcost;
    blokptr = blokcost;

    for(i=0; i<symbmtx->cblknbr; i++, cblkptr++) {
        *cblkptr = sum2dext( f, symbmtx, i, blokptr );

        blokptr += symbmtx->cblktab[i+1].bloknum
            -      symbmtx->cblktab[i  ].bloknum;
    }

    assert( ( cblkptr - cblkcost ) == symbmtx->cblknbr );
    assert( ( blokptr - blokcost ) == symbmtx->bloknbr );
}

/**
 * @}
 */
