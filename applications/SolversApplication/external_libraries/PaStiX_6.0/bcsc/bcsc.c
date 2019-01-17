/**
 *
 * @file bcsc.c
 *
 * @copyright 2004-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.1
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @date 2018-07-16
 *
 **/
#include "common.h"
#include "pastix/order.h"
#include "spm.h"
#include "solver.h"
#include "bcsc.h"

#include "bcsc_z.h"
#include "bcsc_c.h"
#include "bcsc_d.h"
#include "bcsc_s.h"

/**
 *******************************************************************************
 *
 * @ingroup bcsc_internal
 *
 * @brief Initialize the coltab of a block csc matrix.
 *
 *******************************************************************************
 *
 * @param[in] solvmtx
 *          The solver matrix structure that describe the data distribution.
 *
 * @param[in] newcoltab
 *          Array of size spm->gN+1. This array is global coltab with -1 for non
 *          local indexes.
 *
 * @param[in] dof
 *          The degree of freedom of each unknown.
 *
 * @param[inout] bcsc
 *          On entry, the pointer to an allocated bcsc.
 *          On exit, the bcsc stores the ininitialized coltab split per block.
 *
 *******************************************************************************
 *
 * @return The number of non zero unknowns in the matrix.
 *
 *******************************************************************************/
static inline pastix_int_t
bcsc_init_coltab( const SolverMatrix  *solvmtx,
                  const pastix_int_t  *newcoltab,
                        pastix_int_t   dof,
                        pastix_bcsc_t *bcsc )
{
    bcsc_cblk_t *blockcol;
    pastix_int_t index, iter, idxcol, nodeidx, colsize;

    bcsc->cscfnbr = solvmtx->cblknbr;
    MALLOC_INTERN( bcsc->cscftab, bcsc->cscfnbr, bcsc_cblk_t );

    idxcol = 0;
    blockcol = bcsc->cscftab;
    for (index=0; index<bcsc->cscfnbr; index++, blockcol++)
    {
        pastix_int_t fcolnum = solvmtx->cblktab[index].fcolnum;
        pastix_int_t lcolnum = solvmtx->cblktab[index].lcolnum;

        blockcol->colnbr = (lcolnum - fcolnum + 1);
        MALLOC_INTERN( blockcol->coltab, blockcol->colnbr + 1, pastix_int_t );

        /* Works only for DoF constant */
        assert( fcolnum % dof == 0 );

        blockcol->coltab[0] = idxcol;
        for (iter=0; iter < blockcol->colnbr; iter++)
        {
            nodeidx = ( fcolnum + (iter-iter%dof) ) / dof;

            colsize = (newcoltab[nodeidx+1] - newcoltab[nodeidx]) * dof;
            blockcol->coltab[iter+1] = blockcol->coltab[iter] + colsize;
        }

        idxcol = blockcol->coltab[blockcol->colnbr];
    }

    if ( idxcol > 0 ) {
        MALLOC_INTERN( bcsc->rowtab,  idxcol, pastix_int_t);
        MALLOC_INTERN( bcsc->Lvalues, idxcol * pastix_size_of( bcsc->flttype ), char );
    }
    else {
        bcsc->rowtab  = NULL;
        bcsc->Lvalues = NULL;
    }
    bcsc->Uvalues = NULL;

    return idxcol;
}

/**
 *******************************************************************************
 *
 * @ingroup bcsc_internal
 *
 * @brief Restore the coltab array
 *
 * Function to restore the coltab array when it has been modified to initialize
 * the row and values arrays.
 *
 *******************************************************************************
 *
 * @param[inout] bcsc
 *          On entry, the bcsc to restore.
 *          On exit, the coltab array of the bcsc is restored to the correct
 *          indexes.
 *
 *******************************************************************************/
void
bcsc_restore_coltab( pastix_bcsc_t *bcsc )
{
    bcsc_cblk_t *blockcol;
    pastix_int_t index, iter, idxcol, idxcoltmp;

    idxcol = 0;
    blockcol = bcsc->cscftab;
    for (index=0; index<bcsc->cscfnbr; index++, blockcol++)
    {
        for (iter=0; iter <= blockcol->colnbr; iter++)
        {
            idxcoltmp = blockcol->coltab[iter];
            blockcol->coltab[iter] = idxcol;
            idxcol = idxcoltmp;
        }
    }
    return;
}

/**
 *******************************************************************************
 *
 * @brief Initialize the coltab of a centralized block csc matrix.
 *
 *******************************************************************************
 *
 * @param[in] spm
 *          The initial sparse matrix in the spm format.
 *
 * @param[in] ord
 *          The ordering that needs to be applied on the spm to generate the
 *          block csc.
 *
 * @param[in] solvmtx
 *          The solver matrix structure that describe the data distribution.
 *
 * @param[inout] bcsc
 *          On entry, the pointer to an allocated bcsc.
 *          On exit, the bcsc stores the input spm with the permutation applied
 *          and grouped accordingly to the distribution described in solvmtx.
 *
 *******************************************************************************
 *
 * @return The number of non zero unknowns in the matrix.
 *
 *******************************************************************************/
pastix_int_t
bcsc_init_centralized_coltab( const spmatrix_t     *spm,
                              const pastix_order_t *ord,
                              const SolverMatrix   *solvmtx,
                                    pastix_bcsc_t  *bcsc )
{
    pastix_int_t  valuesize, baseval;
    pastix_int_t *globcol  = NULL;
    pastix_int_t *colptr = spm->colptr;
    pastix_int_t *rowptr = spm->rowptr;
    int dof = spm->dof;
    int sym = (spm->mtxtype == SpmSymmetric) || (spm->mtxtype == SpmHermitian);

    bcsc->mtxtype = spm->mtxtype;
    baseval = spm->colptr[0];

    /*
     * Allocate and initialize globcol that contains the number of elements in
     * each column of the input matrix
     * Globcol is equivalent to the classic colptr for the internal blocked
     * csc. The blocked csc integrate the perumtation computed within order
     * structure.
     */
    MALLOC_INTERN( globcol, spm->gN+1, pastix_int_t );
    memset( globcol, 0, (spm->gN+1) * sizeof(pastix_int_t) );

    assert( spm->loc2glob == NULL );

    {
        pastix_int_t itercol, newcol;

        for (itercol=0; itercol<spm->gN; itercol++)
        {
            pastix_int_t frow = colptr[itercol]   - baseval;
            pastix_int_t lrow = colptr[itercol+1] - baseval;
            newcol = ord->permtab[itercol];
            globcol[newcol] += lrow - frow;

            assert( (lrow - frow) >= 0 );
            if (sym) {
                pastix_int_t iterrow, newrow;

                for (iterrow=frow; iterrow<lrow; iterrow++)
                {
                    pastix_int_t tmprow = rowptr[iterrow] - baseval;
                    if (tmprow != itercol) {
                        newrow = ord->permtab[tmprow];
                        globcol[newrow]++;
                    }
                }
            }
        }

        /* Compute displacements to update the colptr array */
        {
            pastix_int_t tmp, idx;

            idx = 0;
            for (itercol=0; itercol<=spm->gN; itercol++)
            {
                tmp = globcol[itercol];
                globcol[itercol] = idx;
                idx += tmp;
            }
        }
    }

    valuesize = bcsc_init_coltab( solvmtx, globcol, dof, bcsc );
    memFree_null( globcol );

    return valuesize;
}

/**
 *******************************************************************************
 *
 * @brief Initialize a centralized block csc when no MPI processes are involved.
 *
 *******************************************************************************
 *
 * @param[in] spm
 *          The initial sparse matrix in the spm format.
 *
 * @param[in] ord
 *          The ordering that needs to be applied on the spm to generate the
 *          block csc.
 *
 * @param[in] solvmtx
 *          The solver matrix structure that describe the data distribution.
 *
 * @param[in] initAt
 *          A flag to enable/disable the initialization of A'
 *
 * @param[inout] bcsc
 *          On entry, the pointer to an allocated bcsc.
 *          On exit, the bcsc stores the input spm with the permutation applied
 *          and grouped accordingly to the distribution described in solvmtx.
 *
 *******************************************************************************/
void
bcsc_init_centralized( const spmatrix_t     *spm,
                       const pastix_order_t *ord,
                       const SolverMatrix   *solvmtx,
                             pastix_int_t    initAt,
                             pastix_bcsc_t  *bcsc )
{
    pastix_int_t  itercol, itercblk;
    pastix_int_t  cblknbr  = solvmtx->cblknbr;
    pastix_int_t  eltnbr   = spm->gNexp;
    pastix_int_t *col2cblk = NULL;

    bcsc->mtxtype = spm->mtxtype;
    bcsc->flttype = spm->flttype;
    bcsc->gN      = spm->gN;
    bcsc->n       = spm->n;

    assert( spm->loc2glob == NULL );

    /*
     * Initialize the col2cblk array. col2cblk[i] contains the cblk index of the
     * i-th column. col2cblk[i] = -1 if not local.
     */
    {
        SolverCblk *cblk = solvmtx->cblktab;

        MALLOC_INTERN( col2cblk, eltnbr, pastix_int_t );
        for (itercol=0; itercol<eltnbr; itercol++)
        {
            col2cblk[itercol] = -1;
        }

        for (itercblk=0; itercblk<cblknbr; itercblk++, cblk++)
        {
            for (itercol  = cblk->fcolnum;
                 itercol <= cblk->lcolnum;
                 itercol++ )
            {
                col2cblk[itercol] = itercblk;
            }
        }
    }

    /*
     * Fill in the lower triangular part of the blocked csc with values and
     * rows. The upper triangular part is done later if required through LU
     * factorization.
     */
    switch( spm->flttype ) {
    case SpmFloat:
        bcsc_sinit_centralized( spm, ord, solvmtx, col2cblk, initAt, bcsc );
        break;
    case SpmDouble:
        bcsc_dinit_centralized( spm, ord, solvmtx, col2cblk, initAt, bcsc );
        break;
    case SpmComplex32:
        bcsc_cinit_centralized( spm, ord, solvmtx, col2cblk, initAt, bcsc );
        break;
    case SpmComplex64:
        bcsc_zinit_centralized( spm, ord, solvmtx, col2cblk, initAt, bcsc );
        break;
    case SpmPattern:
    default:
        fprintf(stderr, "bcsc_init_centralized: Error unknown floating type for input spm\n");
    }

    memFree_null(col2cblk);
}

/**
 *******************************************************************************
 *
 * @brief Initialize the block csc matrix.
 *
 * The block csc matrix is used to initialize the factorized matrix, and to
 * perform the matvec operations in refinement.
 *
 *******************************************************************************
 *
 * @param[in] spm
 *          The initial sparse matrix in the spm format.
 *
 * @param[in] ord
 *          The ordering that needs to be applied on the spm to generate the
 *          block csc.
 *
 * @param[in] solvmtx
 *          The solver matrix structure that describe the data distribution.
 *
 * @param[in] initAt
 *          A flag to enable/disable the initialization of A'
 *
 * @param[inout] bcsc
 *          On entry, the pointer to an allocated bcsc.
 *          On exit, the bcsc stores the input spm with the permutation applied
 *          and grouped accordingly to the distribution described in solvmtx.
 *
 *******************************************************************************
 *
 * @return The time spent to initialize the bcsc structure.
 *
 *******************************************************************************/
double
bcscInit( const spmatrix_t     *spm,
          const pastix_order_t *ord,
          const SolverMatrix   *solvmtx,
                pastix_int_t    initAt,
                pastix_bcsc_t  *bcsc )
{
    assert( ord->baseval == 0 );
    assert( ord->vertnbr == spm->n );

    double time = 0.;
    clockStart(time);

    if ( spm->loc2glob == NULL ) {
        bcsc_init_centralized( spm, ord, solvmtx, initAt, bcsc );
    }
    else {
        fprintf(stderr, "bcscInit: Distributed SPM not yet supported");
    }

    clockStop(time);
    return time;
}

/**
 *******************************************************************************
 *
 * @brief Cleanup the block csc structure but do not free the bcsc pointer.
 *
 *******************************************************************************
 *
 * @param[inout] bcsc
 *          The block csc matrix to free.
 *
 *******************************************************************************/
void
bcscExit( pastix_bcsc_t *bcsc )
{
    bcsc_cblk_t *cblk;
    pastix_int_t i;

    if ( bcsc->cscftab == NULL ) {
        return;
    }

    for (i=0, cblk=bcsc->cscftab; i < bcsc->cscfnbr; i++, cblk++ ) {
        memFree_null( cblk->coltab );
    }

    memFree_null( bcsc->cscftab );
    memFree_null( bcsc->rowtab );

    if ( (bcsc->Uvalues != NULL) &&
         (bcsc->Uvalues != bcsc->Lvalues) ) {
        memFree_null( bcsc->Uvalues );
    }

    memFree_null( bcsc->Lvalues );
}
