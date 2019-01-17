/**
 *
 * @file z_spm_expand.c
 *
 * SParse Matrix package random multi-dof spm generator.
 *
 * @copyright 2016-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Alban Bellot
 * @date 2015-01-01
 *
 * @precisions normal z -> c d s p
 **/
#include "common.h"
#include "z_spm.h"

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_dof
 *
 * @brief Expand a single dof CSC to a multi-dofs CSC.
 *
 * Each element matrix is fully initialized with the same element as the
 * original one.
 *
 *******************************************************************************
 *
 * @param[in] spm
 *           The original sparse CSC matrix used as the template for the
 *           muti-dof matrix.
 *
 *******************************************************************************
 *
 * @return The expanded CSC matrix according to the dofs properties previously
 * set.
 *
 *******************************************************************************/
spmatrix_t *
z_spmCSCExpand(const spmatrix_t *spm)
{
    spmatrix_t       *newspm;
    spm_int_t        i, j, k, ii, jj, dofi, dofj, col, row, baseval, lda;
    spm_int_t        diag, height;
    spm_int_t       *newcol, *newrow, *oldcol, *oldrow, *dofs;
#if !defined(PRECISION_p)
    spm_complex64_t *newval = NULL;
#endif
    spm_complex64_t *oldval2, *oldval = NULL;

    assert( spm->fmttype == SpmCSC );
    assert( spm->flttype == SpmComplex64 );

    if ( spm->dof == 1 ) {
        return (spmatrix_t*)spm;
    }

    if ( spm->layout != SpmColMajor ) {
        fprintf( stderr, "Unsupported layout\n" );
        return NULL;
    }

    newspm = malloc( sizeof(spmatrix_t) );
    memcpy( newspm, spm, sizeof(spmatrix_t) );

    baseval = spmFindBase( spm );
    oldcol = spm->colptr;
    oldrow = spm->rowptr;
    dofs   = spm->dofs;
#if !defined(PRECISION_p)
    oldval = oldval2 = (spm_complex64_t*)(spm->values);
#endif

    newspm->n = spm->nexp;
    newspm->colptr = newcol = malloc(sizeof(spm_int_t)*(spm->nexp+1));

    /**
     * First loop to compute the new colptr
     */
    *newcol = baseval;
    for(j=0; j<spm->n; j++, oldcol++)
    {
        diag = 0;
        dofj = (spm->dof > 0 ) ? spm->dof : dofs[j+1] - dofs[j];

        /* Sum the heights of the elements in the column */
        newcol[1] = newcol[0];
        for(k=oldcol[0]; k<oldcol[1]; k++)
        {
            i = oldrow[k-baseval] - baseval;
            dofi = (spm->dof > 0 ) ? spm->dof : dofs[i+1] - dofs[i];
            newcol[1] += dofi;

            diag = (diag || (i == j));
        }

        diag = (diag & (spm->mtxtype != SpmGeneral));
        height = newcol[1] - newcol[0];
        newcol++;

        /* Add extra columns */
        for(jj=1; jj<dofj; jj++, newcol++)
        {
            newcol[1] = newcol[0] + height;

            if ( diag ) {
                newcol[1] -= jj;
            }
        }
    }
    assert( ((spm->mtxtype == SpmGeneral) && ((newcol[0]-baseval) == spm->nnzexp)) ||
            ((spm->mtxtype != SpmGeneral) && ((newcol[0]-baseval) <= spm->nnzexp)) );

    newspm->nnz = newcol[0] - baseval;
    newspm->rowptr = newrow = malloc(sizeof(spm_int_t)*newspm->nnz);
#if !defined(PRECISION_p)
    newspm->values = newval = malloc(sizeof(spm_complex64_t)*newspm->nnz);
#endif

    /**
     * Second loop to compute the new rowptr and valptr
     */
    oldcol = spm->colptr;
    oldrow = spm->rowptr;
    newcol = newspm->colptr;
    for(j=0, col=0; j<spm->n; j++, oldcol++)
    {
        /**
         * Backup current position in oldval because we will pick
         * interleaved data inside the buffer
         */
        lda = newcol[1] - newcol[0];
        oldval2 = oldval;

        if ( spm->dof > 0 ) {
            dofj = spm->dof;
            assert( col == spm->dof * j );
        }
        else {
            dofj = dofs[j+1] - dofs[j];
            assert( col == (dofs[j] - baseval) );
        }

        for(jj=0; jj<dofj; jj++, col++, newcol++)
        {
            assert( ((spm->mtxtype == SpmGeneral) && (lda == (newcol[1] - newcol[0]))) ||
                    ((spm->mtxtype != SpmGeneral) && (lda >= (newcol[1] - newcol[0]))) );

            /* Move to the top of the column jj in coming element (i,j) */
            oldval = oldval2;

            for(k=oldcol[0]; k<oldcol[1]; k++)
            {
                i = oldrow[k-baseval] - baseval;

                if ( spm->dof > 0 ) {
                    dofi = spm->dof;
                    row  = spm->dof * i;
                }
                else {
                    dofi = dofs[i+1] - dofs[i];
                    row  = dofs[i] - baseval;
                }

                /* Move to the top of the jj column in the current element */
                oldval += dofi * jj;

                for(ii=0; ii<dofi; ii++, row++)
                {
                    if ( (spm->mtxtype == SpmGeneral) ||
                         (i != j) ||
                         ((i == j) && (row >= col)) )
                    {
                        (*newrow) = row + baseval;
                        newrow++;
#if !defined(PRECISION_p)
                        (*newval) = *oldval;
                        newval++;
#endif
                    }
                    oldval++;
                }
                /* Move to the top of the next element */
                oldval += dofi * (dofj-jj-1);
            }
        }
    }

    newspm->gN      = newspm->n;
    newspm->gnnz    = newspm->nnz;

    newspm->gNexp   = newspm->gN;
    newspm->nexp    = newspm->n;
    newspm->gnnzexp = newspm->gnnz;
    newspm->nnzexp  = newspm->nnz;

    newspm->dof     = 1;
    newspm->dofs    = NULL;
    newspm->layout  = SpmColMajor;

    assert(spm->loc2glob == NULL);

    (void)lda;
    return newspm;
}

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_dof
 *
 * @brief Expand a single dof CSR to a multi-dofs CSR.
 *
 * Each element matrix is fully initialized with the same element as the
 * original one.
 *
 *******************************************************************************
 *
 * @param[in] spm
 *           The original sparse CSR matrix used as the template for the
 *           muti-dof matrix.
 *
 *******************************************************************************
 *
 * @return The expanded CSR matrix according to the dofs properties previously
 * set.
 *
 *******************************************************************************/
spmatrix_t *
z_spmCSRExpand(const spmatrix_t *spm)
{
    spmatrix_t       *newspm;
    spm_int_t        i, j, k, ii, jj, dofi, dofj, col, row, baseval, lda;
    spm_int_t        diag, height;
    spm_int_t       *newcol, *newrow, *oldcol, *oldrow, *dofs;
    spm_complex64_t *newval = NULL;
    spm_complex64_t *oldval2, *oldval = NULL;

    assert( spm->fmttype == SpmCSR );
    assert( spm->flttype == SpmComplex64 );

    if ( spm->dof == 1 ) {
        return (spmatrix_t*)spm;
    }

    if ( spm->layout != SpmColMajor ) {
        fprintf( stderr, "Unsupported layout\n" );
        return NULL;
    }

    newspm = malloc( sizeof(spmatrix_t) );
    memcpy( newspm, spm, sizeof(spmatrix_t) );

    baseval = spmFindBase( spm );
    oldcol = spm->colptr;
    oldrow = spm->rowptr;
    dofs   = spm->dofs;
#if !defined(PRECISION_p)
    oldval = oldval2 = (spm_complex64_t*)(spm->values);
#endif

    newspm->n = spm->nexp;
    newspm->rowptr = newrow = malloc(sizeof(spm_int_t)*(spm->nexp+1));

    /**
     * First loop to compute the new rowptr
     */
    *newrow = baseval;
    for(i=0; i<spm->n; i++, oldrow++)
    {
        diag = 0;
        dofi = (spm->dof > 0 ) ? spm->dof : dofs[i+1] - dofs[i];

        /* Sum the width of the elements in the row */
        newrow[1] = newrow[0];
        for(k=oldrow[0]; k<oldrow[1]; k++)
        {
            j = oldcol[k-baseval] - baseval;
            dofj = (spm->dof > 0 ) ? spm->dof : dofs[j+1] - dofs[j];
            newrow[1] += dofj;

            diag = (diag || (i == j));
        }

        diag = (diag & (spm->mtxtype != SpmGeneral));
        height = newrow[1] - newrow[0];
        newrow++;

        /* Add extra rows */
        for(ii=1; ii<dofi; ii++, newrow++)
        {
            newrow[1] = newrow[0] + height;

            if ( diag ) {
                newrow[1] -= ii;
            }
        }
    }
    assert( ((spm->mtxtype == SpmGeneral) && ((newrow[0]-baseval) == spm->nnzexp)) ||
            ((spm->mtxtype != SpmGeneral) && ((newrow[0]-baseval) <= spm->nnzexp)) );

    newspm->nnz = newrow[0] - baseval;
    newspm->colptr = newcol = malloc(sizeof(spm_int_t)*newspm->nnz);
#if !defined(PRECISION_p)
    newspm->values = newval = malloc(sizeof(spm_complex64_t)*newspm->nnz);
#endif

    /**
     * Second loop to compute the new colptr and valptr
     */
    oldcol = spm->colptr;
    oldrow = spm->rowptr;
    newrow = newspm->rowptr;
    for(i=0, row=0; i<spm->n; i++, oldrow++)
    {
        /**
         * Backup current position in oldval because we will pick
         * interleaved data inside the buffer
         */
        lda = newrow[1] - newrow[0];
        oldval2 = oldval;

        if ( spm->dof > 0 ) {
            dofi = spm->dof;
            assert( row == spm->dof * i );
        }
        else {
            dofi = dofs[i+1] - dofs[i];
            assert( row == dofs[i] - baseval );
        }

        for(ii=0; ii<dofi; ii++, row++, newrow++)
        {
            assert( ((spm->mtxtype == SpmGeneral) && (lda == (newrow[1] - newrow[0]))) ||
                    ((spm->mtxtype != SpmGeneral) && (lda >= (newrow[1] - newrow[0]))) );

            /* Move to the beginning of the row ii in coming element (i,j) */
            oldval = oldval2 + ii;

            for(k=oldrow[0]; k<oldrow[1]; k++)
            {
                j = oldcol[k-baseval] - baseval;

                if ( spm->dof > 0 ) {
                    dofj = spm->dof;
                    col  = spm->dof * j;
                }
                else {
                    dofj = dofs[j+1] - dofs[j];
                    col  = dofs[j] - baseval;
                }

                for(jj=0; jj<dofj; jj++, col++)
                {
                    if ( (spm->mtxtype == SpmGeneral) ||
                         (i != j) ||
                         ((i == j) && (row <= col)) )
                    {
                        (*newcol) = col + baseval;
                        newcol++;
#if !defined(PRECISION_p)
                        (*newval) = *oldval;
                        newval++;
#endif
                    }
                    /* Move to next value in row ii */
                    oldval += dofi;
                }
            }
        }
        /* Move to the begining of the next row of elements */
        oldval -= (dofi-1);
    }

    newspm->gN      = newspm->n;
    newspm->gnnz    = newspm->nnz;

    newspm->gNexp   = newspm->gN;
    newspm->nexp    = newspm->n;
    newspm->gnnzexp = newspm->gnnz;
    newspm->nnzexp  = newspm->nnz;

    newspm->dof     = 1;
    newspm->dofs    = NULL;
    newspm->layout  = SpmColMajor;

    assert(spm->loc2glob == NULL);

    (void)newval;
    (void)lda;
    return newspm;
}

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_dof
 *
 * @brief Expand a single dof IJV to a multi-dofs IJV.
 *
 * Each element matrix is fully initialized with the same element as the
 * original one.
 *
 *******************************************************************************
 *
 * @param[in] spm
 *           The original sparse IJV matrix used as the template for the
 *           muti-dof matrix.
 *
 *******************************************************************************
 *
 * @return The expanded IJV matrix according to the dofs properties previously
 * set.
 *
 *******************************************************************************/
spmatrix_t *
z_spmIJVExpand(const spmatrix_t *spm)
{
    spmatrix_t       *newspm;
    spm_int_t        i, j, k, ii, jj, dofi, dofj, col, row, baseval;
    spm_int_t       *newcol, *newrow, *oldcol, *oldrow, *dofs;
#if !defined(PRECISION_p)
    spm_complex64_t *newval = NULL;
#endif
    spm_complex64_t *oldval = NULL;
    assert( spm->fmttype == SpmIJV );
    assert( spm->flttype == SpmComplex64 );

    if ( spm->dof == 1 ) {
        return (spmatrix_t*)spm;
    }

    newspm = malloc( sizeof(spmatrix_t) );
    memcpy( newspm, spm, sizeof(spmatrix_t) );

    baseval = spmFindBase( spm );
    oldcol = spm->colptr;
    oldrow = spm->rowptr;
    dofs   = spm->dofs;
#if !defined(PRECISION_p)
    oldval = (spm_complex64_t*)(spm->values);
#endif

    /**
     * First loop to compute the size of the vectores
     */
    newspm->n = spm->nexp;
    if (spm->mtxtype == SpmGeneral) {
        newspm->nnz = spm->nnzexp;
    }
    else {
        newspm->nnz = 0;
        for(k=0; k<spm->nnz; k++, oldrow++, oldcol++)
        {
            i = *oldrow - baseval;
            j = *oldcol - baseval;

            if ( spm->dof > 0 ) {
                dofi = spm->dof;
                dofj = spm->dof;
            }
            else {
                dofi = dofs[i+1] - dofs[i];
                dofj = dofs[j+1] - dofs[j];
            }

            if ( i != j ) {
                newspm->nnz += dofi * dofj;
            }
            else {
                assert( dofi == dofj );
                newspm->nnz += (dofi * (dofi+1)) / 2;
            }
        }
        assert( newspm->nnz <= spm->nnzexp );
    }

    newspm->rowptr = newrow = malloc(sizeof(spm_int_t)*newspm->nnz);
    newspm->colptr = newcol = malloc(sizeof(spm_int_t)*newspm->nnz);
#if !defined(PRECISION_p)
    newspm->values = newval = malloc(sizeof(spm_complex64_t)*newspm->nnz);
#endif

    /**
     * Second loop to compute the new rowptr, colptr and valptr
     */
    oldrow = spm->rowptr;
    oldcol = spm->colptr;
    for(k=0; k<spm->nnz; k++, oldrow++, oldcol++)
    {
        i = *oldrow - baseval;
        j = *oldcol - baseval;

        if ( spm->dof > 0 ) {
            dofi = spm->dof;
            row  = spm->dof * i;
            dofj = spm->dof;
            col  = spm->dof * j;
        }
        else {
            dofi = dofs[i+1] - dofs[i];
            row  = dofs[i] - baseval;
            dofj = dofs[j+1] - dofs[j];
            col  = dofs[j] - baseval;
        }

        if ( spm->layout == SpmColMajor ) {
            for(jj=0; jj<dofj; jj++)
            {
                for(ii=0; ii<dofi; ii++, oldval++)
                {
                    if ( (spm->mtxtype == SpmGeneral) ||
                         (i != j) ||
                         ((i == j) && (row+ii >= col+jj)) )
                    {
                        assert( row + ii < newspm->n );
                        assert( col + jj < newspm->n );
                        (*newrow) = row + ii + baseval;
                        (*newcol) = col + jj + baseval;
                        newrow++;
                        newcol++;
#if !defined(PRECISION_p)
                        (*newval) = *oldval;
                        newval++;
#endif
                    }
                }
            }
        }
        else {
            for(ii=0; ii<dofi; ii++)
            {
                for(jj=0; jj<dofj; jj++, oldval++)
                {
                    if ( (spm->mtxtype == SpmGeneral) ||
                         (i != j) ||
                         ((i == j) && (row+ii >= col+jj)) )
                    {
                        assert( row + ii < newspm->n );
                        assert( col + jj < newspm->n );
                        (*newrow) = row + ii + baseval;
                        (*newcol) = col + jj + baseval;
                        newrow++;
                        newcol++;
#if !defined(PRECISION_p)
                        (*newval) = *oldval;
                        newval++;
#endif
                    }
                }
            }
        }
    }
    assert( newcol - newspm->colptr == newspm->nnz );
    assert( newrow - newspm->rowptr == newspm->nnz );

    newspm->gN      = newspm->n;
    newspm->gnnz    = newspm->nnz;

    newspm->gNexp   = newspm->gN;
    newspm->nexp    = newspm->n;
    newspm->gnnzexp = newspm->gnnz;
    newspm->nnzexp  = newspm->nnz;

    newspm->dof     = 1;
    newspm->dofs    = NULL;
    newspm->layout  = SpmColMajor;

    assert(spm->loc2glob == NULL);

    return newspm;
}

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_dof
 *
 * @brief Expand a single dof sparse matrix to a multi-dofs sparse matrix.
 *
 * The original value of the element is replicated through the entire element
 * matrix that is generated in this routine.
 *
 *******************************************************************************
 *
 * @param[in] spm
 *           The original sparse IJV matrix used as the template for the
 *           muti-dof matrix.
 *
 *******************************************************************************
 *
 * @return The expanded matrix according to the dofs properties previously set.
 *
 *******************************************************************************/
spmatrix_t *
z_spmExpand( const spmatrix_t *spm )
{
    switch (spm->fmttype) {
    case SpmCSC:
        return z_spmCSCExpand( spm );
    case SpmCSR:
        return z_spmCSRExpand( spm );
    case SpmIJV:
        return z_spmIJVExpand( spm );
    }
    return NULL;
}
