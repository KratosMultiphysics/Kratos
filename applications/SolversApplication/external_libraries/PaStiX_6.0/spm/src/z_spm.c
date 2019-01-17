/**
 *
 * @file z_spm.c
 *
 * SParse Matrix package precision dependent routines.
 *
 * @copyright 2016-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @date 2013-06-24
 *
 * @precisions normal z -> c d s p
 *
 **/
#include "common.h"
#include "z_spm.h"
#include <string.h>

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_check
 *
 * @brief This routine sorts the spm matrix.
 *
 * For the CSC and CSR formats, the subarray of edges for each vertex are sorted.
 * For the IJV format, the edges are storted first by column indexes, and then
 * by row indexes. To perform a sort first by row, second by column, please swap
 * the colptr and rowptr of the structure before calling the subroutine.
 *
 * @warning This function should NOT be called if dof is greater than 1.
 *
 *******************************************************************************
 *
 * @param[inout] spm
 *          On entry, the pointer to the sparse matrix structure.
 *          On exit, the same sparse matrix with subarrays of edges sorted by
 *          ascending order.
 *
 *******************************************************************************/
void
z_spmSort( spmatrix_t *spm )
{
    spm_int_t       *colptr = spm->colptr;
    spm_int_t       *rowptr = spm->rowptr;
    spm_complex64_t *values = spm->values;
    void *sortptr[3];
    spm_int_t n = spm->n;
    spm_int_t i, size;
    (void)sortptr;

#if !defined(PRECISION_p)
    assert( spm->dof == 1 );
#endif

    /* Sort in place each subset */
    if ( spm->fmttype == SpmCSC ) {
        for (i=0; i<n; i++, colptr++)
        {
            size = colptr[1] - colptr[0];

#if defined(PRECISION_p)
            spmIntSort1Asc1( rowptr, size );
#else
            sortptr[0] = rowptr;
            sortptr[1] = values;
            z_spmIntFltSortAsc( sortptr, size );
#endif
            rowptr += size;
            values += size;
        }
    }
    else if ( spm->fmttype == SpmCSR ) {
        for (i=0; i<n; i++, rowptr++)
        {
            size = rowptr[1] - rowptr[0];

#if defined(PRECISION_p)
            spmIntSort1Asc1( rowptr, size );
#else
            sortptr[0] = colptr;
            sortptr[1] = values;
            z_spmIntFltSortAsc( sortptr, size );
#endif
            colptr += size;
            values += size;
        }
    }
    else if ( spm->fmttype == SpmIJV ) {
        size = spm->nnz;

        sortptr[0] = colptr;
        sortptr[1] = rowptr;

#if defined(PRECISION_p)
        spmIntMSortIntAsc( sortptr, size );
#else
        sortptr[2] = values;
        z_spmIntIntFltSortAsc( sortptr, size );
#endif
    }
}

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_check
 *
 * @brief This routine merge the multiple entries in a sparse
 * matrix by suming their values together.
 *
 * The sparse matrix needs to be sorted first (see z_spmSort()).
 *
 * @warning This function should NOT be called if dof is greater than 1.
 *
 *******************************************************************************
 *
 * @param[inout] spm
 *          On entry, the pointer to the sparse matrix structure.
 *          On exit, the reduced sparse matrix of multiple entries were present
 *          in it. The multiple values for a same vertex are sum up together.
 *
 ********************************************************************************
 *
 * @return The number of vertices that were merged. -1 on error.
 *
 *******************************************************************************/
spm_int_t
z_spmMergeDuplicate( spmatrix_t *spm )
{
    spm_int_t       *colptr = spm->colptr;
    spm_int_t       *oldrow = spm->rowptr;
    spm_int_t       *newrow = spm->rowptr;
    spm_complex64_t *newval = spm->values;
    spm_complex64_t *oldval = spm->values;
    spm_int_t n       = spm->n;
    spm_int_t baseval = spm->colptr[0];
    spm_int_t dof2    = spm->dof * spm->dof;
    spm_int_t idx, i, j, size;
    spm_int_t merge = 0;
#if !defined(PRECISION_p)
    spm_int_t d;
#endif

    if ( spm->fmttype == SpmCSC ) {
        idx = 0;
        for (i=0; i<n; i++, colptr++)
        {
            size = colptr[1] - colptr[0];
            for (j = 0; j < size;
                 j++, oldrow++, oldval+=dof2, newrow++, newval+=dof2, idx++)
            {
                if ( newrow != oldrow ) {
                    newrow[0] = oldrow[0];
#if !defined(PRECISION_p)
                    memcpy( newval, oldval, dof2 * sizeof(spm_complex64_t) );
#endif
                }

                while( ((j+1) < size) && (newrow[0] == oldrow[1]) ) {
                    j++;
                    oldrow++;
                    oldval += dof2;
#if !defined(PRECISION_p)
                    /* Merge the two sets of values */
                    for (d = 0; d < dof2; d++) {
                        newval[d] += oldval[d];
                    }
#endif
                    merge++;
                }
            }
            assert( ((merge == 0) && ( colptr[1] == idx+baseval)) ||
                    ((merge != 0) && ( colptr[1] >  idx+baseval)) );

            colptr[1] = idx + baseval;
        }
        assert( ((merge == 0) && (spm->nnz         == idx)) ||
                ((merge != 0) && (spm->nnz - merge == idx)) );

        if (merge > 0) {
            spm->nnz = spm->nnz - merge;
            spm->gnnz = spm->nnz;

            newrow = malloc( spm->nnz * sizeof(spm_int_t) );
            memcpy( newrow, spm->rowptr, spm->nnz * sizeof(spm_int_t) );
            free(spm->rowptr);
            spm->rowptr = newrow;

#if !defined(PRECISION_p)
            newval = malloc( spm->nnz * dof2 * sizeof(spm_int_t) );
            memcpy( newval, spm->values, spm->nnz * dof2 * sizeof(spm_complex64_t) );
            free(spm->values);
            spm->values = newval;
#endif
        }
    }

    return merge;
}

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_check
 *
 * @brief This routine corrects the sparse matrix structure if it's
 * pattern is not symmetric.
 *
 * It returns the new symmetric pattern with zeroes on
 * the new entries.
 *
 *******************************************************************************
 *
 * @param[inout] spm
 *          On entry, the pointer to the sparse matrix structure.
 *          On exit, the same sparse matrix with extra entries that makes it
 *          pattern symmetric.
 *
 *******************************************************************************
 *
 * @retval Return the number of entries added to the matrix.
 *
 *******************************************************************************/
spm_int_t
z_spmSymmetrize( spmatrix_t *spm )
{
    spm_complex64_t *oldval, *valtmp, *newval = NULL;
    spm_int_t *oldcol, *coltmp, *newcol = NULL;
    spm_int_t *oldrow, *rowtmp, *newrow = NULL;
    spm_int_t *toaddtab = NULL;
    spm_int_t *toaddtmp, toaddcnt, toaddsze;
    spm_int_t  n       = spm->n;
    spm_int_t  dof2    = spm->dof * spm->dof;
    spm_int_t  i, j, k, r, size;
    spm_int_t  baseval;

    toaddcnt = 0;
    toaddsze = 0;

    if ( (spm->fmttype == SpmCSC) || (spm->fmttype == SpmCSR) ) {
        if (spm->fmttype == SpmCSC) {
            oldcol = spm->colptr;
            coltmp = spm->colptr;
            oldrow = spm->rowptr;
            rowtmp = spm->rowptr;
        }
        else {
            oldcol = spm->rowptr;
            coltmp = spm->rowptr;
            oldrow = spm->colptr;
            rowtmp = spm->colptr;
        }

        baseval  = oldcol[0];
        for (i=0; i<n; i++, coltmp++)
        {
            size = coltmp[1] - coltmp[0];
            for (r=0; r<size; r++, rowtmp++ )
            {
                j = rowtmp[0]-baseval;
                if ( i != j ) {
                    /* Look for the element (j, i) */
                    spm_int_t frow = oldcol[ j ] - baseval;
                    spm_int_t lrow = oldcol[ j+1 ] - baseval;
                    int found = 0;

                    for (k = frow; k < lrow; k++)
                    {
                        if (i == (oldrow[k]-baseval))
                        {
                            found = 1;
                            break;
                        }
                        else if ( i < (oldrow[k]-baseval))
                        {
                            /* The spm is sorted so we will not find it later */
                            break;
                        }
                    }

                    if ( !found ) {
                        if ( newcol == NULL ) {
                            newcol = malloc( (spm->n+1) * sizeof(spm_int_t) );
                            for (k=0; k<spm->n; k++) {
                                newcol[k] = oldcol[k+1] - oldcol[k];
                            }

                            /* Let's start with a buffer that will contain 5% of extra elements */
                            toaddsze = spm_imax( 0.05 * (double)spm->nnz, 1 );
                            toaddtab = malloc( 2*toaddsze * sizeof(spm_int_t) );
                        }

                        if (toaddcnt >= toaddsze) {
                            toaddsze *= 2;
                            toaddtab = (spm_int_t*)realloc( toaddtab, 2*toaddsze*sizeof(spm_int_t) );
                        }

                        /* Newcol stores the number of elements per column */
                        newcol[ j ]++;
                        /* toaddtab stores the couples (j, i) to be added in the final sparse matrix */
                        toaddtab[ toaddcnt * 2     ] = j;
                        toaddtab[ toaddcnt * 2 + 1 ] = i;
                        toaddcnt++;
                    }
                }
            }
        }

        if (toaddcnt > 0) {
            spm_int_t newsize, oldsize;

            /* Sort the array per column */
            spmIntSort2Asc1(toaddtab, toaddcnt);

            spm->nnz = spm->nnz + toaddcnt;
            spm->gnnz = spm->nnz;

            newrow = malloc( spm->nnz        * sizeof(spm_int_t) );
#if !defined(PRECISION_p)
            newval = malloc( spm->nnz * dof2 * sizeof(spm_complex64_t) );
#endif
            coltmp = newcol;
            rowtmp = newrow;
            valtmp = newval;
            oldval = spm->values;
            toaddtmp = toaddtab;

            newsize = coltmp[0];
            coltmp[0] = baseval;
            for (i=0; i<n; i++, coltmp++, oldcol++)
            {
                /* Copy the existing elements */
                oldsize = oldcol[1] - oldcol[0];
                memcpy( rowtmp, oldrow, oldsize * sizeof(spm_int_t) );
#if !defined(PRECISION_p)
                memcpy( valtmp, oldval, oldsize * dof2 * sizeof(spm_complex64_t) );
#endif

                oldrow += oldsize;
                rowtmp += oldsize;
                oldval += oldsize * dof2;
                valtmp += oldsize * dof2;

                /* Some elements have been added */
                assert( newsize >= oldsize );
                if ( newsize > oldsize ) {
                    int added = 0;
                    int tobeadded = newsize - oldsize;

                    /* At least one element is equal to i */
                    assert( toaddtmp[0] == i );

                    /* Let's add the new vertices */
                    while( (added < tobeadded) && (toaddtmp[0] == i) )
                    {
                        rowtmp[0] = toaddtmp[1] + baseval;
                        rowtmp++;
                        toaddtmp+=2;
                        added++;
                    }
                    assert( added == tobeadded );

#if !defined(PRECISION_p)
                    /* Add 0.s for the associated values */
                    memset( valtmp, 0, added * dof2 * sizeof(spm_complex64_t) );
                    valtmp += added * dof2;
#endif
                }

                /* Use oldsize as temporary variable to update the new colptr */
                oldsize = newsize;
                newsize = coltmp[1];
                coltmp[1] = coltmp[0] + oldsize;
            }

            assert( coltmp[0]-baseval == spm->nnz );
            free( toaddtab );
            free( spm->colptr );
            free( spm->rowptr );
            free( spm->values );
            if (spm->fmttype == SpmCSC) {
                spm->colptr = newcol;
                spm->rowptr = newrow;
            }
            else {
                spm->colptr = newrow;
                spm->rowptr = newcol;
            }
            spm->values = newval;
        }
    }

    return toaddcnt;
}
