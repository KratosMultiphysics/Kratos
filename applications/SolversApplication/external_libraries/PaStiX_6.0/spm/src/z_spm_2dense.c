/**
 *
 * @file z_spm_2dense.c
 *
 * SParse Matrix package conversion to dense routine.
 *
 * @copyright 2016-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Theophile Terraz
 * @author Alban Bellot
 * @date 2015-01-01
 *
 * @precisions normal z -> c d s
 *
 **/
#include "common.h"
#include "z_spm.h"

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_convert
 *
 * @brief Convert a CSC matrix into a dense matrix.
 *
 * The denses matrix is initialized with zeroes and filled with the spm matrix
 * values. When the matrix is hermitian or symmetric, both sides (upper and
 * lower) of the dense matrix are initialized.
 *
 *******************************************************************************
 *
 * @param[in] spm
 *          The sparse matrix in the CSC format.
 *
 *******************************************************************************
 *
 * @return A dense matrix in Lapack layout format
 *
 *******************************************************************************/
spm_complex64_t *
z_spmCSC2dense( const spmatrix_t *spm )
{
    spm_int_t i, j, k, lda, baseval;
    spm_complex64_t *A, *valptr;
    spm_int_t *colptr, *rowptr;

    assert( spm->fmttype == SpmCSC );
    assert( spm->flttype == SpmComplex64 );

    lda = spm->gNexp;
    A = (spm_complex64_t*)malloc(lda * lda * sizeof(spm_complex64_t));
    memset( A, 0, lda * lda * sizeof(spm_complex64_t));

    baseval = spmFindBase( spm );
    i = 0;
    j = 0;

    colptr = spm->colptr;
    rowptr = spm->rowptr;
    valptr = (spm_complex64_t*)(spm->values);

    /**
     * Constant degree of fredom of 1
     */
    if ( spm->dof == 1 ) {
        switch( spm->mtxtype ){
#if defined(PRECISION_z) || defined(PRECISION_c)
        case SpmHermitian:
            for(j=0; j<spm->n; j++, colptr++)
            {
                for(k=colptr[0]; k<colptr[1]; k++, rowptr++, valptr++)
                {
                    i = (*rowptr-baseval);
                    if( i == j ) {
                        /* Make sure the matrix is hermitian */
                        A[ j * lda + i ] = creal(*valptr) + I * 0.;
                    }
                    else {
                        A[ j * lda + i ] = *valptr;
                        A[ i * lda + j ] = conj(*valptr);
                    }
                }
            }
            break;
#endif
        case SpmSymmetric:
            for(j=0; j<spm->n; j++, colptr++)
            {
                for(k=colptr[0]; k<colptr[1]; k++, rowptr++, valptr++)
                {
                    i = (*rowptr-baseval);
                    A[ j * lda + i ] = *valptr;
                    A[ i * lda + j ] = *valptr;
                }
            }
            break;
        case SpmGeneral:
        default:
            for(j=0; j<spm->n; j++, colptr++)
            {
                for(k=colptr[0]; k<colptr[1]; k++, rowptr++, valptr++)
                {
                    i = (*rowptr-baseval);
                    A[ j * lda + i ] = *valptr;
                }
            }
        }
    }
    /**
     * General degree of freedom (constant or variable)
     */
    else {
        spm_int_t  k, ii, jj, dofi, dofj, col, row;
        spm_int_t *dofs = spm->dofs;

        switch( spm->mtxtype ){
#if defined(PRECISION_z) || defined(PRECISION_c)
        case SpmHermitian:
            for(j=0; j<spm->n; j++, colptr++)
            {
                dofj = ( spm->dof > 1 ) ?  spm->dof      : dofs[j+1] - dofs[j];
                col  = ( spm->dof > 1 ) ? (spm->dof * j) : dofs[j] - baseval;

                for(k=colptr[0]; k<colptr[1]; k++, rowptr++)
                {
                    i = (*rowptr - baseval);
                    dofi = ( spm->dof > 1 ) ?  spm->dof      : dofs[i+1] - dofs[i];
                    row  = ( spm->dof > 1 ) ? (spm->dof * i) : dofs[i] - baseval;

                    for(jj=0; jj<dofj; jj++)
                    {
                        for(ii=0; ii<dofi; ii++, valptr++)
                        {
                            if( col+jj == row+ii ) {
                                /* Make sure the matrix is hermitian */
                                A[ (col + jj) * lda + (row + ii) ] = creal(*valptr) + I * 0.;
                            }
                            else {
                                A[ (col + jj) * lda + (row + ii) ] = *valptr;
                                A[ (row + ii) * lda + (col + jj) ] = conj(*valptr);
                            }
                        }
                    }
                }
            }
            break;
#endif
        case SpmSymmetric:
            for(j=0; j<spm->n; j++, colptr++)
            {
                dofj = ( spm->dof > 1 ) ? spm->dof     : dofs[j+1] - dofs[j];
                col  = ( spm->dof > 1 ) ? spm->dof * j : dofs[j] - baseval;

                for(k=colptr[0]; k<colptr[1]; k++, rowptr++)
                {
                    i = (*rowptr - baseval);
                    dofi = ( spm->dof > 1 ) ? spm->dof     : dofs[i+1] - dofs[i];
                    row  = ( spm->dof > 1 ) ? spm->dof * i : dofs[i] - baseval;

                    for(jj=0; jj<dofj; jj++)
                    {
                        for(ii=0; ii<dofi; ii++, valptr++)
                        {
                            A[ (col + jj) * lda + (row + ii) ] = *valptr;
                            A[ (row + ii) * lda + (col + jj) ] = *valptr;
                        }
                    }
                }
            }
            break;
        case SpmGeneral:
        default:
            for(j=0; j<spm->n; j++, colptr++)
            {
                dofj = ( spm->dof > 1 ) ? spm->dof     : dofs[j+1] - dofs[j];
                col  = ( spm->dof > 1 ) ? spm->dof * j : dofs[j] - baseval;

                for(k=colptr[0]; k<colptr[1]; k++, rowptr++)
                {
                    i = (*rowptr - baseval);
                    dofi = ( spm->dof > 1 ) ? spm->dof     : dofs[i+1] - dofs[i];
                    row  = ( spm->dof > 1 ) ? spm->dof * i : dofs[i] - baseval;

                    for(jj=0; jj<dofj; jj++)
                    {
                        for(ii=0; ii<dofi; ii++, valptr++)
                        {
                            A[ (col + jj) * lda + (row + ii) ] = *valptr;
                        }
                    }
                }
            }
        }
    }
    return A;
}

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_convert
 *
 * @brief Convert a CSR matrix into a dense matrix.
 *
 * The denses matrix is initialized with zeroes and filled with the spm matrix
 * values. When the matrix is hermitian or symmetric, both sides (upper and
 * lower) of the dense matrix are initialized.
 *
 *******************************************************************************
 *
 * @param[in] spm
 *          The sparse matrix in the CSR format.
 *
 *******************************************************************************
 *
 * @return A dense matrix in Lapack layout format
 *
 *******************************************************************************/
spm_complex64_t *
z_spmCSR2dense( const spmatrix_t *spm )
{
    spm_int_t i, j, k, lda, baseval;
    spm_complex64_t *A, *valptr;
    spm_int_t *colptr, *rowptr;

    assert( spm->fmttype == SpmCSR );
    assert( spm->flttype == SpmComplex64 );

    lda = spm->gNexp;
    A = (spm_complex64_t*)malloc(lda * lda * sizeof(spm_complex64_t));
    memset( A, 0, lda * lda * sizeof(spm_complex64_t));

    baseval = spmFindBase( spm );
    i = 0;
    j = 0;

    colptr = spm->colptr;
    rowptr = spm->rowptr;
    valptr = (spm_complex64_t*)(spm->values);

    /**
     * Constant degree of fredom of 1
     */
    if ( spm->dof == 1 ) {
        switch( spm->mtxtype ){
#if defined(PRECISION_z) || defined(PRECISION_c)
        case SpmHermitian:
            for(i=0; i<spm->n; i++, rowptr++)
            {
                for(k=rowptr[0]; k<rowptr[1]; k++, colptr++, valptr++)
                {
                    j = (*colptr-baseval);
                    if( i == j ) {
                        /* Make sure the matrix is hermitian */
                        A[ j * lda + i ] = creal(*valptr) + I * 0.;
                    }
                    else {
                        A[ j * lda + i ] = *valptr;
                        A[ i * lda + j ] = conj(*valptr);
                    }
                }
            }
            break;
#endif
        case SpmSymmetric:
            for(i=0; i<spm->n; i++, rowptr++)
            {
                for(k=rowptr[0]; k<rowptr[1]; k++, colptr++, valptr++)
                {
                    j = (*colptr-baseval);
                    A[ j * lda + i ] = *valptr;
                    A[ i * lda + j ] = *valptr;
                }
            }
            break;
        case SpmGeneral:
        default:
            for(i=0; i<spm->n; i++, rowptr++)
            {
                for(k=rowptr[0]; k<rowptr[1]; k++, colptr++, valptr++)
                {
                    j = (*colptr-baseval);
                    A[ j * lda + i ] = *valptr;
                }
            }
        }
    }
    /**
     * General degree of freedom (constant or variable)
     */
    else {
        spm_int_t  k, ii, jj, dofi, dofj, col, row;
        spm_int_t *dofs = spm->dofs;

        switch( spm->mtxtype ){
#if defined(PRECISION_z) || defined(PRECISION_c)
        case SpmHermitian:
            for(i=0; i<spm->n; i++, rowptr++)
            {
                dofi = ( spm->dof > 1 ) ? spm->dof     : dofs[i+1] - dofs[i];
                row  = ( spm->dof > 1 ) ? spm->dof * i : dofs[i] - baseval;

                for(k=rowptr[0]; k<rowptr[1]; k++, colptr++)
                {
                    j = (*colptr - baseval);
                    dofj = ( spm->dof > 1 ) ? spm->dof     : dofs[j+1] - dofs[j];
                    col  = ( spm->dof > 1 ) ? spm->dof * j : dofs[j] - baseval;

                    for(jj=0; jj<dofj; jj++)
                    {
                        for(ii=0; ii<dofi; ii++, valptr++)
                        {
                            if( col+jj == row+ii ) {
                                /* Make sure the matrix is hermitian */
                                A[ (col + jj) * lda + (row + ii) ] = creal(*valptr) + I * 0.;
                            }
                            else {
                                A[ (col + jj) * lda + (row + ii) ] = *valptr;
                                A[ (row + ii) * lda + (col + jj) ] = conj(*valptr);
                            }
                        }
                    }
                }
            }
            break;
#endif
        case SpmSymmetric:
            for(i=0; i<spm->n; i++, rowptr++)
            {
                dofi = ( spm->dof > 1 ) ? spm->dof     : dofs[i+1] - dofs[i];
                row  = ( spm->dof > 1 ) ? spm->dof * i : dofs[i] - baseval;

                for(k=rowptr[0]; k<rowptr[1]; k++, colptr++)
                {
                    j = (*colptr - baseval);
                    dofj = ( spm->dof > 1 ) ? spm->dof     : dofs[j+1] - dofs[j];
                    col  = ( spm->dof > 1 ) ? spm->dof * j : dofs[j] - baseval;

                    for(jj=0; jj<dofj; jj++)
                    {
                        for(ii=0; ii<dofi; ii++, valptr++)
                        {
                            A[ (col + jj) * lda + (row + ii) ] = *valptr;
                            A[ (row + ii) * lda + (col + jj) ] = *valptr;
                        }
                    }
                }
            }
            break;
        case SpmGeneral:
        default:
            for(i=0; i<spm->n; i++, rowptr++)
            {
                dofi = ( spm->dof > 1 ) ? spm->dof     : dofs[i+1] - dofs[i];
                row  = ( spm->dof > 1 ) ? spm->dof * i : dofs[i] - baseval;

                for(k=rowptr[0]; k<rowptr[1]; k++, colptr++)
                {
                    j = (*colptr - baseval);
                    dofj = ( spm->dof > 1 ) ? spm->dof     : dofs[j+1] - dofs[j];
                    col  = ( spm->dof > 1 ) ? spm->dof * j : dofs[j] - baseval;

                    for(jj=0; jj<dofj; jj++)
                    {
                        for(ii=0; ii<dofi; ii++, valptr++)
                        {
                            A[ (col + jj) * lda + (row + ii) ] = *valptr;
                        }
                    }
                }
            }
        }
    }
    return A;
}

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_convert
 *
 * @brief Convert a IJV matrix into a dense matrix.
 *
 * The denses matrix is initialized with zeroes and filled with the spm matrix
 * values. When the matrix is hermitian or symmetric, both sides (upper and
 * lower) of the dense matrix are initialized.
 *
 *******************************************************************************
 *
 * @param[in] spm
 *          The sparse matrix in the IJV format.
 *
 *******************************************************************************
 *
 * @return A dense matrix in Lapack layout format
 *
 *******************************************************************************/
spm_complex64_t *
z_spmIJV2dense( const spmatrix_t *spm )
{
    spm_int_t i, j, k, lda, baseval;
    spm_complex64_t *A, *valptr;
    spm_int_t *colptr, *rowptr;

    assert( spm->fmttype == SpmIJV );
    assert( spm->flttype == SpmComplex64 );

    lda = spm->gNexp;
    A = (spm_complex64_t*)malloc(lda * lda * sizeof(spm_complex64_t));
    memset( A, 0, lda * lda * sizeof(spm_complex64_t));

    baseval = spmFindBase( spm );
    i = 0;
    j = 0;

    colptr = spm->colptr;
    rowptr = spm->rowptr;
    valptr = (spm_complex64_t*)(spm->values);

    /**
     * Constant degree of fredom of 1
     */
    if ( spm->dof == 1 ) {
        switch( spm->mtxtype ){
#if defined(PRECISION_z) || defined(PRECISION_c)
        case SpmHermitian:
            for(k=0; k<spm->nnz; k++, rowptr++, colptr++, valptr++)
            {
                i = *rowptr - baseval;
                j = *colptr - baseval;

                if( i == j ) {
                    /* Make sure the matrix is hermitian */
                    A[ i * lda + i ] = creal(*valptr) + I * 0.;
                }
                else {
                    A[ j * lda + i ] = *valptr;
                    A[ i * lda + j ] = conj(*valptr);
                }
            }
            break;
#endif
        case SpmSymmetric:
            for(k=0; k<spm->nnz; k++, rowptr++, colptr++, valptr++)
            {
                i = *rowptr - baseval;
                j = *colptr - baseval;

                A[ j * lda + i ] = *valptr;
                A[ i * lda + j ] = *valptr;
            }
            break;
        case SpmGeneral:
        default:
            for(k=0; k<spm->nnz; k++, rowptr++, colptr++, valptr++)
            {
                i = *rowptr - baseval;
                j = *colptr - baseval;

                A[ j * lda + i ] = *valptr;
            }
        }
    }
    /**
     * General degree of freedom (constant or variable)
     */
    else {
        spm_int_t  k, ii, jj, dofi, dofj, col, row;
        spm_int_t *dofs = spm->dofs;

        switch( spm->mtxtype ){
#if defined(PRECISION_z) || defined(PRECISION_c)
        case SpmHermitian:
            for(k=0; k<spm->nnz; k++, rowptr++, colptr++)
            {
                i = *rowptr - baseval;
                j = *colptr - baseval;

                if ( spm->dof > 1 ) {
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

                for(jj=0; jj<dofj; jj++)
                {
                    for(ii=0; ii<dofi; ii++, valptr++)
                    {
                        if( col+jj == row+ii ) {
                            /* Make sure the matrix is hermitian */
                            A[ (col+jj) * lda + (row+ii) ] = creal(*valptr) + I * 0.;
                        }
                        else {
                            A[ (col+jj) * lda + (row+ii) ] = *valptr;
                            A[ (row+ii) * lda + (col+jj) ] = conj(*valptr);
                        }
                    }
                }
            }
            break;
#endif
        case SpmSymmetric:
            for(k=0; k<spm->nnz; k++, rowptr++, colptr++)
            {
                i = *rowptr - baseval;
                j = *colptr - baseval;

                if ( spm->dof > 1 ) {
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

                for(jj=0; jj<dofj; jj++)
                {
                    for(ii=0; ii<dofi; ii++, valptr++)
                    {
                        A[ (col+jj) * lda + (row+ii) ] = *valptr;
                        A[ (row+ii) * lda + (col+jj) ] = *valptr;
                    }
                }
            }
            break;

        case SpmGeneral:
        default:
            for(k=0; k<spm->nnz; k++, rowptr++, colptr++)
            {
                i = *rowptr - baseval;
                j = *colptr - baseval;

                if ( spm->dof > 1 ) {
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

                for(jj=0; jj<dofj; jj++)
                {
                    for(ii=0; ii<dofi; ii++, valptr++)
                    {
                        A[ (col+jj) * lda + (row+ii) ] = *valptr;
                    }
                }
            }
        }
    }
    return A;
}

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_convert
 *
 * @brief Convert a sparse matrix into a dense matrix.
 *
 * The denses matrix is initialized with zeroes and filled with the spm matrix
 * values. When the matrix is hermitian or symmetric, both sides (upper and
 * lower) of the dense matrix are initialized.
 *
 *******************************************************************************
 *
 * @param[in] spm
 *          The sparse matrix to convert in any format.
 *
 *******************************************************************************
 *
 * @return A dense matrix in Lapack layout format
 *
 *******************************************************************************/
spm_complex64_t *
z_spm2dense( const spmatrix_t *spm )
{
    switch (spm->fmttype) {
    case SpmCSC:
        return z_spmCSC2dense( spm );
    case SpmCSR:
        return z_spmCSR2dense( spm );
    case SpmIJV:
        return z_spmIJV2dense( spm );
    }
    return NULL;
}

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_convert
 *
 * @brief Print a dense matrix to the given file
 *
 *******************************************************************************
 *
 * @param[in] f
 *          Open file descriptor on which to write the matrix.
 *
 * @param[in] m
 *          Number of rows of the matrix A.
 *
 * @param[in] n
 *          Number of columns of the matrix A.
 *
 *
 * @param[in] A
 *          The matrix to print of size lda -by- n
 *
 * @param[in] lda
 *          the leading dimension of the matrix A. lda >= m
 *
 *******************************************************************************/
void
z_spmDensePrint( FILE *f, spm_int_t m, spm_int_t n, const spm_complex64_t *A, spm_int_t lda )
{
    spm_int_t i, j;

    for(j=0; j<n; j++)
    {
        for(i=0; i<m; i++)
        {
            if ( cabs( A[ j * lda + i ] ) != 0. ) {
                z_spmPrintElt( f, i, j, A[lda * j + i] );
            }
        }
    }
    return;
}
