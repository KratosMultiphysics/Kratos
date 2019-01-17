/**
 *
 * @file z_spm_print.c
 *
 * SParse Matrix package printing routines.
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
 * @precisions normal z -> c d s p
 *
 **/
#include "common.h"
#include "z_spm.h"

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_print
 *
 * @brief Write CSC matrix in a file
 *
 *******************************************************************************
 *
 * @param[in] f
 *          Output file
 *
 * @param[in] spm
 *          The spm structure describing the matrix.
 *
 *******************************************************************************/
void
z_spmCSCPrint( FILE *f, const spmatrix_t *spm )
{
    spm_int_t i, j, baseval;
    spm_int_t k, ii, jj, dofi, dofj, col, row;
    spm_complex64_t *valptr;
    spm_int_t *colptr, *rowptr, *dofs;

    assert( spm->fmttype == SpmCSC );
    assert( spm->flttype == SpmComplex64 );

    baseval = spmFindBase( spm );
    i = 0;
    j = 0;

    colptr = spm->colptr;
    rowptr = spm->rowptr;
    valptr = (spm_complex64_t*)(spm->values);
    dofs   = spm->dofs;

    switch( spm->mtxtype ){
#if defined(PRECISION_z) || defined(PRECISION_c)
    case SpmHermitian:
        for(j=0; j<spm->n; j++, colptr++)
        {
            dofj = ( spm->dof > 0 ) ? spm->dof     : dofs[j+1] - dofs[j];
            col  = ( spm->dof > 0 ) ? spm->dof * j : dofs[j] - baseval;

            for(k=colptr[0]; k<colptr[1]; k++, rowptr++)
            {
                i = (*rowptr - baseval);
                dofi = ( spm->dof > 0 ) ? spm->dof     : dofs[i+1] - dofs[i];
                row  = ( spm->dof > 0 ) ? spm->dof * i : dofs[i] - baseval;

                if ( spm->layout == SpmColMajor ) {
                    for(jj=0; jj<dofj; jj++)
                    {
                        for(ii=0; ii<dofi; ii++, valptr++)
                        {
                            if ( row == col ) {
                                if (row+ii >= col+jj) {
                                    z_spmPrintElt( f, row + ii, col + jj, *valptr );
                                    if (row+ii > col+jj) {
                                        z_spmPrintElt( f, col + jj, row + ii, conj(*valptr) );
                                    }
                                }
                            }
                            else {
                                z_spmPrintElt( f, row + ii, col + jj, *valptr );
                                z_spmPrintElt( f, col + jj, row + ii, conj(*valptr) );
                            }
                        }
                    }
                }
                else {
                    for(ii=0; ii<dofi; ii++)
                    {
                        for(jj=0; jj<dofj; jj++, valptr++)
                        {
                            if ( row == col ) {
                                if (row+ii >= col+jj) {
                                    z_spmPrintElt( f, row + ii, col + jj, *valptr );
                                    if (row+ii > col+jj) {
                                        z_spmPrintElt( f, col + jj, row + ii, conj(*valptr) );
                                    }
                                }
                            }
                            else {
                                z_spmPrintElt( f, row + ii, col + jj, *valptr );
                                z_spmPrintElt( f, col + jj, row + ii, conj(*valptr) );
                            }
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
            dofj = ( spm->dof > 0 ) ? spm->dof     : dofs[j+1] - dofs[j];
            col  = ( spm->dof > 0 ) ? spm->dof * j : dofs[j] - baseval;

            for(k=colptr[0]; k<colptr[1]; k++, rowptr++)
            {
                i = (*rowptr - baseval);
                dofi = ( spm->dof > 0 ) ? spm->dof     : dofs[i+1] - dofs[i];
                row  = ( spm->dof > 0 ) ? spm->dof * i : dofs[i] - baseval;

                if ( spm->layout == SpmColMajor ) {
                    for(jj=0; jj<dofj; jj++)
                    {
                        for(ii=0; ii<dofi; ii++, valptr++)
                        {
                            if ( row == col ) {
                                if (row+ii >= col+jj) {
                                    z_spmPrintElt( f, row + ii, col + jj, *valptr );
                                    if (row+ii > col+jj) {
                                        z_spmPrintElt( f, col + jj, row + ii, *valptr );
                                    }
                                }
                            }
                            else {
                                z_spmPrintElt( f, row + ii, col + jj, *valptr );
                                z_spmPrintElt( f, col + jj, row + ii, *valptr );
                            }
                        }
                    }
                }
                else {
                    for(ii=0; ii<dofi; ii++)
                    {
                        for(jj=0; jj<dofj; jj++, valptr++)
                        {
                            if ( row == col ) {
                                if (row+ii >= col+jj) {
                                    z_spmPrintElt( f, row + ii, col + jj, *valptr );
                                    if (row+ii > col+jj) {
                                        z_spmPrintElt( f, col + jj, row + ii, *valptr );
                                    }
                                }
                            }
                            else {
                                z_spmPrintElt( f, row + ii, col + jj, *valptr );
                                z_spmPrintElt( f, col + jj, row + ii, *valptr );
                            }
                        }
                    }
                }
            }
        }
        break;
    case SpmGeneral:
    default:
        for(j=0; j<spm->n; j++, colptr++)
        {
            dofj = ( spm->dof > 0 ) ? spm->dof     : dofs[j+1] - dofs[j];
            col  = ( spm->dof > 0 ) ? spm->dof * j : dofs[j] - baseval;

            for(k=colptr[0]; k<colptr[1]; k++, rowptr++)
            {
                i = (*rowptr - baseval);
                dofi = ( spm->dof > 0 ) ? spm->dof     : dofs[i+1] - dofs[i];
                row  = ( spm->dof > 0 ) ? spm->dof * i : dofs[i] - baseval;

                if ( spm->layout == SpmColMajor ) {
                    for(jj=0; jj<dofj; jj++)
                    {
                        for(ii=0; ii<dofi; ii++, valptr++)
                        {
                            z_spmPrintElt( f, row + ii, col + jj, *valptr );
                        }
                    }
                }
                else {
                    for(ii=0; ii<dofi; ii++)
                    {
                        for(jj=0; jj<dofj; jj++, valptr++)
                        {
                            z_spmPrintElt( f, row + ii, col + jj, *valptr );
                        }
                    }
                }
            }
        }
    }
    return;
}

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_print
 *
 * @brief Write CSR matrix in a file
 *
 *******************************************************************************
 *
 * @param[in] f
 *          Output file
 *
 * @param[in] spm
 *          The spm structure describing the matrix.
 *
 *******************************************************************************/
void
z_spmCSRPrint( FILE *f, const spmatrix_t *spm )
{
    spm_int_t i, j, baseval;
    spm_int_t k, ii, jj, dofi, dofj, col, row;
    spm_complex64_t *valptr;
    spm_int_t *colptr, *rowptr, *dofs;

    assert( spm->fmttype == SpmCSR );
    assert( spm->flttype == SpmComplex64 );

    baseval = spmFindBase( spm );
    i = 0;
    j = 0;

    colptr = spm->colptr;
    rowptr = spm->rowptr;
    valptr = (spm_complex64_t*)(spm->values);
    dofs   = spm->dofs;

    switch( spm->mtxtype ){
#if defined(PRECISION_z) || defined(PRECISION_c)
    case SpmHermitian:
        for(i=0; i<spm->n; i++, rowptr++)
        {
            dofi = ( spm->dof > 0 ) ? spm->dof     : dofs[i+1] - dofs[i];
            row  = ( spm->dof > 0 ) ? spm->dof * i : dofs[i] - baseval;

            for(k=rowptr[0]; k<rowptr[1]; k++, colptr++)
            {
                j = (*colptr - baseval);
                dofj = ( spm->dof > 0 ) ? spm->dof     : dofs[j+1] - dofs[j];
                col  = ( spm->dof > 0 ) ? spm->dof * j : dofs[j] - baseval;

                if ( spm->layout == SpmColMajor ) {
                    for(jj=0; jj<dofj; jj++)
                    {
                        for(ii=0; ii<dofi; ii++, valptr++)
                        {
                            if ( row == col ) {
                                if (row+ii >= col+jj) {
                                    z_spmPrintElt( f, row + ii, col + jj, *valptr );
                                    if (row+ii > col+jj) {
                                        z_spmPrintElt( f, col + jj, row + ii, conj(*valptr) );
                                    }
                                }
                            }
                            else {
                                z_spmPrintElt( f, row + ii, col + jj, *valptr );
                                z_spmPrintElt( f, col + jj, row + ii, conj(*valptr) );
                            }
                        }
                    }
                }
                else {
                    for(ii=0; ii<dofi; ii++)
                    {
                        for(jj=0; jj<dofj; jj++, valptr++)
                        {
                            if ( row == col ) {
                                if (row+ii >= col+jj) {
                                    z_spmPrintElt( f, row + ii, col + jj, *valptr );
                                    if (row+ii > col+jj) {
                                        z_spmPrintElt( f, col + jj, row + ii, conj(*valptr) );
                                    }
                                }
                            }
                            else {
                                z_spmPrintElt( f, row + ii, col + jj, *valptr );
                                z_spmPrintElt( f, col + jj, row + ii, conj(*valptr) );
                            }
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
            dofi = ( spm->dof > 0 ) ? spm->dof     : dofs[i+1] - dofs[i];
            row  = ( spm->dof > 0 ) ? spm->dof * i : dofs[i] - baseval;

            for(k=rowptr[0]; k<rowptr[1]; k++, colptr++)
            {
                j = (*colptr - baseval);
                dofj = ( spm->dof > 0 ) ? spm->dof     : dofs[j+1] - dofs[j];
                col  = ( spm->dof > 0 ) ? spm->dof * j : dofs[j] - baseval;

                if ( spm->layout == SpmColMajor ) {
                    for(jj=0; jj<dofj; jj++)
                    {
                        for(ii=0; ii<dofi; ii++, valptr++)
                        {
                            if ( row == col ) {
                                if (row+ii >= col+jj) {
                                    z_spmPrintElt( f, row + ii, col + jj, *valptr );
                                    if (row+ii > col+jj) {
                                        z_spmPrintElt( f, col + jj, row + ii, *valptr );
                                    }
                                }
                            }
                            else {
                                z_spmPrintElt( f, row + ii, col + jj, *valptr );
                                z_spmPrintElt( f, col + jj, row + ii, *valptr );
                            }
                        }
                    }
                }
                else {
                    for(ii=0; ii<dofi; ii++)
                    {
                        for(jj=0; jj<dofj; jj++, valptr++)
                        {
                            if ( row == col ) {
                                if (row+ii >= col+jj) {
                                    z_spmPrintElt( f, row + ii, col + jj, *valptr );
                                    if (row+ii > col+jj) {
                                        z_spmPrintElt( f, col + jj, row + ii, *valptr );
                                    }
                                }
                            }
                            else {
                                z_spmPrintElt( f, row + ii, col + jj, *valptr );
                                z_spmPrintElt( f, col + jj, row + ii, *valptr );
                            }
                        }
                    }
                }
            }
        }
        break;
    case SpmGeneral:
    default:
        for(i=0; i<spm->n; i++, rowptr++)
        {
            dofi = ( spm->dof > 0 ) ? spm->dof     : dofs[i+1] - dofs[i];
            row  = ( spm->dof > 0 ) ? spm->dof * i : dofs[i] - baseval;

            for(k=rowptr[0]; k<rowptr[1]; k++, colptr++)
            {
                j = (*colptr - baseval);
                dofj = ( spm->dof > 0 ) ? spm->dof     : dofs[j+1] - dofs[j];
                col  = ( spm->dof > 0 ) ? spm->dof * j : dofs[j] - baseval;

                if ( spm->layout == SpmColMajor ) {
                    for(jj=0; jj<dofj; jj++)
                    {
                        for(ii=0; ii<dofi; ii++, valptr++)
                        {
                            z_spmPrintElt( f, row + ii, col + jj, *valptr );
                        }
                    }
                }
                else {
                    for(ii=0; ii<dofi; ii++)
                    {
                        for(jj=0; jj<dofj; jj++, valptr++)
                        {
                            z_spmPrintElt( f, row + ii, col + jj, *valptr );
                        }
                    }
                }
            }
        }
    }
    return;
}

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_print
 *
 * @brief Write IJV matrix in a file
 *
 *******************************************************************************
 *
 * @param[in] f
 *          Output file
 *
 * @param[in] spm
 *          The spm structure describing the matrix.
 *
 *******************************************************************************/
void
z_spmIJVPrint( FILE *f, const spmatrix_t *spm )
{
    spm_int_t i, j, baseval;
    spm_int_t k, ii, jj, dofi, dofj, col, row;
    spm_complex64_t *valptr;
    spm_int_t *colptr, *rowptr, *dofs;

    assert( spm->fmttype == SpmIJV );
    assert( spm->flttype == SpmComplex64 );

    baseval = spmFindBase( spm );
    i = 0;
    j = 0;

    colptr = spm->colptr;
    rowptr = spm->rowptr;
    valptr = (spm_complex64_t*)(spm->values);
    dofs   = spm->dofs;

    switch( spm->mtxtype ){
#if defined(PRECISION_z) || defined(PRECISION_c)
    case SpmHermitian:
        for(k=0; k<spm->nnz; k++, rowptr++, colptr++)
        {
            i = *rowptr - baseval;
            j = *colptr - baseval;

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
                    for(ii=0; ii<dofi; ii++, valptr++)
                    {
                        if ( row == col ) {
                            if (row+ii >= col+jj) {
                                z_spmPrintElt( f, row + ii, col + jj, *valptr );
                                if (row+ii > col+jj) {
                                    z_spmPrintElt( f, col + jj, row + ii, conj(*valptr) );
                                }
                            }
                        }
                        else {
                            z_spmPrintElt( f, row + ii, col + jj, *valptr );
                            z_spmPrintElt( f, col + jj, row + ii, conj(*valptr) );
                        }
                    }
                }
            }
            else {
                for(ii=0; ii<dofi; ii++)
                {
                    for(jj=0; jj<dofj; jj++, valptr++)
                    {
                        if ( row == col ) {
                            if (row+ii >= col+jj) {
                                z_spmPrintElt( f, row + ii, col + jj, *valptr );
                                if (row+ii > col+jj) {
                                    z_spmPrintElt( f, col + jj, row + ii, conj(*valptr) );
                                }
                            }
                        }
                        else {
                            z_spmPrintElt( f, row + ii, col + jj, *valptr );
                            z_spmPrintElt( f, col + jj, row + ii, conj(*valptr) );
                        }
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
                    for(ii=0; ii<dofi; ii++, valptr++)
                    {
                        if ( row == col ) {
                            if (row+ii >= col+jj) {
                                z_spmPrintElt( f, row + ii, col + jj, *valptr );
                                if (row+ii > col+jj) {
                                    z_spmPrintElt( f, col + jj, row + ii, *valptr );
                                }
                            }
                        }
                        else {
                            z_spmPrintElt( f, row + ii, col + jj, *valptr );
                            z_spmPrintElt( f, col + jj, row + ii, *valptr );
                        }
                    }
                }
            }
            else {
                for(ii=0; ii<dofi; ii++)
                {
                    for(jj=0; jj<dofj; jj++, valptr++)
                    {
                        if ( row == col ) {
                            if (row+ii >= col+jj) {
                                z_spmPrintElt( f, row + ii, col + jj, *valptr );
                                if (row+ii > col+jj) {
                                    z_spmPrintElt( f, col + jj, row + ii, *valptr );
                                }
                            }
                        }
                        else {
                            z_spmPrintElt( f, row + ii, col + jj, *valptr );
                            z_spmPrintElt( f, col + jj, row + ii, *valptr );
                        }
                    }
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
                    for(ii=0; ii<dofi; ii++, valptr++)
                    {
                        z_spmPrintElt( f, row + ii, col + jj, *valptr );
                    }
                }
            }
            else {
                for(ii=0; ii<dofi; ii++)
                {
                    for(jj=0; jj<dofj; jj++, valptr++)
                    {
                        z_spmPrintElt( f, row + ii, col + jj, *valptr );
                    }
                }
            }
        }
    }
    return;
}

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_print
 *
 * @brief Write a spm matrix in a file
 *
 *******************************************************************************
 *
 * @param[in] f
 *          Output file
 *
 * @param[in] spm
 *          The spm structure describing the matrix.
 *
 *******************************************************************************/
void
z_spmPrint( FILE *f, const spmatrix_t *spm )
{
    switch (spm->fmttype) {
    case SpmCSC:
        z_spmCSCPrint( f, spm );
        break;
    case SpmCSR:
        z_spmCSRPrint( f, spm );
        break;
    case SpmIJV:
        z_spmIJVPrint( f, spm );
    }
    return;
}
