/**
 * @file spm_gen_fake_values.c
 *
 * SParse Matrix generic laplacian value generator routines.
 *
 * @copyright 2016-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Xavier Lacoste
 * @author Theophile Terraz
 * @date 2015-01-01
 *
 **/
#include "common.h"

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_driver
 *
 * @brief Compute the degree of each vertex.
 *
 *******************************************************************************
 *
 * @param[in] spm
 *          The spm to study.
 *
 * @param[inout] degrees
 *          Array of size spm->n allocated on entry. On exit, contains the
 *          degree of each vertex in the spm matrix.
 *
 *******************************************************************************
 *
 * @return the number of diagonal elements found during the computation.
 *
 *******************************************************************************/
static inline spm_int_t
spm_compute_degrees( const spmatrix_t *spm,
                     spm_int_t *degrees )
{
    spm_int_t i, j, k;
    spm_int_t *colptr = spm->colptr;
    spm_int_t *rowptr = spm->rowptr;
    spm_int_t baseval;
    spm_int_t diagval = 0;

    baseval = spmFindBase( spm );
    memset( degrees, 0, spm->n * sizeof(spm_int_t) );

    switch(spm->fmttype)
    {
    case SpmCSR:
        /* Swap pointers to call CSC */
        colptr = spm->rowptr;
        rowptr = spm->colptr;

        spm_attr_fallthrough;

    case SpmCSC:
        for(j=0; j<spm->n; j++, colptr++) {
            for(k=colptr[0]; k<colptr[1]; k++, rowptr++) {
                i = *rowptr - baseval;

                if ( i != j ) {
                    degrees[j] += 1;
                    if ( spm->mtxtype != SpmGeneral ) {
                        degrees[i] += 1;
                    }
                }
                else {
                    diagval++;
                }
            }
        }
        break;
    case SpmIJV:
        for(k=0; k<spm->nnz; k++, rowptr++, colptr++)
        {
            i = *rowptr - baseval;
            j = *colptr - baseval;

            if ( i != j ) {
                degrees[j] += 1;

                if ( spm->mtxtype != SpmGeneral ) {
                    degrees[i] += 1;
                }
            }
            else {
                diagval++;
            }
        }
    }

    return diagval;
}

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_driver
 *
 * @brief Insert diagonal elements to the graph to have a full Laplacian
 * generated
 *
 *******************************************************************************
 *
 * @param[in] spm
 *          At start, the initial spm structure with missing diagonal elements.
 *          At exit, contains the same sparse matrix with diagonal elements added.
 *
 * @param[in] diagval
 *          The number of diagonal elements already present in the matrix.
 *
 *******************************************************************************/
static inline void
spm_add_diag( spmatrix_t *spm,
              spm_int_t  diagval )
{
    spmatrix_t oldspm;
    spm_int_t i, j, k;
    spm_int_t *oldcol = spm->colptr;
    spm_int_t *oldrow = spm->rowptr;
    spm_int_t *newrow, *newcol;
    spm_int_t baseval;

    baseval = spmFindBase( spm );

    memcpy( &oldspm, spm, sizeof(spmatrix_t));

    spm->nnz = oldspm.nnz + (spm->n - diagval);
    newrow = malloc( spm->nnz * sizeof(spm_int_t) );

    switch(spm->fmttype)
    {
    case SpmCSR:
        /* Swap pointers to call CSC */
        oldcol = spm->rowptr;
        oldrow = spm->colptr;
        spm->colptr = newrow;

        spm_attr_fallthrough;

    case SpmCSC:
        newcol = oldcol;
        if ( spm->fmttype == SpmCSC ) {
            spm->rowptr = newrow;
        }
        diagval = 0;
        for(j=0; j<spm->n; j++, newcol++) {
            spm_int_t nbelt = newcol[1] - newcol[0];
            int diag = 0;

            memcpy( newrow, oldrow, nbelt * sizeof(spm_int_t) );
            newrow += nbelt;

            for(k=0; k<nbelt; k++, oldrow++) {
                i = *oldrow - baseval;

                if ( i == j ) {
                    diag = 1;
                }
            }
            newcol[0] += diagval;
            if ( !diag ) {
                *newrow = j + baseval;
                newrow++;
                diagval++;
            }
        }
        newcol[0] += diagval;

        if ( spm->fmttype == SpmCSC ) {
            free( oldspm.rowptr );
        }
        else {
            free( oldspm.colptr );
        }
        assert( diagval == spm->n );
        break;

    case SpmIJV:
        newcol = malloc( spm->nnz * sizeof(spm_int_t) );
        spm->colptr = newcol;
        spm->rowptr = newrow;

        for(k=0; k<spm->n; k++, newrow++, newcol++)
        {
            *newrow = k + baseval;
            *newcol = k + baseval;
        }

        for(k=0; k<spm->nnz; k++, oldrow++, oldcol++)
        {
            i = *oldrow - baseval;
            j = *oldcol - baseval;

            if ( i == j ) {
                continue;
            }

            *newrow = i + baseval;
            *newcol = j + baseval;
            newrow++;
            newcol++;
        }
        free( oldspm.colptr );
        free( oldspm.rowptr );
    }
}

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_driver
 *
 * @brief Generate the fake values array such that \[ M =  \alpha * D - \beta * A \]
 *
 * D is the degree matrix, and A the adjacency matrix.
 *
 *******************************************************************************
 *
 * @param[inout] spm
 *          The spm structure for which the values array must be generated.
 *
 * @param[in] degrees
 *          Array of size spm->n that containes the degree of each vertex in the
 *          spm structure.
 *
 *******************************************************************************/
static inline void
spm_generate_fake_values( spmatrix_t *spm,
                          const spm_int_t *degrees,
                          double alpha, double beta )
{
    double *values;
    spm_int_t i, j, k;
    spm_int_t *colptr = spm->colptr;
    spm_int_t *rowptr = spm->rowptr;
    spm_int_t baseval;

    baseval = spmFindBase( spm );

    spm->values = malloc( spm->nnzexp * sizeof(double) );
    values = spm->values;

    switch(spm->fmttype)
    {
    case SpmCSR:
        /* Swap pointers to call CSC */
        colptr = spm->rowptr;
        rowptr = spm->colptr;

        spm_attr_fallthrough;

    case SpmCSC:
        for(j=0; j<spm->n; j++, colptr++) {
            for(k=colptr[0]; k<colptr[1]; k++, rowptr++, values++) {
                i = *rowptr - baseval;

                if ( i == j ) {
                    *values = alpha * degrees[j];
                }
                else {
                    *values = - beta;
                }
            }
        }
        break;
    case SpmIJV:
        for(k=0; k<spm->nnz; k++, rowptr++, colptr++, values++)
        {
            i = *rowptr - baseval;
            j = *colptr - baseval;

            if ( i == j ) {
                *values = alpha * degrees[j];
            }
            else {
                *values = - beta;
            }
        }
    }

    spm->flttype = SpmDouble;
    if ( spm->mtxtype == SpmHermitian ) {
        spm->mtxtype = SpmSymmetric;
    }
}

/**
 *******************************************************************************
 *
 * @ingroup spm
 *
 * @brief Generate the fake values array such that \[ M =  \alpha * D - \beta * A \]
 *
 * D is the degree matrix, and A the adjacency matrix. The resulting matrix uses
 * real double.
 *
 *******************************************************************************
 *
 * @param[inout] spm
 *          The spm structure for which the values array must be generated.
 *
 *******************************************************************************/
void
spmGenFakeValues( spmatrix_t *spm )
{
    spm_int_t *degrees, diagval;
    double alpha = 10.;
    double beta  = 1.;

    assert( spm->flttype == SpmPattern );
    assert( spm->values == NULL );
    assert( spm->dof == 1 );

    /*
     * Read environment values for alpha/beta
     */
    {
        char *str = spm_getenv( "SPM_FAKE_ALPHA" );
        double value;

        if ( str != NULL ) {
            value = strtod( str, NULL );
            if ( (value != HUGE_VAL) && (value != 0.) &&
                 !isnan(value) && !isinf(value) )
            {
                alpha = value;
            }
            spm_cleanenv( str );
        }

        str = spm_getenv( "SPM_FAKE_BETA" );
        if ( str != NULL ) {
            value = strtod( str, NULL );
            if ( (value != HUGE_VAL) && (value != 0.) &&
                 !isnan(value) && !isinf(value) )
            {
                beta = value;
            }
            spm_cleanenv( str );
        }
    }

    degrees = malloc( spm->n * sizeof(spm_int_t));
    diagval = spm_compute_degrees( spm, degrees );
    if ( diagval != spm->n ) {
        /* Diagonal elements must be added to the sparse matrix */
        spm_add_diag( spm, diagval );
        spmUpdateComputedFields( spm );
    }
    spm_generate_fake_values( spm, degrees, alpha, beta );
    free( degrees );

    return;
}
