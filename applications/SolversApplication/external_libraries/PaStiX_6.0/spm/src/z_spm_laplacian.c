/**
 * @file z_spm_laplacian.c
 *
 * SParse Matrix package laplacian generator routines.
 *
 * @copyright 2016-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Xavier Lacoste
 * @author Theophile Terraz
 * @date 2015-01-01
 * @precisions normal z -> c d s p
 *
 **/
#include "common.h"
#include "drivers/laplacian.h"

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_driver
 *
 * @brief Generate a laplacian matrix for a 3D 7-points stencil
 * \f[ M = \alpha * D - \beta * A \f]
 *
 * Example:
 * >  3 -1 -1  0 -1  0  0  0
 * > -1  3  0 -1  0 -1  0  0
 * > -1  0  3 -1  0  0 -1  0
 * >  0 -1 -1  3  0  0  0 -1
 * > -1  0  0  0  3 -1 -1  0
 * >  0 -1  0  0 -1  3  0 -1
 * >  0  0 -1  0 -1  0  3 -1
 * >  0  0  0 -1  0 -1 -1  3
 *
 *******************************************************************************
 *
 * @param[inout] spm
 *          At start, an allocated spm structure.
 *          Contains the size of the laplacian in spm->n.
 *          At exit, contains the matrix in csc format.
 *
 * @param[in] dim1
 *          contains the first dimension of the grid of the laplacian.
 *
 * @param[in] dim2
 *          contains the second dimension of the grid of the laplacian.
 *
 * @param[in] dim3
 *          contains the third dimension of the grid of the laplacian.
 *
 * @param[in] alpha
 *          The alpha coefficient for the degree matrix
 *
 * @param[in] beta
 *          The beta coefficient for the adjacency matrix
 *
 * @remark: In complex, the Laplacian is set to hermitian. See
 * z_spmLaplacian_27points() to get a symmetric Laplacian, or change the
 * mtxtype field by hand.
 *
 *******************************************************************************/
void
z_spmLaplacian_7points( spmatrix_t   *spm,
                        spm_int_t    dim1,
                        spm_int_t    dim2,
                        spm_int_t    dim3,
                        spm_fixdbl_t alpha,
                        spm_fixdbl_t beta )
{

    spm_complex64_t *valptr;
    spm_complex64_t  lalpha = (spm_complex64_t)alpha;
    spm_complex64_t  lbeta  = (spm_complex64_t)beta;
    spm_int_t *colptr, *rowptr;
    spm_int_t i, j, k, l, degree = 0;
    spm_int_t nnz = (2*(dim1)-1)*dim2*dim3 + (dim2-1)*dim1*dim3 + dim2*dim1*(dim3-1);

    spm->mtxtype  = SpmHermitian;
    spm->flttype  = SpmComplex64;
    spm->fmttype  = SpmCSC;
    spm->nnz      = nnz;
    spm->dof      = 1;

    assert( spm->n == dim1*dim2*dim3 );

    /* Allocating */
    spm->colptr = malloc((spm->n+1)*sizeof(spm_int_t));
    spm->rowptr = malloc(nnz       *sizeof(spm_int_t));
    assert( spm->colptr );
    assert( spm->rowptr );

#if !defined(PRECISION_p)
    spm->values = malloc(nnz       *sizeof(spm_complex64_t));
    assert( spm->values );
#endif

    /* Building ia, ja and values*/
    colptr = spm->colptr;
    rowptr = spm->rowptr;
    valptr = (spm_complex64_t*)(spm->values);

    /* Building ia, ja and values*/
    *colptr = 1;
    l = 1; /* Column index in the matrix ((i-1) * dim1 * dim2 + (j-1) * dim1 + k-1) */

    /* Start with one for each dimension (top corner) */
    degree = 3;
    for(i=0; i<dim3; i++)
    {
        /* +1 at the second range */
        if ( i == 1 ) {
            degree++;
        }
        if ( i == (dim3-1) ) {
            degree--;
        }

        for(j=0; j<dim2; j++)
        {
            /* +1 at the second range */
            if ( j == 1 ) {
                degree++;
            }
            if ( j == (dim2-1) ) {
                degree--;
            }

            for(k=0; k<dim1; k++)
            {
                colptr[1] = colptr[0];

                /* +1 at the second range */
                if ( k == 1 ) {
                    degree++;
                }
                if ( k == (dim1-1) ) {
                    degree--;
                }

                /* Diagonal value */
                *rowptr = l;
#if !defined(PRECISION_p)
                *valptr = (spm_complex64_t)degree * lalpha;
#endif
                valptr++;
                rowptr++;
                colptr[1]++;

                /* Connexion along dimension 1 */
                if (k < (dim1-1)) {
                    *rowptr = l+1;
#if !defined(PRECISION_p)
                    *valptr = -lbeta;
#endif
                    valptr++;
                    rowptr++;
                    colptr[1]++;
                }

                /* Connexion along dimension 2 */
                if (j < (dim2-1)) {
                    *rowptr = l+dim1;
#if !defined(PRECISION_p)
                    *valptr = -lbeta;
#endif
                    valptr++;
                    rowptr++;
                    colptr[1]++;
                }

                /* Connexion along dimension 3 */
                if (i < (dim3-1)) {
                    *rowptr = l+dim1*dim2;
#if !defined(PRECISION_p)
                    *valptr = -lbeta;
#endif
                    valptr++;
                    rowptr++;
                    colptr[1]++;
                }

                colptr++;
                l++;
            }
        }
    }

    assert( (spm->colptr[ spm->n ] - spm->colptr[0]) == nnz );
    (void)lalpha;
    (void)lbeta;
    (void)degree;
}

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_driver
 *
 * @brief Generate an extended laplacian matrix for a 3D 27-points stencil with
 * \f[ M = \alpha * D - \beta * A \f], where D is the matrix of degrees, and A
 * the matrix of adjacency with coefficients of 1 for B connexions, 1/ sqrt(2)
 * for X connexions, and 1/sqrt(3) for D connexions.
 *
 *      D-------X-------D
 *     /|      /|      /|
 *    X-------B-------X |
 *   /| |    /| |    /| |
 *  D-------X-|-----X | |
 *  | | X---|-|-B---|-|-X
 *  | |/|   | |/    | |/|
 *  | B-----|-A-----|-B |
 *  |/| |   |/|     |/| |
 *  X-------B-------X | |
 *  | | D---|-|-X---|-|-D
 *  | |/    | |/    | |/
 *  | X-----|-B-----|-X
 *  |/      |/      |/
 *  D-------X-------D
 *
 * @remark: In complex, the Laplacian is set to symmetric. See
 * z_spmLaplacian_7points() to get an hermitian Laplacian, or change the
 * mtxtype field by hand.
 *
 *******************************************************************************
 *
 * @param[inout] spm
 *          At start, an allocated spm structure.
 *          Contains the size of the laplacian in spm->n.
 *          At exit, contains the matrix in csc format.
 *
 * @param[in] dim1
 *          contains the first dimension of the grid of the laplacian.
 *
 * @param[in] dim2
 *          contains the second dimension of the grid of the laplacian.
 *
 * @param[in] dim3
 *          contains the third dimension of the grid of the laplacian.
 *
 * @param[in] alpha
 *          The alpha coefficient for the degree matrix
 *
 * @param[in] beta
 *          The beta coefficient for the adjacency matrix
 *
 *******************************************************************************/
void
z_spmLaplacian_27points( spmatrix_t   *spm,
                         spm_int_t    dim1,
                         spm_int_t    dim2,
                         spm_int_t    dim3,
                         spm_fixdbl_t alpha,
                         spm_fixdbl_t beta )
{

    spm_complex64_t *valptr;
    /*
     * See https://crd.lbl.gov/assets/pubs_presos/iwapt09-27pt.pdf for the
     * meaning of alpha, beta, gamma, and delta.
     * "Auto-tuning the 27-point Stencil for Multicore", K. Datta, S. Williams,
     * V. Volkov, J. Carter, L. Oliker, J. Shalf, and K. Yelick
     */
    spm_complex64_t  lalpha = (spm_complex64_t)alpha;
    spm_complex64_t  lbeta  = (spm_complex64_t)beta;
    spm_complex64_t  lgamma = (spm_complex64_t)beta / sqrt(2.);
    spm_complex64_t  ldelta = (spm_complex64_t)beta / sqrt(3.);
    spm_int_t *colptr, *rowptr;
    spm_int_t i, j, k, l, degree, d;
    spm_int_t nnz = (2*dim1-1) *  dim2    * dim3
        +           (3*dim1-2) * (dim2-1) * dim3
        +          ((3*dim1-2) *  dim2 + 2 * (3*dim1-2) *(dim2-1)) * (dim3-1);

    spm->mtxtype  = SpmSymmetric;
    spm->flttype  = SpmComplex64;
    spm->fmttype  = SpmCSC;
    spm->nnz      = nnz;
    spm->dof      = 1;

    assert( spm->n == dim1*dim2*dim3 );

    /* Allocating */
    spm->colptr = malloc((spm->n+1)*sizeof(spm_int_t));
    spm->rowptr = malloc(nnz       *sizeof(spm_int_t));
    assert( spm->colptr );
    assert( spm->rowptr );

#if !defined(PRECISION_p)
    spm->values = malloc(nnz       *sizeof(spm_complex64_t));
    assert( spm->values );
#endif

    /* Building ia, ja and values*/
    colptr = spm->colptr;
    rowptr = spm->rowptr;
    valptr = (spm_complex64_t*)(spm->values);

    /* Building ia, ja and values*/
    *colptr = 1;
    l = 1; /* Column index in the matrix ((i-1) * dim1 * dim2 + (j-1) * dim1 + k-1) */
    for(i=1; i<=dim3; i++)
    {
        for(j=1; j<=dim2; j++)
        {
            for(k=1; k<=dim1; k++)
            {
                colptr[1] = colptr[0];

                /* Diagonal value */
                *rowptr = l;
                degree = 1;
                d = 1;

#if !defined(PRECISION_p)
                if (k > 1) {
                    d++;
                }
                if (k < dim1) {
                    d++;
                }
                degree = degree * d;
                d = 1;

                if (j > 1) {
                    d++;
                }
                if (j < dim2) {
                    d++;
                }
                degree = degree * d;
                d = 1;

                if (i > 1) {
                    d++;
                }
                if (i < dim3) {
                    d++;
                }
                degree = degree * d - 1;

                *valptr = (spm_complex64_t)degree * lalpha;
#endif
                valptr++;
                rowptr++;
                colptr[1]++;

                /* Connexion along dimension 1 */
                if (k < dim1) {
                    *rowptr = l+1;
#if !defined(PRECISION_p)
                    *valptr = -lbeta;
#endif
                    valptr++;
                    rowptr++;
                    colptr[1]++;
                }

                /* Connexion along dimension 2 */
                if (j < dim2)
                {
                    if (k > 1)
                    {
                        *rowptr = l+dim1-1;
#if !defined(PRECISION_p)
                        *valptr = -lgamma;
#endif
                        valptr++;
                        rowptr++;
                        colptr[1]++;
                    }

                    *rowptr = l+dim1;
#if !defined(PRECISION_p)
                    *valptr = -lbeta;
#endif
                    valptr++;
                    rowptr++;
                    colptr[1]++;

                    if (k < dim1)
                    {
                        *rowptr = l+dim1+1;
#if !defined(PRECISION_p)
                        *valptr = -lgamma;
#endif
                        valptr++;
                        rowptr++;
                        colptr[1]++;
                    }
                }

                /* Connexion along dimension 3 */
                if (i < dim3) {
                    if( j > 1 )
                    {
                        if (k > 1)
                        {
                            *rowptr = l+dim1*dim2-dim1-1;
#if !defined(PRECISION_p)
                            *valptr = -ldelta;
#endif
                            valptr++;
                            rowptr++;
                            colptr[1]++;
                        }

                        *rowptr = l+dim1*dim2-dim1;
#if !defined(PRECISION_p)
                        *valptr = -lgamma;
#endif
                        valptr++; rowptr++; colptr[1]++;

                        if (k < dim1)
                        {
                            *rowptr = l+dim1*dim2-dim1+1;
#if !defined(PRECISION_p)
                            *valptr = -ldelta;
#endif
                            valptr++;
                            rowptr++;
                            colptr[1]++;
                        }
                    }
                    if (k > 1)
                    {
                        *rowptr = l+dim1*dim2-1;
#if !defined(PRECISION_p)
                        *valptr = -lgamma;
#endif
                        valptr++;
                        rowptr++;
                        colptr[1]++;
                    }

                    *rowptr = l+dim1*dim2;
#if !defined(PRECISION_p)
                    *valptr = -lbeta;
#endif
                    valptr++;
                    rowptr++;
                    colptr[1]++;

                    if (k < dim1)
                    {
                        *rowptr = l+dim1*dim2+1;
#if !defined(PRECISION_p)
                        *valptr = -lgamma;
#endif
                        valptr++;
                        rowptr++;
                        colptr[1]++;
                    }

                    if( j < dim2 )
                    {
                        if (k > 1)
                        {
                            *rowptr = l+dim1*dim2+dim1-1;
#if !defined(PRECISION_p)
                            *valptr = -ldelta;
#endif
                            valptr++;
                            rowptr++;
                            colptr[1]++;
                        }

                        *rowptr = l+dim1*dim2+dim1;
#if !defined(PRECISION_p)
                        *valptr = -lgamma;
#endif
                        valptr++;
                        rowptr++;
                        colptr[1]++;

                        if (k < dim1)
                        {
                            *rowptr = l+dim1*dim2+dim1+1;
#if !defined(PRECISION_p)
                            *valptr = -ldelta;
#endif
                            valptr++;
                            rowptr++;
                            colptr[1]++;
                        }
                    }
                }

                colptr++; l++;
            }
        }
    }

    assert( (spm->colptr[ spm->n ] - spm->colptr[0]) == nnz );
    (void)lalpha;
    (void)lbeta;
    (void)lgamma;
    (void)ldelta;
    (void)degree;
    (void)d;
}
