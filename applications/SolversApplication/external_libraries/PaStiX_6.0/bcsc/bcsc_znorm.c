/**
 *
 * @file bcsc_znorm.c
 *
 * @copyright 2004-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.1
 * @author Mathieu Faverge
 * @author Théophile terraz
 * @date 2018-07-16
 * @precisions normal z -> c d s
 *
 **/
#include "common.h"
#include "bcsc.h"
#include <math.h>
#include "frobeniusupdate.h"

/**
 *******************************************************************************
 *
 * @ingroup bcsc_internal
 *
 * @brief Compute the max norm of a bcsc matrix.
 *
 *******************************************************************************
 *
 * @param[in] values
 *          The values array of the matrix.
 *
 *******************************************************************************
 *
 * @retval The norm of the matrix.
 *
 *******************************************************************************/
double
bcsc_znorm_max( const pastix_bcsc_t *bcsc )
{
    double temp;
    double norm = 0.;
    pastix_complex64_t *valptr = (pastix_complex64_t*)bcsc->Lvalues;
    pastix_int_t i, j, bloc;

    for( bloc=0; bloc < bcsc->cscfnbr; bloc++ )
    {
        for( j=0; j < bcsc->cscftab[bloc].colnbr; j++ )
        {
            for( i = bcsc->cscftab[bloc].coltab[j]; i < bcsc->cscftab[bloc].coltab[j+1]; i++ )
            {
                temp = cabs(valptr[i]);
                if(norm < temp)
                {
                    norm = temp;
                }
            }
        }
    }

    return norm;
}

/**
 *******************************************************************************
 *
 * @ingroup bcsc_internal
 *
 * @brief Compute the infinity norm of a bcsc matrix.
 * The infinity norm is equal to the maximum value of the sum of the
 * absolute values of the elements of each rows.
 *
 *******************************************************************************
 *
 * @param[in] bcsc
 *          The Pastix bcsc.
 *
 *******************************************************************************
 *
 * @retval The norm of the matrix.
 *
 *******************************************************************************/
double
bcsc_znorm_inf( const pastix_bcsc_t *bcsc )
{
    double norm = 0.;
    pastix_complex64_t *valptr;
    int i, j, bloc;

    if( bcsc->Uvalues != NULL )
    {
        double sum;

        valptr = (pastix_complex64_t*)(bcsc->Uvalues);
        for( bloc=0; bloc<bcsc->cscfnbr; bloc++ )
        {
            for( j=0; j<bcsc->cscftab[bloc].colnbr; j++ )
            {
                sum = 0.;
                for( i = bcsc->cscftab[bloc].coltab[j]; i < bcsc->cscftab[bloc].coltab[j+1]; i++ )
                {
                    sum += cabs(valptr[i]);
                }
                if( sum > norm ) {
                    norm = sum;
                }
            }
        }
    }
    else {
        double *sumrow;
        valptr = (pastix_complex64_t*)bcsc->Lvalues;

        MALLOC_INTERN( sumrow, bcsc->n, double);
        memset( sumrow, 0, bcsc->n * sizeof(double) );

        for( bloc=0; bloc < bcsc->cscfnbr; bloc++ )
        {
            for( j=0; j < bcsc->cscftab[bloc].colnbr; j++ )
            {
                for( i = bcsc->cscftab[bloc].coltab[j]; i < bcsc->cscftab[bloc].coltab[j+1]; i++ )
                {
                    sumrow[ bcsc->rowtab[i] ] += cabs(valptr[i]);
                }
            }
        }

        for( i=0; i<bcsc->n; i++)
        {
            if(norm < sumrow[i])
            {
                norm = sumrow[i];
            }
        }
        memFree_null( sumrow );
    }

    return norm;
}

/**
 *******************************************************************************
 *
 * @ingroup bcsc_internal
 *
 * @brief Compute the norm 1 of a bcsc matrix.
 * Norm 1 is equal to the maximum value of the sum of the
 * absolute values of the elements of each columns.
 *
 *******************************************************************************
 *
 * @param[in] bcsc
 *          The Pastix bcsc.
 *
 *******************************************************************************
 *
 * @retval The norm of the matrix.
 *
 *******************************************************************************/
double
bcsc_znorm_one( const pastix_bcsc_t *bcsc )
{
    pastix_complex64_t *valptr = (pastix_complex64_t*)bcsc->Lvalues;
    double sum, norm = 0.;
    int i, j, bloc;

    for( bloc=0; bloc<bcsc->cscfnbr; bloc++ )
    {
        for( j=0; j<bcsc->cscftab[bloc].colnbr; j++ )
        {
            sum = 0.;
            for( i = bcsc->cscftab[bloc].coltab[j]; i < bcsc->cscftab[bloc].coltab[j+1]; i++ )
            {
                sum += cabs(valptr[i]);
            }
            if( sum > norm ) {
                norm = sum;
            }
        }
    }

    return norm;
}

/**
 *******************************************************************************
 *
 * @ingroup bcsc_internal
 *
 * @brief Compute the frobenius norm of a bcsc matrix.
 *
 *******************************************************************************
 *
 * @param[in] bcsc
 *          The Pastix bcsc.
 *
 *******************************************************************************
 *
 *          The norm of the matrix
 *
 *******************************************************************************/
double
bcsc_znorm_frobenius( const pastix_bcsc_t *bcsc )
{
    double scale = 0.;
    double sum = 1.;
    double norm;
    double *valptr = (double*)bcsc->Lvalues;
    pastix_int_t i, j, bloc;

    for( bloc=0; bloc < bcsc->cscfnbr; bloc++ )
    {
        for( j=0; j < bcsc->cscftab[bloc].colnbr; j++ )
        {
            for( i = bcsc->cscftab[bloc].coltab[j]; i < bcsc->cscftab[bloc].coltab[j+1]; i++, valptr++ )
            {
                frobenius_update( 1, &scale, &sum, valptr);
#if defined(PRECISION_z) || defined(PRECISION_c)
                valptr++;
                frobenius_update( 1, &scale, &sum, valptr);
#endif
            }
        }
    }

    norm = scale*sqrt(sum);

    return norm;
}

/**
 *******************************************************************************
 *
 * @brief Compute the norm of an bcsc matrix
 *
 *******************************************************************************
 *
 * @param[in] ntype
 *          = PastixMaxNorm: Max norm
 *          = PastixOneNorm: One norm
 *          = PastixInfNorm: Infinity norm
 *          = PastixFrobeniusNorm: Frobenius norm
 *
 * @param[in] bcsc
 *          The bcsc structure describing the matrix.
 *
 *******************************************************************************
 *
 * @retval The norm of the matrix.
 *
 *******************************************************************************/
double
bcsc_znorm( pastix_normtype_t    ntype,
            const pastix_bcsc_t *bcsc )
{
    double norm = 0.;

    if(bcsc == NULL)
    {
        return -1.;
    }

    switch( ntype ) {
    case PastixMaxNorm:
        norm = bcsc_znorm_max( bcsc );
        break;

    case PastixInfNorm:
        norm = bcsc_znorm_inf( bcsc );
        break;

    case PastixOneNorm:
        norm = bcsc_znorm_one( bcsc );
        break;

    case PastixFrobeniusNorm:
        norm = bcsc_znorm_frobenius( bcsc );
        break;

    default:
        fprintf(stderr, "bcsc_znorm: invalid norm type\n");
        return -1.;
    }

    return norm;
}
