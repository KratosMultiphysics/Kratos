/**
 * @file z_spm_norm.c
 *
 * SParse Matrix package norm routine.
 *
 * @copyright 2016-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @author Theophile Terraz
 * @date 2015-06-01
 * @precisions normal z -> c d s
 *
 * @addtogroup spm_dev_norm
 * @{
 *
 **/
#include "common.h"
#include "z_spm.h"
#include "frobeniusupdate.h"

/**
 *******************************************************************************
 *
 * @brief Compute the Frobenius norm of the non distributed given
 * spm structure.
 *
 *  ||A|| = sqrt( sum( a_ij ^ 2 ) )
 *
 *******************************************************************************
 *
 * @param[in] spm
 *           The spm from which the norm need to be computed.
 *
 *******************************************************************************
 *
 * @return The computed frobenius norm
 *
 *******************************************************************************/
double
z_spmFrobeniusNorm( const spmatrix_t *spm )
{
    spm_int_t i, j, baseval;
    double *valptr = (double*)spm->values;
    double scale = 1.;
    double sumsq = 0.;

    if (spm->mtxtype == SpmGeneral) {
        for(i=0; i <spm->nnzexp; i++, valptr++) {
            frobenius_update( 1, &scale, &sumsq, valptr );

#if defined(PRECISION_z) || defined(PRECISION_c)
            valptr++;
            frobenius_update( 1, &scale, &sumsq, valptr );
#endif
        }
    }
    else {
        spm_int_t *colptr = spm->colptr;
        spm_int_t *rowptr = spm->rowptr;
        int nb;
        baseval = spmFindBase( spm );

        switch( spm->fmttype ){
        case SpmCSC:
            for(i=0; i<spm->n; i++, colptr++) {
                for(j=colptr[0]; j<colptr[1]; j++, rowptr++, valptr++) {
                    nb = ( i == (*rowptr-baseval) ) ? 1 : 2;
                    frobenius_update( nb, &scale, &sumsq, valptr );

#if defined(PRECISION_z) || defined(PRECISION_c)
                    valptr++;
                    frobenius_update( nb, &scale, &sumsq, valptr );
#endif
                }
            }
            break;
        case SpmCSR:
            for(i=0; i<spm->n; i++, rowptr++) {
                for(j=rowptr[0]; j<rowptr[1]; j++, colptr++, valptr++) {
                    nb = ( i == (*colptr-baseval) ) ? 1 : 2;
                    frobenius_update( nb, &scale, &sumsq, valptr );

#if defined(PRECISION_z) || defined(PRECISION_c)
                    valptr++;
                    frobenius_update( nb, &scale, &sumsq, valptr );
#endif
                }
            }
            break;
        case SpmIJV:
            for(i=0; i <spm->nnz; i++, valptr++, colptr++, rowptr++) {
                nb = ( *rowptr == *colptr ) ? 1 : 2;
                frobenius_update( nb, &scale, &sumsq, valptr );

#if defined(PRECISION_z) || defined(PRECISION_c)
                valptr++;
                frobenius_update( nb, &scale, &sumsq, valptr );
#endif
            }
            break;
        default:
            fprintf(stderr, "z_spmFrobeniusNorm: Unknown Format\n");
        }
    }

    return scale * sqrt( sumsq );
}

/**
 *******************************************************************************
 *
 * @brief Compute the Max norm of the non distributed given spm
 * structure.
 *
 *  ||A|| = max( abs(a_ij) )
 *
 *******************************************************************************
 *
 * @param[in] spm
 *           The spm from which the norm need to be computed.
 *
 *******************************************************************************
 *
 * @return The computed max norm
 *
 *******************************************************************************/
double
z_spmMaxNorm( const spmatrix_t *spm )
{
    spm_int_t i;
    spm_complex64_t *valptr = (spm_complex64_t*)spm->values;
    double tmp, norm = 0.;

    for(i=0; i <spm->nnzexp; i++, valptr++) {
        tmp = cabs( *valptr );
        norm = norm > tmp ? norm : tmp;
    }

    return norm;
}

/**
 *******************************************************************************
 *
 * @brief Compute the Infinity norm of the non distributed given spm
 * structure given by the maximum column sum.
 *
 *  ||A|| = max_i( sum_j(|a_ij|) )
 *
 *******************************************************************************
 *
 * @param[in] spm
 *           The spm from which the norm need to be computed.
 *
 *******************************************************************************
 *
 * @return The computed infinity norm
 *
 *******************************************************************************/
double
z_spmInfNorm( const spmatrix_t *spm )
{
    spm_int_t col, row, i, baseval;
    spm_complex64_t *valptr = (spm_complex64_t*)spm->values;
    double norm = 0.;
    double *sumrow;

    sumrow = malloc( spm->gN * sizeof(double) );
    memset( sumrow, 0, spm->gN * sizeof(double) );
    baseval = spmFindBase( spm );

    switch( spm->fmttype ){
    case SpmCSC:
        for( col=0; col < spm->gN; col++ )
        {
            for( i=spm->colptr[col]-baseval; i<spm->colptr[col+1]-baseval; i++ )
            {
                row = spm->rowptr[i] - baseval;
                sumrow[row] += cabs( valptr[i] );

                /* Add the symmetric/hermitian part */
                if ( ( spm->mtxtype != SpmGeneral ) &&
                     ( row != col ) )
                {
                    sumrow[col] += cabs( valptr[i] );
                }
            }
        }
        break;

    case SpmCSR:
        for( row=0; row < spm->gN; row++ )
        {
            for( i=spm->rowptr[row]-baseval; i<spm->rowptr[row+1]-baseval; i++ )
            {
                sumrow[row] += cabs( valptr[i] );
            }
        }

        /* Add the symmetric/hermitian part */
        if ( spm->mtxtype != SpmGeneral )
        {
            for( row=0; row < spm->gN; row++ )
            {
                for( i=spm->rowptr[row]-baseval; i<spm->rowptr[row+1]-baseval; i++ )
                {
                    col = spm->colptr[i] - baseval;
                    if ( row != col ) {
                        sumrow[col] += cabs( valptr[i] );
                    }
                }
            }
        }
        break;

    case SpmIJV:
        for(i=0; i < spm->nnz; i++)
        {
            row = spm->rowptr[i]-baseval;
            sumrow[row] += cabs( valptr[i] );
        }

        /* Add the symmetric/hermitian part */
        if ( spm->mtxtype != SpmGeneral )
        {
            for(i=0; i < spm->nnz; i++)
            {
                row = spm->rowptr[i]-baseval;
                col = spm->colptr[i]-baseval;
                if( row != col ) {
                    sumrow[col] += cabs( valptr[i] );
                }
            }
        }
        break;

    default:
        free( sumrow );
        return SPM_ERR_BADPARAMETER;
    }

    for( i=0; i<spm->gN; i++)
    {
        if(norm < sumrow[i])
        {
            norm = sumrow[i];
        }
    }
    free( sumrow );

    return norm;
}

/**
 *******************************************************************************
 *
 * @brief  Compute the one norm of the non distributed given spm
 * structure fiven by the maximum row sum
 *
 *  ||A|| = max_j( sum_i(|a_ij|) )
 *
 *******************************************************************************
 *
 * @param[in] spm
 *           The spm from which the norm need to be computed.
 *
 *******************************************************************************
 *
 * @return The computed one norm
 *
 *******************************************************************************/
double
z_spmOneNorm( const spmatrix_t *spm )
{
    spm_int_t col, row, i, baseval;
    spm_complex64_t *valptr = (spm_complex64_t*)spm->values;
    double norm = 0.;
    double *sumcol;

    sumcol = malloc( spm->gN * sizeof(double) );
    memset( sumcol, 0, spm->gN * sizeof(double) );
    baseval = spmFindBase( spm );

    switch( spm->fmttype ){
    case SpmCSC:
        for( col=0; col<spm->gN; col++ )
        {
            for( i=spm->colptr[col]-baseval; i<spm->colptr[col+1]-baseval; i++ )
            {
                sumcol[col] += cabs( valptr[i] );
            }
        }

        /* Add the symmetric/hermitian part */
        if ( spm->mtxtype != SpmGeneral )
        {
            for( col=0; col < spm->gN; col++ )
            {
                for( i=spm->colptr[col]-baseval; i<spm->colptr[col+1]-baseval; i++ )
                {
                    row = spm->rowptr[i] - baseval;
                    if (row != col) {
                        sumcol[row] += cabs( valptr[i] );
                    }
                }
            }
        }
        break;

    case SpmCSR:
        for( row=0; row < spm->gN; row++ )
        {
            for( i=spm->rowptr[row]-baseval; i<spm->rowptr[row+1]-baseval; i++ )
            {
                col = spm->colptr[i] - baseval;
                sumcol[col] += cabs( valptr[i] );

                /* Add the symmetric/hermitian part */
                if ( ( spm->mtxtype != SpmGeneral ) &&
                     ( row != col ) )
                {
                    sumcol[row] += cabs( valptr[i] );
                }
            }
        }
        break;

    case SpmIJV:
        for(i=0; i < spm->nnz; i++)
        {
            sumcol[spm->colptr[i]-baseval] += cabs( valptr[i] );
        }

        /* Add the symmetric/hermitian part */
        if ( spm->mtxtype != SpmGeneral )
        {
            for(i=0; i < spm->nnz; i++)
            {
                if(spm->rowptr[i] != spm->colptr[i])
                    sumcol[spm->rowptr[i]-baseval] += cabs( valptr[i] );
            }
        }
        break;

    default:
        free( sumcol );
        return SPM_ERR_BADPARAMETER;
    }

    for( i=0; i<spm->gN; i++)
    {
        if(norm < sumcol[i])
        {
            norm = sumcol[i];
        }
    }
    free( sumcol );

    return norm;
}

/**
 *******************************************************************************
 *
 * @brief Compute the norm of an spm matrix
 *
 *******************************************************************************
 *
 * @param[in] ntype
 *          = SpmMaxNorm: Max norm
 *          = SpmOneNorm: One norm
 *          = SpmInfNorm: Infinity norm
 *          = SpmFrobeniusNorm: Frobenius norm
 *
 * @param[in] spm
 *          The spm structure describing the matrix.
 *
 *******************************************************************************
 *
 * @return The norm of the spm matrix
 *
 *******************************************************************************/
double
z_spmNorm( int ntype,
           const spmatrix_t *spm )
{
    double norm = 0.;

    if(spm == NULL)
    {
        return -1.;
    }

    switch( ntype ) {
    case SpmMaxNorm:
        norm = z_spmMaxNorm( spm );
        break;

    case SpmInfNorm:
        norm = z_spmInfNorm( spm );
        break;

    case SpmOneNorm:
        norm = z_spmOneNorm( spm );
        break;

    case SpmFrobeniusNorm:
        norm = z_spmFrobeniusNorm( spm );
        break;

    default:
        fprintf(stderr, "z_spmNorm: invalid norm type\n");
        return -1.;
    }

    return norm;
}
/**
 * @}
 */
