/**
 *
 * @file readhb.c
 *
 * @copyright 2014-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.0
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @date 2011-11-11
 *
 **/
#include "common.h"
#include "spm_drivers.h"
#include "drivers/iohb.h"

/**
 *******************************************************************************
 *
 * @ingroup spm_driver
 *
 * @brief Interface to the Harwell-Boeing C driver (iohb.c)
 *
 *******************************************************************************
 *
 * @param[in] filename
 *          The file containing the matrix.
 *
 * @param[in] spm
 *          At exit, contains the matrix in spm format.
 *
 *******************************************************************************
 *
 * @retval SPM_SUCCESS if the matrix has been read successfully
 * @retval SPM_ERR_IO if a problem occured in the Harwell Boeing driver
 * @retval SPM_ERR_BADPARAMETER if the matrix is no in a supported format
 *
 *******************************************************************************/
int
readHB( const char   *filename,
        spmatrix_t *spm )
{
    int M, N, nz, nrhs;

    /* Harwell Boeing is a variant of RSA */
    spm->fmttype = SpmCSC;
    spm->dof     = 1;
    spm->loc2glob= NULL;

    /* Read header informations */
    {
        char *Type;
        Type = malloc(4*sizeof(char));
        Type[0] ='a';

        readHB_info(filename, &M, &N, &nz, &Type, &nrhs);

        if ( M != N ) {
            fprintf(stderr, "readHB: PaStiX does not support non square matrices (m=%d, N=%d\n", M, N);
            return SPM_ERR_BADPARAMETER;
        }

        spm->gN   = M;
        spm->n    = M;
        spm->gnnz = nz;
        spm->nnz  = nz;

        /* Check float type */
        switch( Type[0] ) {
        case 'C':
        case 'c':
            spm->flttype = SpmComplex64;
            break;
        case 'R':
        case 'r':
            spm->flttype = SpmDouble;
            break;
        case 'P':
        case 'p':
            spm->flttype = SpmPattern;
            break;
        default:
            fprintf(stderr, "readhb: Floating type unknown (%c)\n", Type[0]);
            return SPM_ERR_BADPARAMETER;
        }

        /* Check Symmetry */
        switch( Type[1] ) {
        case 'S':
        case 's':
            spm->mtxtype = SpmSymmetric;
            break;
        case 'H':
        case 'h':
            spm->mtxtype = SpmHermitian;
            assert( spm->flttype == SpmDouble );
            break;
        case 'U':
        case 'u':
        default:
            spm->mtxtype = SpmGeneral;
        }
        free(Type);
    }

    /* Read the matrix and its values */
    {
        int    *colptr, *rowind;
        int     rc;

        rc = readHB_newmat_double( filename, &M, &N, &nz,
                                   &colptr, &rowind, (double**)(&(spm->values)) );

        if (rc == 0) {
            fprintf(stderr, "readhb: Error in reading the HB matrix values\n");
            return SPM_ERR_IO;
        }

        /* Move the colptr/rowind from int to spm_int_t if different sizes */
        spm->colptr = spmIntConvert(spm->n+1, colptr);
        spm->rowptr = spmIntConvert(spm->nnz, rowind);
    }
    return SPM_SUCCESS;
}
