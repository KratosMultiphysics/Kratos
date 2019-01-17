/**
 *
 * @file readmm.c
 *
 * @copyright 2011-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.0
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @author Theophile Terraz
 * @date 2011-11-11
 *
 **/
#include "common.h"
#include "spm_drivers.h"
#include "drivers/mmio.h"

/**
 *******************************************************************************
 *
 * @ingroup spm_driver
 *
 * @brief Read the data part of a complex matrix in Matrix Market file.
 *
 * For more information about matrix market format see mmio.c/mmio.h
 *
 *******************************************************************************
 *
 * @param[in] file
 *          The file opened in readMM which contains the matrix stored in Matrix
 *          Market format.
 *
 * @param[inout] spm
 *          At exit, the data of the matrix are stored in the spm structure.
 *
 *******************************************************************************
 *
 * @retval SPM_SUCCESS if the matrix has been read successfully
 * @retval SPM_ERR_IO if a problem occured in the RSA driver
 *
 *******************************************************************************/
int
z_readMM( FILE *file,
          spmatrix_t *spm )
{
    spm_complex64_t *valptr;
    spm_int_t *colptr;
    spm_int_t *rowptr;
    spm_int_t i;
    long row, col;
    double re, im;

    spm->values = malloc( spm->nnz * sizeof(spm_complex64_t) );

    colptr = spm->colptr;
    rowptr = spm->rowptr;
    valptr = (spm_complex64_t*)(spm->values);

    for (i=0; i<spm->nnz; i++, colptr++, rowptr++, valptr++)
    {
        if (4 != fscanf(file,"%ld %ld %lg %lg\n", &row, &col, &re, &im))
        {
            fprintf(stderr, "readmm: erro while reading matrix file (line %ld)\n", (long)i);
            return SPM_ERR_IO;
        }

        *rowptr = (spm_int_t)row;
        *colptr = (spm_int_t)col;
        *valptr = (spm_complex64_t)(re + im * I);
    }

    return SPM_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup spm_driver
 *
 * @brief Read the data part of a real matrix in Matrix Market file.
 * For more information about matrix market format see mmio.c/mmio.h
 *
 *******************************************************************************
 *
 * @param[in] file
 *          The file opened in readMM which contains the matrix stored in Matrix
 *          Market format.
 *
 * @param[inout] spm
 *          At exit, the data of the matrix are stored in the spm structure.
 *
 *******************************************************************************
 *
 * @retval SPM_SUCCESS if the matrix has been read successfully
 * @retval SPM_ERR_IO if a problem occured in the MM driver
 *
 *******************************************************************************/
int
d_readMM( FILE *file,
          spmatrix_t *spm )
{
    double       *valptr;
    spm_int_t *colptr;
    spm_int_t *rowptr;
    spm_int_t i;
    long row, col;
    double re;

    spm->values = malloc( spm->nnz * sizeof(double) );

    colptr = spm->colptr;
    rowptr = spm->rowptr;
    valptr = (double*)(spm->values);

    for (i=0; i<spm->nnz; i++, colptr++, rowptr++, valptr++)
    {
        if (3 != fscanf(file,"%ld %ld %lg\n", &row, &col, &re))
        {
            fprintf(stderr, "readmm: erro while reading matrix file (line %ld)\n", (long)i);
            return SPM_ERR_IO;
        }

        *rowptr = (spm_int_t)row;
        *colptr = (spm_int_t)col;
        *valptr = re;
    }

    return SPM_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup spm_driver
 *
 * @brief Read the data part of a pattern matrix in Matrix Market file.
 * For more information about matrix market format see mmio.c/mmio.h
 *
 *******************************************************************************
 *
 * @param[in] file
 *          The file opened in readMM which contains the matrix stored in Matrix
 *          Market format.
 *
 * @param[inout] spm
 *          At exit, the data of the matrix are stored in the spm structure.
 *
 *******************************************************************************
 *
 * @retval SPM_SUCCESS if the matrix has been read successfully
 * @retval SPM_ERR_IO if a problem occured in the MM driver
 *
 *******************************************************************************/
int
p_readMM( FILE *file,
          spmatrix_t *spm )
{
    spm_int_t *colptr;
    spm_int_t *rowptr;
    spm_int_t i;
    long row, col;

    spm->values = NULL;

    colptr = spm->colptr;
    rowptr = spm->rowptr;

    for (i=0; i<spm->nnz; i++, colptr++, rowptr++)
    {
        if (2 != fscanf(file,"%ld %ld\n", &row, &col))
        {
            fprintf(stderr, "readmm: erro while reading matrix file (line %ld)\n", (long)i);
            return SPM_ERR_IO;
        }

        *rowptr = (spm_int_t)row;
        *colptr = (spm_int_t)col;
    }

    return SPM_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup spm_driver
 *
 * @brief Read a matrix in Matrix Market fill. This corresponds to
 * IJV format with (%d %d[ %lf[ %lf]]) format per line.
 * For more information about matrix market format see mmio.c/mmio.h
 *
 *******************************************************************************
 *
 * @param[in] filename
 *          The filename that contains the matrix stored in Matrix Market format.
 *
 * @param[in] spm
 *          At exit, contains the matrix in spm format.
 *
 *******************************************************************************
 *
 * @retval SPM_SUCCESS if the matrix has been read successfully
 * @retval SPM_ERR_IO if a problem occured in the MM driver
 * @retval SPM_ERR_BADPARAMETER if the matrix is no in a supported format
 *
 *******************************************************************************/
int
readMM( const char   *filename,
        spmatrix_t *spm )
{
    MM_typecode matcode;
    FILE *file;
    int rc;

    file = fopen(filename,"r");
    if (file == NULL)
    {
        fprintf(stderr,"readmm: Cannot open the file (%s)\n", filename);
        return SPM_ERR_BADPARAMETER;
    }

    if (mm_read_banner(file, &matcode) != 0)
    {
        fprintf(stderr,"readmm: Could not process Matrix Market banner.\n");
        return SPM_ERR_IO;
    }

    /* Float values type */

    if (mm_is_complex(matcode)) {
        spm->flttype = SpmComplex64;
    }
    else if (mm_is_real(matcode)) {
        spm->flttype = SpmDouble;
    }
    else if (mm_is_pattern(matcode)) {
        spm->flttype = SpmPattern;
    }
    else {
        fprintf(stderr,"readmm: Unsupported type of matrix.\n");
        return SPM_ERR_BADPARAMETER;
    }

    /* Matrix structure */
    if (mm_is_general(matcode)) {
        spm->mtxtype = SpmGeneral;
    }
    else if (mm_is_symmetric(matcode)) {
        spm->mtxtype = SpmSymmetric;
    }
    else if (mm_is_hermitian(matcode)) {
        spm->mtxtype = SpmHermitian;
    }
    else {
        fprintf(stderr,"readmm: Unsupported type of matrix.\n");
        return SPM_ERR_BADPARAMETER;
    }

    spm->fmttype = SpmIJV;
    spm->dof     = 1;
    spm->loc2glob= NULL;

    /* Read the size */
    {
        int m, n, nnz;
        if (mm_read_mtx_crd_size(file, &m, &n, &nnz) != 0) {
            fprintf(stderr, "readmm: error while reading matrix sizes\n");
            return SPM_ERR_IO;
        }

        spm->gN   = n;
        spm->n    = n;
        spm->gnnz = nnz;
        spm->nnz  = nnz;
    }

    spm->colptr = (spm_int_t*)malloc(spm->nnz * sizeof(spm_int_t));
    spm->rowptr = (spm_int_t*)malloc(spm->nnz * sizeof(spm_int_t));

    switch( spm->flttype ) {
    case SpmComplex64:
        rc = z_readMM(file, spm);
        break;

    case SpmDouble:
        rc = d_readMM(file, spm);
        break;

    case SpmPattern:
    default:
        rc = p_readMM(file, spm);
    }

    fclose(file);
    return rc;
}
