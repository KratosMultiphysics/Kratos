/**
 *
 * @file spm_io.c
 *
 * SParse Matrix package I/O routines.
 *
 * @copyright 2016-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.0.0
 * @author Xavier Lacoste
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @date 2013-06-24
 *
 **/
#include "common.h"

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_io
 *
 * @brief Read an array of integer.
 *
 *******************************************************************************
 *
 * @param[in] stream
 *          The opened file in which the spm is stored.
 *
 * @param[in] n
 *          Number of elements.
 *
 * @param[out] array
 *          Array of results.
 *
 ********************************************************************************
 *
 * @retval SPM_SUCCESS if the read happened successfully,
 * @retval SPM_ERR_FILE if the input format is incorrect.
 *
 *******************************************************************************/
static inline int
readArrayOfInteger( FILE         *stream,
                    spm_int_t  n,
                    spm_int_t *array )
{
    long tmp1, tmp2, tmp3, tmp4;
    spm_int_t i;

    /* Read 4 by 4 */
    for (i=0; i<(n-3); i+=4)
    {
        if (4 != fscanf(stream, "%ld %ld %ld %ld", &tmp1, &tmp2, &tmp3, &tmp4)){
            spm_print_error("spmLoad: Wrong input format");
            return SPM_ERR_FILE;
        }

        array[i  ] = (spm_int_t)tmp1;
        array[i+1] = (spm_int_t)tmp2;
        array[i+2] = (spm_int_t)tmp3;
        array[i+3] = (spm_int_t)tmp4;
    }

    assert( n-i < 4 );
    switch ( n - i )
    {
    case 3:
        if (3 != fscanf(stream, "%ld %ld %ld", &tmp1, &tmp2, &tmp3)){
            spm_print_error("spmLoad: Wrong input format");
            return SPM_ERR_FILE;
        }

        array[i  ] = (spm_int_t)tmp1;
        array[i+1] = (spm_int_t)tmp2;
        array[i+2] = (spm_int_t)tmp3;
        break;
    case 2:
        if (2 != fscanf(stream, "%ld %ld", &tmp1, &tmp2)){
            spm_print_error("spmLoad: Wrong input format");
            return SPM_ERR_FILE;
        }

        array[i  ] = (spm_int_t)tmp1;
        array[i+1] = (spm_int_t)tmp2;
        break;
    case 1:
        if (1 != fscanf(stream, "%ld", &tmp1)){
            spm_print_error("spmLoad: Wrong input format");
            return SPM_ERR_FILE;
        }

        array[i  ] = (spm_int_t)tmp1;
        break;
    }

    return SPM_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_io
 *
 * @brief  Read an array of 64bits complex.
 *
 *******************************************************************************
 *
 * @param[in] stream
 *          The opened file in which the spm is stored.
 *
 * @param[in] n
 *          Number of elements.
 *
 * @param[out] array
 *          Array of results.
 *
 ********************************************************************************
 *
 * @retval SPM_SUCCESS if the read happened successfully,
 * @retval SPM_ERR_FILE if the input format is incorrect.
 *
 *******************************************************************************/
static inline int
readArrayOfComplex64( FILE               *stream,
                      spm_int_t        n,
                      spm_complex64_t *array )
{
    double tmp1, tmp2, tmp3, tmp4;
    double tmp5, tmp6, tmp7, tmp8;
    spm_int_t i;

    /* Read 4 by 4 */
    for (i=0; i<(n-3); i+=4)
    {
        if (8 != fscanf(stream, "%lg %lg %lg %lg %lg %lg %lg %lg",
                        &tmp1, &tmp2, &tmp3, &tmp4,
                        &tmp5, &tmp6, &tmp7, &tmp8)){
            spm_print_error("spmLoad: Wrong input format");
            return SPM_ERR_FILE;
        }
        array[i  ] = (spm_complex64_t)(tmp1 + I * tmp2);
        array[i+1] = (spm_complex64_t)(tmp3 + I * tmp4);
        array[i+2] = (spm_complex64_t)(tmp5 + I * tmp6);
        array[i+3] = (spm_complex64_t)(tmp7 + I * tmp8);
    }

    assert( n-i < 4 );
    switch ( n - i )
    {
    case 3:
        if (6 != fscanf(stream, "%lg %lg %lg %lg %lg %lg",
                        &tmp1, &tmp2, &tmp3, &tmp4,
                        &tmp5, &tmp6)){
            spm_print_error("spmLoad: Wrong input format");
            return SPM_ERR_FILE;
        }
        array[i  ] = (spm_complex64_t)(tmp1 + I * tmp2);
        array[i+1] = (spm_complex64_t)(tmp3 + I * tmp4);
        array[i+2] = (spm_complex64_t)(tmp5 + I * tmp6);
        break;

    case 2:
        if (4 != fscanf(stream, "%lg %lg %lg %lg",
                        &tmp1, &tmp2, &tmp3, &tmp4)){
            spm_print_error("spmLoad: Wrong input format");
            return SPM_ERR_FILE;
        }
        array[i  ] = (spm_complex64_t)(tmp1 + I * tmp2);
        array[i+1] = (spm_complex64_t)(tmp3 + I * tmp4);
        break;

    case 1:
        if (2 != fscanf(stream, "%lg %lg",
                        &tmp1, &tmp2)){
            spm_print_error("spmLoad: Wrong input format");
            return SPM_ERR_FILE;
        }
        array[i  ] = (spm_complex64_t)(tmp1 + I * tmp2);
        break;
    }

    return SPM_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_io
 *
 * @brief  Read an array of 32bits complex.
 *
 *******************************************************************************
 *
 * @param[in] stream
 *          The opened file in which the spm is stored.
 *
 * @param[in] n
 *          Number of elements.
 *
 * @param[out] array
 *          Array of results.
 *
 ********************************************************************************
 *
 * @retval SPM_SUCCESS if the read happened successfully,
 * @retval SPM_ERR_FILE if the input format is incorrect.
 *
 *******************************************************************************/
static inline int
readArrayOfComplex32( FILE               *stream,
                      spm_int_t        n,
                      spm_complex32_t *array )
{
    float tmp1, tmp2, tmp3, tmp4;
    float tmp5, tmp6, tmp7, tmp8;
    spm_int_t i;

    /* Read 4 by 4 */
    for (i=0; i<(n-3); i+=4)
    {
        if (8 != fscanf(stream, "%g %g %g %g %g %g %g %g",
                        &tmp1, &tmp2, &tmp3, &tmp4,
                        &tmp5, &tmp6, &tmp7, &tmp8)){
            spm_print_error("spmLoad: Wrong input format");
            return SPM_ERR_FILE;
        }
        array[i  ] = (spm_complex32_t)(tmp1 + I * tmp2);
        array[i+1] = (spm_complex32_t)(tmp3 + I * tmp4);
        array[i+2] = (spm_complex32_t)(tmp5 + I * tmp6);
        array[i+3] = (spm_complex32_t)(tmp7 + I * tmp8);
    }

    assert( n-i < 4 );
    switch ( n - i )
    {
    case 3:
        if (6 != fscanf(stream, "%g %g %g %g %g %g",
                        &tmp1, &tmp2, &tmp3, &tmp4,
                        &tmp5, &tmp6)){
            spm_print_error("spmLoad: Wrong input format");
            return SPM_ERR_FILE;
        }
        array[i  ] = (spm_complex32_t)(tmp1 + I * tmp2);
        array[i+1] = (spm_complex32_t)(tmp3 + I * tmp4);
        array[i+2] = (spm_complex32_t)(tmp5 + I * tmp6);
        break;

    case 2:
        if (4 != fscanf(stream, "%g %g %g %g",
                        &tmp1, &tmp2, &tmp3, &tmp4)){
            spm_print_error("spmLoad: Wrong input format");
            return SPM_ERR_FILE;
        }
        array[i  ] = (spm_complex32_t)(tmp1 + I * tmp2);
        array[i+1] = (spm_complex32_t)(tmp3 + I * tmp4);
        break;

    case 1:
        if (2 != fscanf(stream, "%g %g",
                        &tmp1, &tmp2)){
            spm_print_error("spmLoad: Wrong input format");
            return SPM_ERR_FILE;
        }
        array[i  ] = (spm_complex32_t)(tmp1 + I * tmp2);
        break;
    }

    return SPM_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_io
 *
 * @brief  Read an array of double.
 *
 *******************************************************************************
 *
 * @param[in] stream
 *          The opened file in which the spm is stored.
 *
 * @param[in] n
 *          Number of elements.
 *
 * @param[out] array
 *          Array of results.
 *
 ********************************************************************************
 *
 * @retval SPM_SUCCESS if the read happened successfully,
 * @retval SPM_ERR_FILE if the input format is incorrect.
 *
 *******************************************************************************/
static inline int
readArrayOfDouble( FILE         *stream,
                   spm_int_t  n,
                   double       *array )
{
    double tmp1, tmp2, tmp3, tmp4;
    spm_int_t i;

    /* Read 4 by 4 */
    for (i=0; i<(n-3); i+=4)
    {
        if (4 != fscanf(stream, "%lg %lg %lg %lg",
                        &tmp1, &tmp2, &tmp3, &tmp4)){
            spm_print_error("spmLoad: Wrong input format");
            return SPM_ERR_FILE;
        }
        array[i  ] = (double)(tmp1);
        array[i+1] = (double)(tmp2);
        array[i+2] = (double)(tmp3);
        array[i+3] = (double)(tmp4);
    }

    assert( n-i < 4 );
    switch ( n - i )
    {
    case 3:
        if (1 != fscanf(stream, "%lg %lg %lg",
                        &tmp1, &tmp2, &tmp3)){
            spm_print_error("spmLoad: Wrong input format");
            return SPM_ERR_FILE;
        }
        array[i  ] = (double)(tmp1);
        array[i+1] = (double)(tmp2);
        array[i+2] = (double)(tmp3);
        break;

    case 2:
        if (2 != fscanf(stream, "%lg %lg",
                        &tmp1, &tmp2)){
            spm_print_error("spmLoad: Wrong input format");
            return SPM_ERR_FILE;
        }
        array[i  ] = (double)(tmp1);
        array[i+1] = (double)(tmp2);
        break;

    case 1:
        if (1 != fscanf(stream, "%lg",
                        &tmp1)){
            spm_print_error("spmLoad: Wrong input format");
            return SPM_ERR_FILE;
        }
        array[i  ] = (double)(tmp1);
        break;
    }

    return SPM_SUCCESS;
}


/**
 *******************************************************************************
 *
 * @ingroup spm_dev_io
 *
 * @brief  Read an array of float.
 *
 *******************************************************************************
 *
 * @param[in] stream
 *          The opened file in which the spm is stored.
 *
 * @param[in] n
 *          Number of elements.
 *
 * @param[out] array
 *          Array of results.
 *
 ********************************************************************************
 *
 * @retval SPM_SUCCESS if the read happened successfully,
 * @retval SPM_ERR_FILE if the input format is incorrect.
 *
 *******************************************************************************/
static inline int
readArrayOfFloat( FILE         *stream,
                  spm_int_t  n,
                  float        *array )
{
    float tmp1, tmp2, tmp3, tmp4;
    spm_int_t i;

    /* Read 4 by 4 */
    for (i=0; i<(n-3); i+=4)
    {
        if (4 != fscanf(stream, "%g %g %g %g",
                        &tmp1, &tmp2, &tmp3, &tmp4)){
            spm_print_error("spmLoad: Wrong input format");
            return SPM_ERR_FILE;
        }
        array[i  ] = (float)(tmp1);
        array[i+1] = (float)(tmp2);
        array[i+2] = (float)(tmp3);
        array[i+3] = (float)(tmp4);
    }

    assert( n-i < 4 );
    switch ( n - i )
    {
    case 3:
        if (3 != fscanf(stream, "%g %g %g",
                        &tmp1, &tmp2, &tmp3)){
            spm_print_error("spmLoad: Wrong input format");
            return SPM_ERR_FILE;
        }
        array[i  ] = (float)(tmp1);
        array[i+1] = (float)(tmp2);
        array[i+2] = (float)(tmp3);
        break;

    case 2:
        if (2 != fscanf(stream, "%g %g",
                        &tmp1, &tmp2)){
            spm_print_error("spmLoad: Wrong input format");
            return SPM_ERR_FILE;
        }
        array[i  ] = (float)(tmp1);
        array[i+1] = (float)(tmp2);
        break;

    case 1:
        if (1 != fscanf(stream, "%g", &tmp1)){
            spm_print_error("spmLoad: Wrong input format");
            return SPM_ERR_FILE;
        }
        array[i  ] = (float)(tmp1);
        break;
    }

    return SPM_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup spm
 *
 * @brief Load the spm structure from a file (internal format).
 *
 *
 *******************************************************************************
 *
 * @param[inout] spm
 *          On entry, an allocated spm data structure.
 *          On exit, the spm filled with the information read in the file.
 *
 * @param[in] infile
 *          The opened file in which the spm is stored. If infile == NULL,
 *          matrix.spm is opened.
 *
 *******************************************************************************
 *
 * @retval SPM_SUCCESS if the load happened successfully,
 * @retval SPM_ERR_FILE if the input format is incorrect.
 *
 *******************************************************************************/
int
spmLoad( spmatrix_t  *spm,
         FILE          *infile )
{
    spm_int_t colsize=0, rowsize=0;
    char line[256], *test;
    int rc = SPM_SUCCESS;
    int local_stream = 0;

    if ( infile == NULL ) {
        infile = fopen( "matrix.spm", "r" );

        if ( infile == NULL ) {
            spm_print_error( "spmLoad: Impossible to open the file matrix.spm\n");
            return SPM_ERR_FILE;
        }

        local_stream = 1;
    }

    /*
     * Skip comments
     */
    do {
        test = fgets( line, 256, infile );
        if ( test != line ) {
            return SPM_ERR_FILE;
        }
    }
    while( line[0] == '#' );

    /*
     * Read header
     */
    {
        int version, mtxtype, flttype, fmttype, dof, layout;
        long gN, n, nnz, nnzexp;

        if ( 10 != sscanf( line, "%d %d %d %d %ld %ld %ld %d %ld %d\n",
                           &version, &mtxtype, &flttype, &fmttype,
                           &gN, &n, &nnz, &dof, &nnzexp, &layout ) ) {
            return SPM_ERR_FILE;
        }

        /* Handle only version 1 for now */
        if (version != 1) {
            return SPM_ERR_BADPARAMETER;
        }

        spm->mtxtype = (spm_mtxtype_t)mtxtype;
        spm->flttype = (spm_coeftype_t)flttype;
        spm->fmttype = (spm_fmttype_t)fmttype;
        spm->gN      = gN;
        spm->n       = n;
        spm->nnz     = nnz;
        spm->dof     = dof;
        spm->layout  = (spm_layout_t)layout;

        spmUpdateComputedFields( spm );

        assert( nnzexp == spm->nnzexp );
        assert( spm->gN == gN );
    }

    switch(spm->fmttype){
    case SpmCSC:
        colsize = spm->n + 1;
        rowsize = spm->nnz;
        break;
    case SpmCSR:
        colsize = spm->nnz;
        rowsize = spm->n + 1;
        break;
    case SpmIJV:
        colsize = spm->nnz;
        rowsize = spm->nnz;
        break;
    }

    /*
     * Read colptr
     */
    spm->colptr = malloc( colsize * sizeof(spm_int_t) );
    rc = readArrayOfInteger( infile, colsize, spm->colptr );
    if (rc != SPM_SUCCESS ) {
        return rc;
    }

    /*
     * Read rowptr
     */
    spm->rowptr = malloc( rowsize * sizeof(spm_int_t) );
    rc = readArrayOfInteger( infile, rowsize, spm->rowptr );
    if (rc != SPM_SUCCESS ) {
        return rc;
    }

    /*
     * Read dofs
     */
    if ( spm->dof > 0 ) {
        spm->dofs = NULL;
    }
    else {
        spm->dofs = malloc( (spm->n+1) * sizeof(spm_int_t) );
        rc = readArrayOfInteger( infile, spm->n+1, spm->dofs );
        if (rc != SPM_SUCCESS ) {
            return rc;
        }
    }

    /*
     * Read dofs
     */
    if ( spm->gN == spm->n ) {
        spm->loc2glob = NULL;
    }
    else {
        spm->loc2glob = malloc( spm->n * sizeof(spm_int_t) );
        rc = readArrayOfInteger( infile, spm->n, spm->dofs );
        if (rc != SPM_SUCCESS ) {
            return rc;
        }
    }

    /*
     * Read values
     */
    if (spm->flttype == SpmPattern ) {
        spm->values = NULL;
    }
    else {
        spm->values = malloc( spm->nnzexp * spm_size_of( spm->flttype ) );
    }

    switch( spm->flttype ) {
    case SpmPattern:
        break;
    case SpmFloat:
        rc = readArrayOfFloat( infile, spm->nnzexp, spm->values );
        break;
    case SpmDouble:
        rc = readArrayOfDouble( infile, spm->nnzexp, spm->values );
        break;
    case SpmComplex32:
        rc = readArrayOfComplex32( infile, spm->nnzexp, spm->values );
        break;
    case SpmComplex64:
        rc = readArrayOfComplex64( infile, spm->nnzexp, spm->values );
        break;
    }

    if (local_stream) {
        fclose(infile);
    }

    return rc;
}

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_io
 *
 * @brief write an array of 64bits complex.
 *
 *******************************************************************************
 *
 * @param[in] outfile
 *          TODO
 *
 * @param[in] n
 *          numbers of elements.
 *
 * @param[in] array
 *          array to write.
 *
 *******************************************************************************
 *
 * @return   SPM_SUCCESS if the write happened successfully.
 *
 *******************************************************************************/
static inline int
writeArrayOfComplex64( FILE                     *outfile,
                       spm_int_t              n,
                       const spm_complex64_t *array )
{
    spm_int_t i;

    /* Write 4 by 4 */
    for (i=0; i<n; i++)
    {
        fprintf(outfile, "%e %e ", creal(array[i]), cimag(array[i]));
        if (i%4 == 3) fprintf(outfile, "\n");
    }
    if ((i-1)%4 !=3) fprintf(outfile, "\n");

    return SPM_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_io
 *
 * @brief write an array of 32bits complex.
 *
 *******************************************************************************
 *
 * @param[in] outfile
 *          TODO
 *
 * @param[in] n
 *          numbers of elements.
 *
 * @param[in] array
 *          array to write.
 *
 *******************************************************************************
 *
 * @return   SPM_SUCCESS if the write happened successfully.
 *
 *******************************************************************************/
static inline int
writeArrayOfComplex32( FILE                     *outfile,
                       spm_int_t              n,
                       const spm_complex32_t *array )
{
    spm_int_t i;

    /* Write 4 by 4 */
    for (i=0; i<n; i++)
    {
        fprintf(outfile, "%e %e ", crealf(array[i]), cimagf(array[i]));
        if (i%4 == 3) fprintf(outfile, "\n");
    }
    if ((i-1)%4 !=3) fprintf(outfile, "\n");

    return SPM_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_io
 *
 * @brief write an array of double.
 *
 *******************************************************************************
 *
 * @param[in] outfile
 *          TODO
 *
 * @param[in] n
 *          numbers of elements.
 *
 * @param[in] array
 *          array to write.
 *
 *******************************************************************************
 *
 * @return   SPM_SUCCESS if the write happened successfully.
 *
 *******************************************************************************/
static inline int
writeArrayOfDouble( FILE         *outfile,
                    spm_int_t  n,
                    const double *array )
{
    spm_int_t i;

    /* Write 4 by 4 */
    for (i=0; i<n; i++)
    {
        fprintf(outfile, "%e ", array[i]);
        if (i%4 == 3) fprintf(outfile, "\n");
    }
    if ((i-1)%4 !=3) fprintf(outfile, "\n");

    return SPM_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_io
 *
 * @brief write an array of float.
 *
 *******************************************************************************
 *
 * @param[in] outfile
 *          TODO
 *
 * @param[in] n
 *          numbers of elements.
 *
 * @param[in] array
 *          array to write.
 *
 *******************************************************************************
 *
 * @return   SPM_SUCCESS if the write happened successfully.
 *
 *
 *******************************************************************************/
static inline int
writeArrayOfFloat( FILE         *outfile,
                   spm_int_t  n,
                   const float  *array )
{
    spm_int_t i;

    /* Write 4 by 4 */
    for (i=0; i<n; i++)
    {
        fprintf(outfile, "%e ", array[i]);
        if (i%4 == 3) fprintf(outfile, "\n");
    }
    if ((i-1)%4 !=3) fprintf(outfile, "\n");

    return SPM_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup spm
 *
 * @brief Save the spm structure into a file (internal format).
 *
 *
 *******************************************************************************
 *
 * @param[in] spm
 *          The sparse matrix to write into the file.
 *
 * @param[in] outfile
 *          The opened file in which to store the spm. If outfile == NULL, data
 *          is saved into matrix.spm file.
 *
 ********************************************************************************
 *
 * @return  SPM_SUCCESS if the save happened successfully.
 *
 *******************************************************************************/
int
spmSave( const spmatrix_t *spm,
         FILE               *outfile )
{
    spm_int_t i, colsize, rowsize;
    int local_stream = 0;

    if ( outfile == NULL ) {
        outfile = fopen( "matrix.spm", "w" );
        if ( outfile == NULL ) {
            spm_print_error( "spmSave: Impossible to open the file matrix.spm\n");
            return SPM_ERR_FILE;
        }

        local_stream = 1;
    }

    /*
     * Write header
     */
    fprintf( outfile,
             "# version mtxtype flttype fmttype gN n nnz dof nnzexp layout\n"
             "%d %d %d %d %ld %ld %ld %d %ld %d\n",
             1, spm->mtxtype, spm->flttype, spm->fmttype,
             (long)spm->gN, (long)spm->n, (long)spm->nnz,
             (int)spm->dof, (long)spm->nnzexp, spm->layout );

    switch(spm->fmttype){
    case SpmCSC:
        colsize = spm->n + 1;
        rowsize = spm->nnz;
        break;
    case SpmCSR:
        colsize = spm->nnz;
        rowsize = spm->n + 1;
        break;
    case SpmIJV:
        colsize = spm->nnz;
        rowsize = spm->nnz;
        break;
    default:
        colsize = 0;
        rowsize = 0;
    }

    /*
     * Write colptr
     */
    for (i=0; i<colsize; i++)
    {
        fprintf(outfile, "%ld ", (long)spm->colptr[i]);
        if (i%4 == 3) fprintf(outfile, "\n");
    }
    if ((i-1)%4 !=3) fprintf(outfile, "\n");

    /*
     * Write rowptr
     */
    for (i=0; i<rowsize; i++)
    {
        fprintf(outfile, "%ld ", (long)spm->rowptr[i]);
        if (i%4 == 3) fprintf(outfile, "\n");
    }
    if ((i-1)%4 !=3) fprintf(outfile, "\n");

    /*
     * Write dofs
     */
    if ( spm->dof <= 0 ) {
        for (i=0; i<spm->n+1; i++)
        {
            fprintf(outfile, "%ld ", (long)spm->dofs[i]);
            if (i%4 == 3) fprintf(outfile, "\n");
        }
        if ((i-1)%4 !=3) fprintf(outfile, "\n");
    }

    /*
     * Write loc2glob
     */
    if ( spm->n != spm->gN ) {
        for (i=0; i<spm->n; i++)
        {
            fprintf(outfile, "%ld ", (long)spm->loc2glob[i]);
            if (i%4 == 3) fprintf(outfile, "\n");
        }
        if ((i-1)%4 !=3) fprintf(outfile, "\n");
    }

    /*
     * Write values
     */
    switch( spm->flttype ) {
    case SpmPattern:
        break;
    case SpmFloat:
        writeArrayOfFloat( outfile, spm->nnzexp, spm->values );
        break;
    case SpmDouble:
        writeArrayOfDouble( outfile, spm->nnzexp, spm->values );
        break;
    case SpmComplex32:
        writeArrayOfComplex32( outfile, spm->nnzexp, spm->values );
        break;
    case SpmComplex64:
        writeArrayOfComplex64( outfile, spm->nnzexp, spm->values );
        break;
    }

    if (local_stream) {
        fclose(outfile);
    }
    return SPM_SUCCESS;
}
