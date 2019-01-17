/**
 *
 * @file laplacian.c
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
#include "drivers/laplacian.h"

/**
 *******************************************************************************
 *
 * @ingroup spm_driver
 *
 * laplacian_usage - Print the usage information to generate correct Laplacian
 * matrices.
 *
 *******************************************************************************/
static inline void
laplacian_usage(void)
{
    fprintf(stderr,
            "Usage: genLaplacian( \"[<type>:]<dim1>[:<dim2>[:<dim3>[:<alpha>[:<beta>]]]]\" )\n"
            "   Generate a Laplacian matrix M, of the form alpha * D - beta * A,\n"
            "   where D is the degree matrix, and A the adjacency matrix.\n"
            "   <type> p = pattern only\n"
            "          s = real simple\n"
            "          d = real double [default]\n"
            "          c = complex simple\n"
            "          z = complex double\n"
            "   <dim1> size of the first dimension of the laplacian\n"
            "   <dim2> size of the second dimension of the laplacian\n"
            "   <dim3> size of the third dimension of the laplacian\n"
            "   Example:\n"
            "     genLaplacian( \"z:10:20\" )        generates a 2D complex double laplacian matrix of size 200.\n"
            "     genLaplacian( \"10:1:10:2.:0.5\" ) generates a 2D real double laplacian matrix of size 100 where M = 2. * D - 0.5 * A.\n"
            "     genLaplacian( \"s:10\" )           generates a 1D real single laplacian matrix of size 10.\n"
            );
}

/**
 *******************************************************************************
 *
 * @ingroup spm
 *
 * @brief Parse information given through the filename string to configure the
 * laplacian matrix to generate.
 *
 * The laplacian will be of size dim1 * dim2 * dim3, and will be equal to
 *     \[ M = \alpha * D - \beta * A \]
 *
 * where D is the degree matrix, and A the adjacency matrix.
 *
 *******************************************************************************
 *
 * @param[in] filename
 *          Configuration string of the Laplacian. See laplacian_usage() for
 *          more information.
 *
 * @param[out] flttype
 *          The floating type of the elements in the matrix.
 *
 * @param[out] dim1
 *          The first dimension of the laplacian
 *
 * @param[out] dim2
 *          The second dimension of the laplacian
 *
 * @param[out] dim3
 *          The third dimension of the laplacian
 *
 * @param[out] alpha
 *          The alpha coefficient for the degree matrix
 *
 * @param[out] beta
 *          The beta coefficient for the adjacency matrix
 *
 *******************************************************************************
 *
 * @retval SPM_SUCCESS if the matrix has been generated successfully
 * @retval SPM_ERR_BADPARAMETER if the configuration string is incorrect
 *
 *******************************************************************************/
int
spmParseLaplacianInfo( const char     *filename,
                       spm_coeftype_t *flttype,
                       spm_int_t      *dim1,
                       spm_int_t      *dim2,
                       spm_int_t      *dim3,
                       double         *alpha,
                       double         *beta )
{
    double val1, val2;
    long tmp1, tmp2, tmp3;

    *alpha = 1.;
    *beta = 1.;

    /* Look for the datatype */
    {
        char flt;
        char *tmpf = strndup( filename, 256 );

        if ( sscanf( filename, "%c:%254s", &flt, tmpf ) == 2 ) {
            filename += 2;
            switch( flt ){
            case 'Z':
            case 'z':
                *flttype = SpmComplex64;
                break;

            case 'C':
            case 'c':
                *flttype = SpmComplex32;
                break;

            case 'D':
            case 'd':
                *flttype = SpmDouble;
                break;

            case 'S':
            case 's':
                *flttype = SpmFloat;
                break;

            case 'P':
            case 'p':
                *flttype = SpmPattern;
                break;

            case '1':
            case '2':
            case '3':
            case '4':
            case '5':
            case '6':
            case '7':
            case '8':
            case '9':
                *flttype = SpmDouble;
                /*
                 * The first dimension is only one character long so we come
                 * back to the beginning of the string
                 */
                filename -= 2;
                break;

            default:
                laplacian_usage();
                return SPM_ERR_BADPARAMETER;
            }
        }
        else {
            *flttype = SpmDouble;
        }

        free(tmpf);
    }

    /* Scan the dimensions */
    *dim1 = *dim2 = *dim3 = 1;

    if ( sscanf( filename, "%ld:%ld:%ld:%lf:%lf", &tmp1, &tmp2, &tmp3, &val1, &val2 ) == 5 ) {
        *dim1 = (spm_int_t)tmp1;
        *dim2 = (spm_int_t)tmp2;
        *dim3 = (spm_int_t)tmp3;
        *alpha = val1;
        *beta  = val2;
    }
    else if ( sscanf( filename, "%ld:%ld:%ld:%lf", &tmp1, &tmp2, &tmp3, &val1 ) == 4 ) {
        *dim1 = (spm_int_t)tmp1;
        *dim2 = (spm_int_t)tmp2;
        *dim3 = (spm_int_t)tmp3;
        *alpha = val1;
    }
    else if ( sscanf( filename, "%ld:%ld:%ld", &tmp1, &tmp2, &tmp3 ) == 3 ) {
        *dim1 = (spm_int_t)tmp1;
        *dim2 = (spm_int_t)tmp2;
        *dim3 = (spm_int_t)tmp3;
    }
    else if ( sscanf( filename, "%ld:%ld", &tmp1, &tmp2 ) == 2 ) {
        *dim1 = (spm_int_t)tmp1;
        *dim2 = (spm_int_t)tmp2;
    }
    else if ( sscanf( filename, "%ld", &tmp1 ) == 1 ) {
        *dim1 = (spm_int_t)tmp1;
    }
    else {
        laplacian_usage();
        return SPM_ERR_BADPARAMETER;
    }

    /* One of the dimension was set to 0 */
    if ( (*dim1 == 0) || (*dim2 == 0) || (*dim3 == 0) ) {
        laplacian_usage();
        return SPM_ERR_BADPARAMETER;
    }

    return SPM_SUCCESS;
}

static void (*laplacian_7points[6])(spmatrix_t *, spm_int_t, spm_int_t, spm_int_t, spm_fixdbl_t, spm_fixdbl_t) =
{
    p_spmLaplacian_7points,
    NULL,
    s_spmLaplacian_7points,
    d_spmLaplacian_7points,
    c_spmLaplacian_7points,
    z_spmLaplacian_7points
};

static void (*laplacian_27points[6])(spmatrix_t *, spm_int_t, spm_int_t, spm_int_t, spm_fixdbl_t, spm_fixdbl_t) =
{
    p_spmLaplacian_27points,
    NULL,
    s_spmLaplacian_27points,
    d_spmLaplacian_27points,
    c_spmLaplacian_27points,
    z_spmLaplacian_27points
};

/**
 *******************************************************************************
 *
 * @ingroup spm_driver
 *
 * genLaplacian - Generate a Laplacian of size spm->n
 *
 *******************************************************************************
 *
 * @param[in] filename
 *          Configuration string of the Laplacian.
 *          [<type>:]<dim1>[:<dim2>[:<dim3>]]
 *             <type> p = pattern only\n"
 *                    s = real simple\n"
 *                    d = real double [default]\n"
 *                    c = complex simple\n"
 *                    z = complex double\n"
 *             <dim1> size of the first dimension of the 1D|2D|3D laplacian\n"
 *             <dim2> size of the second dimension of the 2D|3D laplacian\n"
 *             <dim3> size of the third dimension of the 3D laplacian\n"
 *
 * @param[inout] spm
 *          At start, an allocated spm structure.
 *          At exit, contains a laplacian matrix in the spm format.
 *
 *******************************************************************************
 *
 * @return
 *      \retval SPM_SUCCESS if the matrix has been generated successfully
 *      \retval SPM_ERR_BADPARAMETER if the configuration string is incorrect
 *
 *******************************************************************************/
int
genLaplacian( const char    *filename,
              spmatrix_t  *spm )
{
    spm_coeftype_t flttype;
    spm_int_t dim1, dim2, dim3;
    double alpha = 1.;
    double beta = 1.;
    int rc;

    rc = spmParseLaplacianInfo(filename, &flttype, &dim1, &dim2, &dim3, &alpha, &beta );
    if (rc != SPM_SUCCESS)
        return rc;

    spm->flttype = flttype;
    spm->n = dim1 * dim2 * dim3;

    laplacian_7points[spm->flttype](spm, dim1, dim2, dim3, alpha, beta);

    return SPM_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup spm_driver
 *
 * genExtendedLaplacian - Generate a extended Laplacian of size spm->n
 *
 *******************************************************************************
 *
 * @param[in] filename
 *          Configuration string of the Laplacian.
 *          [<type>:]<dim1>[:<dim2>[:<dim3>]]
 *             <type> p = pattern only
 *                    s = real simple
 *                    d = real double [default]
 *                    c = complex simple
 *                    z = complex double
 *             <dim1> size of the first dimension of the 1D|2D|3D laplacian
 *             <dim2> size of the second dimension of the 2D|3D laplacian
 *             <dim3> size of the third dimension of the 3D laplacian
 *
 * @param[inout] spm
 *          At start, an allocated spm structure.
 *          At exit, contains a laplacian matrix in the spm format.
 *
 *******************************************************************************
 *
 * @return
 *      \retval SPM_SUCCESS if the matrix has been generated successfully
 *      \retval SPM_ERR_BADPARAMETER if the configuration string is incorrect
 *
 *******************************************************************************/
int
genExtendedLaplacian( const char    *filename,
                      spmatrix_t  *spm )
{
    spm_coeftype_t flttype;
    spm_int_t dim1, dim2, dim3;
    double alpha = 1.;
    double beta = 1.;
    int rc;

    rc = spmParseLaplacianInfo(filename, &flttype, &dim1, &dim2, &dim3, &alpha, &beta);
    if (rc != SPM_SUCCESS)
        return rc;

    spm->flttype = flttype;
    spm->n = dim1 * dim2 * dim3;

    laplacian_27points[spm->flttype](spm, dim1, dim2, dim3, alpha, beta);

    return SPM_SUCCESS;
}
