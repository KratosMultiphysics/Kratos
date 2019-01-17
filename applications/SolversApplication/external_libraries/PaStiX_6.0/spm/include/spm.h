/**
 *
 * @file spm.h
 *
 * SParse Matrix package header.
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
 * @addtogroup spm
 * @{
 *   @brief Describe all the internals routines of the SParse Matrix package.
 *
 *   This library provides a set of subroutines to manipulate sparse matrices in
 *   different format such as compressed sparse column (CSC), compressed sparse
 *   row (CSR), or coordinate (IJV) with single or multiple degrees of freedom
 *   per unknown. It provides basic BLAS 1 and BLAS 2 functions for those
 *   matrices, as well as norms computations and converter tools.
 *
 **/
#ifndef _spm_h_
#define _spm_h_

#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include "spm_config.h"
#include "spm_const.h"
#include "spm_datatypes.h"

/**
 *
 * @brief The sparse matrix data structure
 *
 * This structure describes matrices with different characteristics that can be useful to any
 * solver:
 *     - the storage format (SpmCSC, SpmCSR or SpmIJV)
 *     - the properties (SpmGeneral, SpmHermitian, SpmSymmetric)
 *     - the base value (0 in C or 1 in Fortran)
 *
 * It is also possible to describe a matrix with constant or variable degrees of freedom.
 *
 */
typedef struct spmatrix_s {
    spm_mtxtype_t  mtxtype; /**< Matrix structure: SpmGeneral, SpmSymmetric
                                 or SpmHermitian.                                               */
    spm_coeftype_t flttype; /**< values datatype: SpmPattern, SpmFloat, SpmDouble,
                                 SpmComplex32 or SpmComplex64                                   */
    spm_fmttype_t  fmttype; /**< Matrix storage format: SpmCSC, SpmCSR, SpmIJV                  */

    spm_int_t      gN;      /**< Global number of vertices in the compressed graph (Computed)   */
    spm_int_t      n;       /**< Local number of vertices in the compressed graph               */
    spm_int_t      gnnz;    /**< Global number of non zeroes in the compressed graph (Computed) */
    spm_int_t      nnz;     /**< Local number of non zeroes in the compressed graph             */

    spm_int_t      gNexp;   /**< Global number of vertices in the compressed graph (Computed)   */
    spm_int_t      nexp;    /**< Local number of vertices in the compressed graph (Computed)    */
    spm_int_t      gnnzexp; /**< Global number of non zeroes in the compressed graph (Computed) */
    spm_int_t      nnzexp;  /**< Local number of non zeroes in the compressed graph (Computed)  */

    spm_int_t      dof;     /**< Number of degrees of freedom per unknown,
                                 if > 0, constant degree of freedom
                                 otherwise, irregular degree of freedom (refer to dofs)         */
    spm_int_t     *dofs;    /**< Array of the first column of each element in the
                                 expanded matrix [+baseval]                                     */
    spm_layout_t   layout;  /**< SpmColMajor, or SpmRowMajor                                    */

    spm_int_t     *colptr;  /**< List of indirections to rows for each vertex [+baseval]        */
    spm_int_t     *rowptr;  /**< List of edges for each vertex [+baseval]                       */
    spm_int_t     *loc2glob;/**< Corresponding numbering from local to global [+baseval]        */
    void          *values;  /**< Values stored in the matrix                                    */
} spmatrix_t;

/**
 * @name SPM basic subroutines
 * @{
 */
void spmInit ( spmatrix_t *spm );
void spmAlloc( spmatrix_t *spm );
void spmExit ( spmatrix_t *spm );

spmatrix_t *spmCopy( const spmatrix_t *spm );
void        spmBase( spmatrix_t *spm, int baseval );
spm_int_t   spmFindBase( const spmatrix_t *spm );
int         spmConvert( int ofmttype, spmatrix_t *ospm );
void        spmUpdateComputedFields( spmatrix_t *spm );
void        spmGenFakeValues( spmatrix_t *spm );

/**
 * @}
 * @name SPM BLAS subroutines
 * @{
 */
double spmNorm( spm_normtype_t ntype, const spmatrix_t *spm );
int    spmMatVec( spm_trans_t trans, double alpha, const spmatrix_t *spm, const void *x, double beta, void *y );
int    spmMatMat( spm_trans_t trans, spm_int_t n,
                  double alpha, const spmatrix_t *A,
                                const void *B, spm_int_t ldb,
                  double beta,        void *C, spm_int_t ldc );
void   spmScalMatrix( double alpha, spmatrix_t *spm );
void   spmScalVector( spm_coeftype_t flt, double alpha, spm_int_t n, void *x, spm_int_t incx );

/**
 * @}
 * @name SPM subroutines to check format
 * @{
 */
int       spmSort( spmatrix_t *spm );
spm_int_t spmMergeDuplicate( spmatrix_t *spm );
spm_int_t spmSymmetrize( spmatrix_t *spm );
int       spmCheckAndCorrect( const spmatrix_t *spm_in, spmatrix_t *spm_out );

/**
 * @}
 * @name SPM subroutines to check factorization/solve
 * @{
 */
int spmGenRHS( spm_rhstype_t     type,
               spm_int_t         nrhs,
               const spmatrix_t *spm,
               void * const      x,
               spm_int_t         ldx,
               void * const      b,
               spm_int_t         ldb );
int spmCheckAxb( double            eps,
                 spm_int_t         nrhs,
                 const spmatrix_t *spm,
                 void * const      x0,
                 spm_int_t         ldx0,
                 void * const      b,
                 spm_int_t         ldb,
                 const void *      x,
                 spm_int_t         ldx );

/**
 * @}
 * @name SPM subroutines to manipulate integers arrays
 * @{
 */
spm_int_t *spmIntConvert( spm_int_t n, int *input );
void       spmIntSort1Asc1( void *const pbase, const spm_int_t n );
void       spmIntSort2Asc1( void *const pbase, const spm_int_t n );
void       spmIntSort2Asc2( void *const pbase, const spm_int_t n );
void       spmIntMSortIntAsc( void **const pbase, const spm_int_t n );

/**
 * @}
 * @name SPM IO subroutines
 * @{
 */
int spmLoad( spmatrix_t *spm, FILE *infile );
int spmSave( const spmatrix_t *spm, FILE *outfile );

/**
 * @}
 * @name SPM driver
 * @{
 */
int spmReadDriver( spm_driver_t driver,
                   const char  *filename,
                   spmatrix_t  *spm );

int spmParseLaplacianInfo( const char *    filename,
                           spm_coeftype_t *flttype,
                           spm_int_t *     dim1,
                           spm_int_t *     dim2,
                           spm_int_t *     dim3,
                           double *        alpha,
                           double *        beta );

/**
 * @}
 * @name SPM debug subroutines
 * @{
 */
void *      spm2Dense( const spmatrix_t *spm );
void        spmPrint( const spmatrix_t *spm, FILE *f );
void        spmPrintInfo( const spmatrix_t *spm, FILE *f );
spmatrix_t *spmExpand( const spmatrix_t *spm );
spmatrix_t *spmDofExtend( const spmatrix_t *spm, const int type, const int dof );

/**
 * @}
 */

/**
 * @}
 */

/**
 * @name SPM dev printing subroutines
 * @{
 *
 */

/**
 * @ingroup spm_dev_print
 * @brief Subroutines to print one element of an spm structure
 *
 * @param[in] f Pointer to the file
 * @param[in] i Row index of the element
 * @param[in] j Column index of the element
 * @param[in] A Value of the element A|i,j]
 *
 * Double complex case
 *
 */
static inline void
z_spmPrintElt( FILE *f, spm_int_t i, spm_int_t j, spm_complex64_t A )
{
    fprintf( f, "%ld %ld %e %e\n", (long)i, (long)j, creal( A ), cimag( A ) );
}

/**
 * @copydoc z_spmPrintElt
 * @details Single complex case
 */
static inline void
c_spmPrintElt( FILE *f, spm_int_t i, spm_int_t j, spm_complex32_t A )
{
    fprintf( f, "%ld %ld %e %e\n", (long)i, (long)j, crealf( A ), cimagf( A ) );
}
/**
 * @copydoc z_spmPrintElt
 * @details Double real case
 */
static inline void
d_spmPrintElt( FILE *f, spm_int_t i, spm_int_t j, double A )
{
    fprintf( f, "%ld %ld %e\n", (long)i, (long)j, A );
}
/**
 * @copydoc z_spmPrintElt
 * @details Single real case
 */
static inline void
s_spmPrintElt( FILE *f, spm_int_t i, spm_int_t j, float A )
{
    fprintf( f, "%ld %ld %e\n", (long)i, (long)j, A );
}
/**
 * @copydoc z_spmPrintElt
 * @details Pattern case
 *
 * @remark: uses a macro to avoid accessing A that would generate segfault.
 */
#define p_spmPrintElt( f, i, j, A )                                                                \
    {                                                                                              \
        fprintf( f, "%ld %ld\n", (long)( i ), (long)( j ) );                                       \
    }

/**
 * @}
 */
#endif /* _spm_h_ */
