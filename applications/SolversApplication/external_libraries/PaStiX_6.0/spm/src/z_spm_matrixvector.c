/**
 *
 * @file z_spm_matrixvector.c
 *
 * SParse Matrix package matrix-vector multiplication routines.
 *
 * @copyright 2016-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.0.0
 * @author Matthieu Kuhn
 * @author Mathieu Faverge
 * @author Theophile Terraz
 * @date 2015-01-01
 *
 * @precisions normal z -> c d s
 **/
#include "common.h"
#include <lapacke.h>
#include <cblas.h>
#include "z_spm.h"

struct __spm_zmatvec_s;
typedef struct __spm_zmatvec_s __spm_zmatvec_t;
typedef spm_complex64_t (*__conj_fct_t)( spm_complex64_t );
typedef int (*__loop_fct_t)( const __spm_zmatvec_t * );

static inline spm_complex64_t
__fct_id( spm_complex64_t val ) {
    return val;
}

#if defined(PRECISION_c) || defined(PRECISION_z)
static inline spm_complex64_t
__fct_conj( spm_complex64_t val ) {
    return conj( val );
}
#endif

struct __spm_zmatvec_s {
    int                    follow_x;

    spm_int_t              baseval, n, nnz;

    spm_complex64_t        alpha;
    const spm_int_t       *rowptr;
    const spm_int_t       *colptr;
    const spm_complex64_t *values;

    const spm_complex64_t *x;
    spm_int_t              incx;

    spm_complex64_t       *y;
    spm_int_t              incy;

    __conj_fct_t           conj_fct;
    __loop_fct_t           loop_fct;
};

static inline int
__spm_zmatvec_sy_csr( const __spm_zmatvec_t *args )
{
    spm_int_t              baseval  = args->baseval;
    spm_int_t              n        = args->n;
    spm_complex64_t        alpha    = args->alpha;
    const spm_int_t       *rowptr   = args->rowptr;
    const spm_int_t       *colptr   = args->colptr;
    const spm_complex64_t *values   = args->values;
    const spm_complex64_t *x        = args->x;
    spm_int_t              incx     = args->incx;
    spm_complex64_t       *y        = args->y;
    spm_int_t              incy     = args->incy;
    const __conj_fct_t     conj_fct = args->conj_fct;
    spm_int_t              col, row, i;

    for( col=0; col<n; col++, colptr++ )
    {
        for( i=colptr[0]; i<colptr[1]; i++, rowptr++, values++ )
        {
            row = *rowptr - baseval;

            if ( row != col ) {
                y[ row * incy ] += alpha * conj_fct( *values ) * x[ col * incx ];
                y[ col * incy ] += alpha *         ( *values ) * x[ row * incx ];
            }
            else {
                y[ col * incy ] += alpha *         ( *values ) * x[ row * incx ];
            }
        }
    }
    return SPM_SUCCESS;
}

static inline int
__spm_zmatvec_sy_csc( const __spm_zmatvec_t *args )
{
    spm_int_t              baseval  = args->baseval;
    spm_int_t              n        = args->n;
    spm_complex64_t        alpha    = args->alpha;
    const spm_int_t       *rowptr   = args->rowptr;
    const spm_int_t       *colptr   = args->colptr;
    const spm_complex64_t *values   = args->values;
    const spm_complex64_t *x        = args->x;
    spm_int_t              incx     = args->incx;
    spm_complex64_t       *y        = args->y;
    spm_int_t              incy     = args->incy;
    const __conj_fct_t     conj_fct = args->conj_fct;
    spm_int_t              col, row, i;

    for( col=0; col<n; col++, colptr++ )
    {
        for( i=colptr[0]; i<colptr[1]; i++, rowptr++, values++ )
        {
            row = *rowptr - baseval;

            if ( row != col ) {
                y[ row * incy ] += alpha *         ( *values ) * x[ col * incx ];
                y[ col * incy ] += alpha * conj_fct( *values ) * x[ row * incx ];
            }
            else {
                y[ col * incy ] += alpha *         ( *values ) * x[ row * incx ];
            }
        }
    }
    return SPM_SUCCESS;
}

static inline int
__spm_zmatvec_ge_csc( const __spm_zmatvec_t *args )
{
    spm_int_t              baseval  = args->baseval;
    spm_int_t              n        = args->n;
    spm_complex64_t        alpha    = args->alpha;
    const spm_int_t       *rowptr   = args->rowptr;
    const spm_int_t       *colptr   = args->colptr;
    const spm_complex64_t *values   = args->values;
    const spm_complex64_t *x        = args->x;
    spm_int_t              incx     = args->incx;
    spm_complex64_t       *y        = args->y;
    spm_int_t              incy     = args->incy;
    const __conj_fct_t     conj_fct = args->conj_fct;
    spm_int_t              col, row, i;

    if ( args->follow_x ) {
        for( col=0; col<n; col++, colptr++, x+=incx )
        {
            for( i=colptr[0]; i<colptr[1]; i++, rowptr++, values++ )
            {
                row = *rowptr - baseval;
                y[ row * incy ] += alpha * conj_fct( *values ) * (*x);
            }
        }
    }
    else {
        for( col=0; col<n; col++, colptr++, y+=incy )
        {
            for( i=colptr[0]; i<colptr[1]; i++, rowptr++, values++ )
            {
                row = *rowptr - baseval;
                *y += alpha * conj_fct( *values ) * x[ row * incx ];
            }
        }
    }
    return SPM_SUCCESS;
}

static inline int
__spm_zmatvec_sy_ijv( const __spm_zmatvec_t *args )
{
    spm_int_t              baseval  = args->baseval;
    spm_int_t              nnz      = args->nnz;
    spm_complex64_t        alpha    = args->alpha;
    const spm_int_t       *rowptr   = args->rowptr;
    const spm_int_t       *colptr   = args->colptr;
    const spm_complex64_t *values   = args->values;
    const spm_complex64_t *x        = args->x;
    spm_int_t              incx     = args->incx;
    spm_complex64_t       *y        = args->y;
    spm_int_t              incy     = args->incy;
    const __conj_fct_t     conj_fct = args->conj_fct;
    spm_int_t              col, row, i;

    for( i=0; i<nnz; i++, colptr++, rowptr++, values++ )
    {
        row = *rowptr - baseval;
        col = *colptr - baseval;

        if ( row != col ) {
            y[ row * incy ] += alpha *         ( *values ) * x[ col * incx ];
            y[ col * incy ] += alpha * conj_fct( *values ) * x[ row * incx ];
        }
        else {
            y[ row * incy ] += alpha * ( *values ) * x[ col * incx ];
        }
    }
    return SPM_SUCCESS;
}

static inline int
__spm_zmatvec_ge_ijv( const __spm_zmatvec_t *args )
{
    spm_int_t              baseval  = args->baseval;
    spm_int_t              nnz      = args->nnz;
    spm_complex64_t        alpha    = args->alpha;
    const spm_int_t       *rowptr   = args->rowptr;
    const spm_int_t       *colptr   = args->colptr;
    const spm_complex64_t *values   = args->values;
    const spm_complex64_t *x        = args->x;
    spm_int_t              incx     = args->incx;
    spm_complex64_t       *y        = args->y;
    spm_int_t              incy     = args->incy;
    const __conj_fct_t     conj_fct = args->conj_fct;
    spm_int_t              col, row, i;

    for( i=0; i<nnz; i++, colptr++, rowptr++, values++ )
    {
        row = *rowptr - baseval;
        col = *colptr - baseval;

        y[ row * incy ] += alpha * conj_fct( *values ) * x[ col * incx ];
    }
    return SPM_SUCCESS;
}

#if !defined(LAPACKE_WITH_LASCL)
static inline void
__spm_zlascl( spm_complex64_t  alpha,
              spm_int_t        m,
              spm_int_t        n,
              spm_complex64_t *A,
              spm_int_t        lda )
{
    spm_int_t i, j;

    for( j=0; j<n; j++ ) {
        for( i=0; i<m; i++, A++ ) {
            *A *= alpha;
        }
        A += m - lda;
    }
}

#define LAPACKE_zlascl_work( _dir_, _uplo_, _kl_, _ku_, _cfrom_, _cto_, _m_, _n_, _A_, _lda_ ) \
    __spm_zlascl( (_cto_), (_m_), (_n_), (_A_), (_lda_) )

#endif

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_matvec
 *
 * @brief Compute a matrix-matrix product.
 *
 *    C = alpha * op(A) * op(B) + beta * C
 * or C = alpha * op(B) * op(A) + beta * C
 *
 * where A is a sparse matrix, B and C two dense matrices. And op(A), op(B) are one of:
 *
 *    op( A ) = A  or op( A ) = A' or op( A ) = conjg( A' )
 *
 *  alpha and beta are scalars.
 *
 *******************************************************************************
 *
 * @param[in] side
 *          Specifies whether A * B is computed or B * A
 *          - SpmLeft:  C = alpha * op(A) * op(B) + beta * C
 *          - SpmRight: C = alpha * op(B) * op(A) + beta * C
 *
 * @param[in] transA
 *          Specifies whether the sparse matrix A is not transposed, transposed
 *          or conjugate transposed:
 *          - SpmNoTrans
 *          - SpmTrans
 *          - SpmConjTrans
 *
 * @param[in] transB
 *          Specifies whether the dense matrix B is not transposed, transposed
 *          or conjugate transposed:
 *          - SpmNoTrans
 *          - SpmTrans
 *          - SpmConjTrans
 *
 * @param[in] K
 *          If side == SpmLeft, specifies the number of columns of the matrices op(B) and C.
 *          If side == SpmRight, specifies the number of rows of the matrices op(B) and C.
 *
 * @param[in] alpha
 *          alpha specifies the scalar alpha.
 *
 * @param[in] A
 *          The sparse matrix A
 *
 * @param[in] B
 *          The matrix B of size: ldb-by-Bn, with Bn = (K, A->m or A->n) based
 *          on the configuration of side, transA and transB.
 *
 * @param[in] ldb
 *          The leading dimension of the matrix B. ldb >= (A->m, A->n or K) based
 *          on the configuration of side, transA, and transB
 *
 * @param[in] beta
 *          beta specifies the scalar beta.
 *
 * @param[inout] C
 *          The matrix C of size ldc-by-Cn with Bn = (K, A->m or A->n) based
 *          on the configuration of side, transA and transB.
 *
 * @param[in] ldc
 *          The leading dimension of the matrix C. ldc >= (A->m, A->n or K) based
 *          on the configuration of side, transA, and transB
 *
 *   side  |  transA     | transB      |     B     |     C     |
 *   -----------------------------------------------------------
 *   Left  | NoTrans     | NoTrans     | A->n by K | A->m by K |
 *   Left  | NoTrans     | [Conj]Trans | K by A->n | A->m by K |
 *   Left  | [Conj]Trans | NoTrans     | A->m by K | A->n by K |
 *   Left  | [Conj]Trans | [Conj]Trans | K by A->m | A->n by K |
 *   Right | NoTrans     | NoTrans     | K by A->m | K by A->n |
 *   Right | NoTrans     | [Conj]Trans | A->m by K | K by A->n |
 *   Right | [Conj]Trans | NoTrans     | K by A->n | K by A->m |
 *   Right | [Conj]Trans | [Conj]Trans | A->n by K | K by A->m |
 *
 *******************************************************************************
 *
 * @retval SPM_SUCCESS if the y vector has been computed successfully,
 * @retval SPM_ERR_BADPARAMETER otherwise.
 *
 *******************************************************************************/
int
spm_zspmm( spm_side_t             side,
           spm_trans_t            transA,
           spm_trans_t            transB,
           spm_int_t              K,
           spm_complex64_t        alpha,
           const spmatrix_t      *A,
           const spm_complex64_t *B,
           spm_int_t              ldb,
           spm_complex64_t        beta,
           spm_complex64_t       *C,
           spm_int_t              ldc )
{
    int rc = SPM_SUCCESS;
    spm_int_t M, N, incx, incy, ldx, ldy, r;
    __spm_zmatvec_t args;

    if ( transB != SpmNoTrans ) {
        fprintf(stderr, "transB != SpmNoTrans not supported yet in spmv computations\n");
        assert( transB == SpmNoTrans );
        return SPM_ERR_BADPARAMETER;
    }

    if ( side == SpmLeft ) {
        M = A->n;
        N = K;

        incx = 1;
        incy = 1;
        ldx  = ldb;
        ldy  = ldc;
    }
    else {
        M = K;
        N = A->n;

        incx = ldb;
        incy = ldc;
        ldx  = 1;
        ldy  = 1;
    }

    if ( beta == 0. ) {
        LAPACKE_zlaset_work( LAPACK_COL_MAJOR, 'A', M, N, 0., 0., C, ldc );
    }
    else {
        LAPACKE_zlascl_work( LAPACK_COL_MAJOR, 'G', -1, -1, 1., beta, M, N, C, ldc );
    }

    if ( alpha == 0. ) {
        return SPM_SUCCESS;
    }

    {
        args.follow_x = 0;
        args.baseval  = spmFindBase( A );
        args.n        = A->n;
        args.nnz      = A->nnz;
        args.alpha    = alpha;
        args.rowptr   = A->rowptr;
        args.colptr   = A->colptr;
        args.values   = A->values;
        args.x        = B;
        args.incx     = incx;
        args.y        = C;
        args.incy     = incy;
        args.conj_fct = __fct_id;
        args.loop_fct = NULL;
    }

#if defined(PRECISION_c) || defined(PRECISION_z)
    if ( ( transA == SpmConjTrans ) ||
         ( A->mtxtype == SpmHermitian ) )
    {
        args.conj_fct = __fct_conj;
    }
#endif

    switch( A->fmttype ) {
    case SpmCSC:
    {
        /* Switch pointers and side to get the correct behaviour */
        if ( ((side == SpmLeft)  && (transA == SpmNoTrans)) ||
             ((side == SpmRight) && (transA != SpmNoTrans)) )
        {
            args.follow_x = 1;
        }
        else {
            args.follow_x = 0;
        }
        args.loop_fct = (A->mtxtype == SpmGeneral) ? __spm_zmatvec_ge_csc : __spm_zmatvec_sy_csc;
    }
    break;
    case SpmCSR:
    {
        /* Switch pointers and side to get the correct behaviour */
        if ( ((side == SpmLeft)  && (transA != SpmNoTrans)) ||
             ((side == SpmRight) && (transA == SpmNoTrans)) )
        {
            args.follow_x = 1;
        }
        else {
            args.follow_x = 0;
        }
        args.colptr = A->rowptr;
        args.rowptr = A->colptr;
        args.loop_fct = (A->mtxtype == SpmGeneral) ? __spm_zmatvec_ge_csc : __spm_zmatvec_sy_csr;
    }
    break;
    case SpmIJV:
    {
        if ( ((side == SpmLeft)  && (transA != SpmNoTrans)) ||
             ((side == SpmRight) && (transA == SpmNoTrans)) )
        {
            args.colptr = A->rowptr;
            args.rowptr = A->colptr;
        }
        args.loop_fct = (A->mtxtype == SpmGeneral) ? __spm_zmatvec_ge_ijv : __spm_zmatvec_sy_ijv;
    }
    break;
    default:
        return SPM_ERR_BADPARAMETER;
    }

    for( r=0; (r < N) && (rc == SPM_SUCCESS); r++ ) {
        args.x = B + r * ldx;
        args.y = C + r * ldy;
        rc = args.loop_fct( &args );
    }

    return rc;
}

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_matvec
 *
 * @brief compute the matrix-vector product:
 *          y = alpha * A + beta * y
 *
 * A is a SpmHermitian spm, alpha and beta are scalars, and x and y are
 * vectors, and A a symm.
 *
 *******************************************************************************
 *
 * @param[in] trans
 *          TODO
 *
 * @param[in] alphaptr
 *          alpha specifies the scalar alpha
 *
 * @param[in] spm
 *          The SpmHermitian spm.
 *
 * @param[in] xptr
 *          The vector x.
 *
 * @param[in] betaptr
 *          beta specifies the scalar beta
 *
 * @param[inout] yptr
 *          The vector y.
 *
 *******************************************************************************
 *
 * @retval SPM_SUCCESS if the y vector has been computed succesfully,
 * @retval SPM_ERR_BADPARAMETER otherwise.
 *
 *******************************************************************************/
int
spm_zspmv( spm_trans_t            trans,
           spm_complex64_t        alpha,
           const spmatrix_t      *A,
           const spm_complex64_t *x,
           spm_int_t              incx,
           spm_complex64_t        beta,
           spm_complex64_t       *y,
           spm_int_t              incy )
{
    int rc = SPM_SUCCESS;
    __spm_zmatvec_t args;

    if ( beta == 0. ) {
        memset( y, 0, A->n * sizeof(spm_complex64_t) );
    }
    else {
        cblas_zscal( A->n, CBLAS_SADDR(beta), y, incy );
    }

    if ( alpha == 0. ) {
        return SPM_SUCCESS;
    }

    {
        args.follow_x = 0;
        args.baseval  = spmFindBase( A );
        args.n        = A->n;
        args.nnz      = A->nnz;
        args.alpha    = alpha;
        args.rowptr   = A->rowptr;
        args.colptr   = A->colptr;
        args.values   = A->values;
        args.x        = x;
        args.incx     = incx;
        args.y        = y;
        args.incy     = incy;
        args.conj_fct = __fct_id;
        args.loop_fct = NULL;
    }

#if defined(PRECISION_c) || defined(PRECISION_z)
    if ( ( trans == SpmConjTrans ) ||
         ( A->mtxtype == SpmHermitian ) )
    {
        args.conj_fct = __fct_conj;
    }
#endif

    switch( A->fmttype ) {
    case SpmCSC:
    {
        args.follow_x = (trans == SpmNoTrans) ? 1 : 0;
        args.loop_fct = (A->mtxtype == SpmGeneral) ? __spm_zmatvec_ge_csc : __spm_zmatvec_sy_csc;
    }
    break;
    case SpmCSR:
    {
        /* Switch pointers and side to get the correct behaviour */
        args.follow_x = (trans == SpmNoTrans) ? 0 : 1;
        args.colptr = A->rowptr;
        args.rowptr = A->colptr;
        args.loop_fct = (A->mtxtype == SpmGeneral) ? __spm_zmatvec_ge_csc : __spm_zmatvec_sy_csr;
    }
    break;
    case SpmIJV:
    {
        if ( trans != SpmNoTrans ) {
            args.colptr = A->rowptr;
            args.rowptr = A->colptr;
        }
        args.loop_fct = (A->mtxtype == SpmGeneral) ? __spm_zmatvec_ge_ijv : __spm_zmatvec_sy_ijv;
    }
    break;
    default:
        return SPM_ERR_BADPARAMETER;
    }

    rc = args.loop_fct( &args );

    return rc;
}
