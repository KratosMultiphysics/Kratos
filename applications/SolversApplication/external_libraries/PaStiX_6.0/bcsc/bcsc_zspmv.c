/**
 *
 * @file bcsc_zspmv.c
 *
 * Functions computing matrix-vector products for the BCSC
 *
 * @copyright 2004-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.1
 * @author Mathieu Faverge
 * @author Théophile Terraz
 * @author Vincent Bridonneau
 * @date 2018-07-16
 *
 * @precisions normal z -> c d s
 *
 **/
#include "common.h"
#include <math.h>
#include "bcsc.h"
#include "bcsc_z.h"
#include "solver.h"
#include "pastix/datatypes.h"

#ifndef DOXYGEN_SHOULD_SKIP_THIS

typedef void (*bcsc_zspmv_Ax_fct_t)( const pastix_bcsc_t *,
                                     const bcsc_cblk_t *,
                                     pastix_complex64_t,
                                     const pastix_complex64_t *,
                                     const pastix_complex64_t *,
                                     pastix_complex64_t,
                                     pastix_complex64_t * );
static inline void
__bcsc_zspmv_by( pastix_int_t        n,
                 pastix_complex64_t  beta,
                 pastix_complex64_t *y )
{
    if( beta != (pastix_complex64_t)0.0 )
    {
        pastix_int_t j;
        for( j=0; j<n; j++, y++ )
        {
            (*y) *= beta;
        }
    }
    else
    {
        memset( y, 0, n * sizeof(pastix_complex64_t) );
    }
}

static inline void
__bcsc_zspmv_Ax( const pastix_bcsc_t      *bcsc,
                 const bcsc_cblk_t        *cblk,
                 pastix_complex64_t        alpha,
                 const pastix_complex64_t *A,
                 const pastix_complex64_t *x,
                 pastix_complex64_t        beta,
                 pastix_complex64_t       *y )
{
    pastix_int_t i, j;

    __bcsc_zspmv_by( cblk->colnbr, beta, y );

    for( j=0; j<cblk->colnbr; j++, y++ )
    {
        for( i=cblk->coltab[j]; i< cblk->coltab[j+1]; i++ )
        {
            *y += alpha * A[i] * x[ bcsc->rowtab[i] ];
        }
    }
}

static inline void
__bcsc_zspmv_Ax_ind( const pastix_bcsc_t      *bcsc,
                     pastix_complex64_t        alpha,
                     const pastix_complex64_t *A,
                     const pastix_complex64_t *x,
                     pastix_complex64_t        beta,
                     pastix_complex64_t       *y )
{
    const pastix_complex64_t *xptr = x;
    pastix_int_t bloc, i, j;

    __bcsc_zspmv_by( bcsc->n, beta, y );

    for( bloc=0; bloc<bcsc->cscfnbr; bloc++ )
    {
        for( j=0; j < bcsc->cscftab[bloc].colnbr; j++, xptr++ )
        {
            for( i = bcsc->cscftab[bloc].coltab[j]; i < bcsc->cscftab[bloc].coltab[j+1]; i++ )
            {
                y[ bcsc->rowtab[i] ] += alpha * A[i] * (*xptr);
            }
        }
    }
}

#if defined(PRECISION_z) || defined(PRECISION_c)
static inline void
__bcsc_zspmv_conjAx( const pastix_bcsc_t      *bcsc,
                     const bcsc_cblk_t        *cblk,
                     pastix_complex64_t        alpha,
                     const pastix_complex64_t *A,
                     const pastix_complex64_t *x,
                     pastix_complex64_t        beta,
                     pastix_complex64_t       *y )
{
    pastix_int_t i, j;

    __bcsc_zspmv_by( cblk->colnbr, beta, y );

    for( j=0; j<cblk->colnbr; j++, y++ )
    {
        for( i=cblk->coltab[j]; i< cblk->coltab[j+1]; i++ )
        {
            *y += alpha * conj( A[i] ) * x[ bcsc->rowtab[i] ];
        }
    }
}
#endif

static inline void
__bcsc_zspmv_loop( pastix_trans_t            trans,
                   pastix_complex64_t        alpha,
                   const pastix_bcsc_t      *bcsc,
                   const pastix_complex64_t *x,
                   pastix_complex64_t        beta,
                   pastix_complex64_t       *y,
                   pastix_int_t              rank,
                   pastix_int_t              begin,
                   pastix_int_t              end )
{
    bcsc_zspmv_Ax_fct_t  zspmv_Ax = __bcsc_zspmv_Ax;
    pastix_complex64_t  *valptr = NULL;
    pastix_int_t         bloc;
    bcsc_cblk_t         *cblk;

    /*
     * There are three cases:
     *    We can use the Lvalues pointer directly:
     *          - The matrix is general and we use A^t
     *          - the matrix is symmetric or hermitian
     *    We can use the Uvalues pointer directly
     *          - The matrix is general and we use A
     *    We have to use Lvalues per row (instead of column)
     *          - The matrix A is general and Uvalues is unavailable
     *
     * To this, we have to add the conj call if ConjTrans or Hermitian
     *
     *     Mtxtype   | trans asked | algo applied
     *     ++++++++++++++++++++++++++++++++++++
     +     General   | NoTrans     | U if possible, otherwise indirect L
     +     General   | Trans       | L
     +     General   | ConjTrans   | conj(L)
     +     Symmetric | NoTrans     | L
     +     Symmetric | Trans       | L
     +     Symmetric | ConjTrans   | conj(L)
     +     Hermitian | NoTrans     | conj(L)
     +     Hermitian | Trans       | conj(L)
     +     Hermitian | ConjTrans   | L
     */
    cblk   = bcsc->cscftab + begin;
    valptr = (pastix_complex64_t*)bcsc->Lvalues;

    if ( (bcsc->mtxtype == PastixGeneral) && (trans == PastixNoTrans) )
    {
        /* U */
        if ( bcsc->Uvalues != NULL ) {
            valptr = (pastix_complex64_t*)bcsc->Uvalues;
        }
        /* Indirect L */
        else {
            /* Execute in sequential */
            if ( rank != 0 ) {
                return;
            }
            return __bcsc_zspmv_Ax_ind( bcsc, alpha, valptr, x, beta, y );
        }
    }
#if defined(PRECISION_z) || defined(PRECISION_c)
    /* Conj(L) */
    else if ( ( (bcsc->mtxtype == PastixGeneral  ) && (trans == PastixConjTrans) ) ||
              ( (bcsc->mtxtype == PastixSymmetric) && (trans == PastixConjTrans) ) ||
              ( (bcsc->mtxtype == PastixHermitian) && (trans != PastixConjTrans) ) )
    {
        zspmv_Ax = __bcsc_zspmv_conjAx;
    }
#endif /* defined(PRECISION_z) || defined(PRECISION_c) */

    for( bloc=begin; bloc<end; bloc++ )
    {
        zspmv_Ax( bcsc, cblk, alpha, valptr, x, beta, y );

        y += cblk->colnbr;
        cblk++;
    }
}

#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/**
 *******************************************************************************
 *
 * @ingroup bcsc
 *
 * @brief Compute the matrix-vector product  y = alpha * A * x + beta * y
 * (Sequential version)
 *
 * Where A is given in the bcsc format, x and y are two vectors of size n, and
 * alpha and beta are two scalars.
 *
 *******************************************************************************
 *
 * @param[in] pastix_data
 *          Provide information about bcsc
 *
 * @param[in] trans
 *          Specifies whether the matrix A from the bcsc is transposed, not
 *          transposed or conjugate transposed:
 *            = PastixNoTrans:   A is not transposed;
 *            = PastixTrans:     A is transposed;
 *            = PastixConjTrans: A is conjugate transposed.
 *
 * @param[in] alpha
 *          alpha specifies the scalar alpha
 *
 * @param[in] x
 *          The vector x.
 *
 * @param[in] beta
 *          beta specifies the scalar beta
 *
 * @param[inout] y
 *          The vector y.
 *
 *******************************************************************************/
void
bcsc_zspmv_seq( const pastix_data_t      *pastix_data,
                pastix_trans_t            trans,
                pastix_complex64_t        alpha,
                const pastix_complex64_t *x,
                pastix_complex64_t        beta,
                pastix_complex64_t       *y )
{
    pastix_bcsc_t *bcsc = pastix_data->bcsc;

    if( (bcsc == NULL) || (y == NULL) || (x == NULL) ) {
        return;
    }

    __bcsc_zspmv_loop( trans, alpha, bcsc, x, beta, y,
                       0, 0, bcsc->cscfnbr );
}

/**
 * @brief Data structure for parallel arguments of spmv functions
 */
struct z_argument_spmv_s
{
  pastix_trans_t            trans;
  pastix_complex64_t        alpha;
  const pastix_bcsc_t      *bcsc;
  const pastix_complex64_t *x;
  pastix_complex64_t        beta;
  pastix_complex64_t       *y;
  SolverMatrix             *mtx;
  pastix_int_t             *start_indexes; /* starting position for each thread*/
  pastix_int_t             *start_bloc;
};

/**
 *******************************************************************************
 *
 * @ingroup bcsc_internal
 *
 * @brief Compute the matrix-vector product  y = alpha * op(A) * x + beta * y
 *
 * Where A is given in the bcsc format, x and y are two vectors of size n, and
 * alpha and beta are two scalars.
 * The op function is specified by the trans parameter and performs the
 * operation as follows:
 *              trans = PastixNoTrans   y := alpha*A       *x + beta*y
 *              trans = PastixTrans     y := alpha*A'      *x + beta*y
 *              trans = PastixConjTrans y := alpha*conj(A')*x + beta*y
 *
 *******************************************************************************
 *
 * @param[in] ctx
 *          the context of the current thread
 *
 * @param[inout] args
 *          The parameter as specified in bcsc_zspmv.
 *
 *******************************************************************************/
static inline void
pthread_bcsc_zspmv( isched_thread_t *ctx,
                    void            *args )
{
    struct z_argument_spmv_s *arg    = (struct z_argument_spmv_s*)args;
    const pastix_bcsc_t      *bcsc   = arg->bcsc;
    pastix_int_t              begin, end, size, rank;
    pastix_int_t             *start_indexes = arg->start_indexes;
    pastix_int_t             *start_bloc    = arg->start_bloc;

    rank = (pastix_int_t)ctx->rank;
    size = (pastix_int_t)ctx->global_ctx->world_size;

    begin = start_bloc[rank];
    if ( rank == (size - 1) )
    {
        end = bcsc->cscfnbr;
    }
    else {
        end = start_bloc[rank + 1];
    }

    __bcsc_zspmv_loop( arg->trans, arg->alpha, bcsc, arg->x,
                       arg->beta, arg->y + start_indexes[rank],
                       rank, begin, end );
}

/**
 *******************************************************************************
 *
 * @ingroup bcsc_internal
 *
 * @brief Compute the matrix-vector product  y = alpha * op(A) * x + beta * y
 *
 * Where A is given in the bcsc format, x and y are two vectors of size n, and
 * alpha and beta are two scalars.
 * The op function is specified by the trans parameter and performs the
 * operation as follows:
 *              trans = PastixNoTrans   y := alpha*A       *x + beta*y
 *              trans = PastixTrans     y := alpha*A'      *x + beta*y
 *              trans = PastixConjTrans y := alpha*conj(A')*x + beta*y
 *
 *******************************************************************************
 *
 * @param[in] ctx
 *          the context of the current thread
 *
 * @param[inout] args
 *          The parameter as specified in bcsc_zspmv.
 *
 *******************************************************************************/
static inline void
pthread_bcsc_zspmv_tasktab( isched_thread_t *ctx,
                            void            *args )
{
    bcsc_zspmv_Ax_fct_t       zspmv_Ax = __bcsc_zspmv_Ax;
    struct z_argument_spmv_s *arg    = (struct z_argument_spmv_s*)args;
    pastix_trans_t            trans  = arg->trans;
    pastix_complex64_t        alpha  = arg->alpha;
    const pastix_bcsc_t      *bcsc   = arg->bcsc;
    const pastix_complex64_t *x      = arg->x;
    pastix_complex64_t        beta   = arg->beta;
    pastix_complex64_t       *y      = arg->y;
    pastix_complex64_t       *valptr = NULL;
    pastix_complex64_t       *yptr;
    pastix_int_t              rank;
    SolverMatrix             *mtx = arg->mtx;
    pastix_int_t              tasknbr, *tasktab;
    pastix_int_t              ii, task_id;
    SolverCblk               *solv_cblk;
    bcsc_cblk_t              *bcsc_cblk;
    Task                     *t;

    rank = (pastix_int_t)ctx->rank;

    tasknbr = mtx->ttsknbr[rank];
    tasktab = mtx->ttsktab[rank];

    /*
     * There are three cases:
     *    We can use the Lvalues pointer directly:
     *          - The matrix is general and we use A^t
     *          - The matrix is symmetric or hermitian
     *    We can use the Uvalues pointer directly
     *          - The matrix is general and we use A
     *    We have to use Lvalues per row (instead of column)
     *          - The matrix A is general and Uvalues is unavailable
     *
     * To this, we have to add the conj call if ConjTrans or Hermitian
     *
     *     Mtxtype   | trans asked | algo applied
     *     ++++++++++++++++++++++++++++++++++++
     +     General   | NoTrans     | U if possible, otherwise indirect L
     +     General   | Trans       | L
     +     General   | ConjTrans   | conj(L)
     +     Symmetric | NoTrans     | L
     +     Symmetric | Trans       | L
     +     Symmetric | ConjTrans   | conj(L)
     +     Hermitian | NoTrans     | conj(L)
     +     Hermitian | Trans       | conj(L)
     +     Hermitian | ConjTrans   | L
     */
    valptr = (pastix_complex64_t*)bcsc->Lvalues;

    if ( (bcsc->mtxtype == PastixGeneral) && (trans == PastixNoTrans) )
    {
        /* U */
        if ( bcsc->Uvalues != NULL ) {
            valptr = (pastix_complex64_t*)bcsc->Uvalues;
        }
        /* Indirect L */
        else {
            /* Execute in sequential */
            if ( rank != 0 ) {
                return;
            }
            return __bcsc_zspmv_Ax_ind( bcsc, alpha, valptr, x, beta, y );
        }
    }
#if defined(PRECISION_z) || defined(PRECISION_c)
    /* Conj(L) */
    else if ( ( (bcsc->mtxtype == PastixGeneral  ) && (trans == PastixConjTrans) ) ||
              ( (bcsc->mtxtype == PastixSymmetric) && (trans == PastixConjTrans) ) ||
              ( (bcsc->mtxtype == PastixHermitian) && (trans != PastixConjTrans) ) )
    {
        zspmv_Ax = __bcsc_zspmv_conjAx;
    }
#endif /* defined(PRECISION_z) || defined(PRECISION_c) */

    for (ii=0; ii<tasknbr; ii++)
    {
        task_id = tasktab[ii];
        t = mtx->tasktab + task_id;

        solv_cblk = mtx->cblktab + t->cblknum;
        bcsc_cblk = bcsc->cscftab + t->cblknum;
        yptr = y + solv_cblk->fcolnum;

        zspmv_Ax( bcsc, bcsc_cblk, alpha, valptr, x, beta, yptr );
    }
}

/**
 *******************************************************************************
 *
 * @ingroup bcsc_internal
 *
 * @brief Initialize indexes for vector pointer and bloc indexes
 *   for parallel version of spmv.
 *
 *   This function Initial indexes for each thread
 *   in order to computes it once instead of once per thread. This is a more
 *   sophisticated version trying to balance the load for each thread in terms
 *   of bloc size.
 *
 *******************************************************************************
 *
 * @param[in] pastix_data
 *          The pastix_data structure providing number of threads and holding
 *          the A matrix.
 *
 * @param[out] args
 *          The argument containing arrays to initialise (blocs and indexes).
 *
 *******************************************************************************/
static inline void
bcsc_zspmv_get_balanced_indexes( const pastix_data_t      *pastix_data,
                                 struct z_argument_spmv_s *args )
{
    pastix_int_t rank, bloc, size;
    pastix_int_t ratio, total, load;
    pastix_bcsc_t *bcsc = pastix_data->bcsc;
    bcsc_cblk_t *cblk = bcsc->cscftab;

    if ( bcsc->mtxtype != PastixGeneral ) {
        total = 2 * pastix_data->csc->nnzexp - bcsc->n;
    } else {
        total = pastix_data->csc->nnzexp;
    }
    size  = pastix_data->isched->world_size;
    ratio = pastix_iceil( total, size );
    load  = 0;

    args->start_bloc[0]    = 0;
    args->start_indexes[0] = 0;

    for ( bloc = 0, rank = 1; bloc < bcsc->cscfnbr; ++bloc, ++cblk )
    {
        if ( load >= ratio ) {
            assert( rank < size );

            args->start_bloc[rank]    = bloc;
            args->start_indexes[rank] = pastix_data->solvmatr->cblktab[bloc].fcolnum;

            rank ++;
            total -= load;
            load = 0;
        }
        load += cblk->coltab[cblk->colnbr] - cblk->coltab[0];
    }

    total -= load;
    assert( total == 0 );

    for ( ; rank < size; rank ++ ) {
        args->start_bloc[rank]    = bcsc->cscfnbr;
        args->start_indexes[rank] = bcsc->n;
    }
}

/**
 *******************************************************************************
 *
 * @ingroup bcsc
 *
 * @brief Perform y = alpha A x + beta y (Parallel version)
 *
 * This functions is parallelized through the internal static scheduler.
 *
 *******************************************************************************
 *
 * @param[in] pastix_data
 *          The pastix_data structure that holds the A matrix.
 *
 * @param[in] trans
 *          Specifies whether the matrix A from the bcsc is transposed, not
 *          transposed or conjugate transposed:
 *            = PastixNoTrans:   A is not transposed;
 *            = PastixTrans:     A is transposed;
 *            = PastixConjTrans: A is conjugate transposed.
 *
 * @param[in] alpha
 *          The scalar alpha.
 *
 * @param[in] x
 *          The vector x
 *
 * @param[in] beta
 *          The scalar beta.
 *
 * @param[inout] y
 *          On entry, the vector y
 *          On exit, alpha A x + y
 *
 *******************************************************************************/
void
bcsc_zspmv_smp( const pastix_data_t      *pastix_data,
                pastix_trans_t            trans,
                pastix_complex64_t        alpha,
                const pastix_complex64_t *x,
                pastix_complex64_t        beta,
                pastix_complex64_t       *y )
{
    pastix_bcsc_t *bcsc = pastix_data->bcsc;
    struct z_argument_spmv_s arg = { trans, alpha, bcsc, x, beta, y,
                                     pastix_data->solvmatr, NULL, NULL };

    if( (bcsc == NULL) || (y == NULL) || (x == NULL) ) {
        return;
    }

    isched_parallel_call( pastix_data->isched, pthread_bcsc_zspmv_tasktab, &arg );

#if 0
    /*
     * Version that balances the number of nnz per thread, instead of exploiting
     * the tasktab array.
     */
    {
        MALLOC_INTERN( arg.start_indexes, 2 * pastix_data->isched->world_size, pastix_int_t );
        arg.start_bloc = arg.start_indexes + pastix_data->isched->world_size;

        bcsc_zspmv_get_balanced_indexes( pastix_data, &arg );

        isched_parallel_call ( pastix_data->isched, pthread_bcsc_zspmv, &arg );

        memFree_null( arg.start_indexes );
    }
#endif
}

/**
 *******************************************************************************
 *
 * @brief Compute the matrix-vector product  y = alpha * op(A) * x + beta * y
 *
 * Where A is given in the bcsc format, x and y are two vectors of size n, and
 * alpha and beta are two scalars.
 * The op function is specified by the trans parameter and performs the
 * operation as follows:
 *              trans = PastixNoTrans   y := alpha*A       *x + beta*y
 *              trans = PastixTrans     y := alpha*A'      *x + beta*y
 *              trans = PastixConjTrans y := alpha*conj(A')*x + beta*y
 *
 * This function is used only in testings.
 *
 *******************************************************************************
 *
 * @param[in] pastix_data
 *          Provide information about bcsc, and select the scheduling version
 *          based on iparm[IPARM_SCHEDULER].
 *
 * @param[in] trans
 *          Specifies whether the matrix A from the bcsc is transposed, not
 *          transposed or conjugate transposed:
 *            = PastixNoTrans:   A is not transposed;
 *            = PastixTrans:     A is transposed;
 *            = PastixConjTrans: A is conjugate transposed.
 *
 * @param[in] alpha
 *          alpha specifies the scalar alpha
 *
 * @param[in] x
 *          The vector x.
 *
 * @param[in] beta
 *          beta specifies the scalar beta
 *
 * @param[inout] y
 *          The vector y.
 *
 *******************************************************************************/
void
bcsc_zspmv( const pastix_data_t      *pastix_data,
            pastix_trans_t            trans,
            pastix_complex64_t        alpha,
            const pastix_complex64_t *x,
            pastix_complex64_t        beta,
            pastix_complex64_t       *y )
{
    pastix_int_t *iparm = pastix_data->iparm;

    if ( iparm[IPARM_SCHEDULER] == PastixSchedStatic ) {
        bcsc_zspmv_smp( pastix_data, trans, alpha, x, beta, y );
    }
    else {
        bcsc_zspmv_seq( pastix_data, trans, alpha, x, beta, y );
    }
}
