/**
 *
 * @file z_bvec_tests.c
 *
 * Tests for performance of gemv in parallel and sequential versions.
 * Times is computed on the average of 50 calls to gemv function.
 * Size of the matrix is given by the size of a Laplacian.
 *
 * @copyright 2015-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.1
 * @author Vincent Bridonneau
 * @date 2018-07-16
 *
 * @precisions normal z -> c d s
 *
 **/
#include <pastix.h>
#include <stdlib.h>
#include <string.h>
#include "refinement/z_refine_functions.h"
#include "common/common.h"
#include "common/flops.h"
#include "bcsc/bcsc.h"
#include "bcsc/bcsc_z.h"
#include "blend/solver.h"
#include "kernels/pastix_zcores.h"
#include <time.h>
#include <pastix/order.h>
#include <cblas.h>
#include <lapacke.h>
#if defined(HAVE_GETOPT_H)
#include <getopt.h>
#endif  /* defined(HAVE_GETOPT_H) */

struct z_argument_init_s
{
    pastix_int_t        n;
    pastix_complex64_t *x;
};

static inline void
z_init_smp( isched_thread_t *ctx,
            void            *args )
{
    struct z_argument_init_s *arg = (struct z_argument_init_s*)args;
    pastix_int_t              n   = arg->n;
    pastix_complex64_t       *x   = arg->x;
    pastix_int_t              begin, end, rank, size;

    size  = ctx->global_ctx->world_size;
    rank  = ctx->rank;
    begin = (n / size) * rank;
    if ( rank == (size - 1) ) {
        end = n;
    }
    else {
        end = (n / size) * (rank + 1);
    }

    core_zplrnt( (end-begin), 1, x, n, n, begin, 0, 7213 );
}

static inline void
z_init( pastix_data_t      *pastix_data,
        pastix_int_t        n,
        pastix_complex64_t *x )
{
    pastix_scheduler_t sched = pastix_data->iparm[IPARM_SCHEDULER];

    if ( sched == PastixSchedSequential ) {
        core_zplrnt( n, 1, x, n, n, 0, 0, 7213 );
    }
    else if ( sched >= 1 ) {
        struct z_argument_init_s args = { n, x };
        isched_parallel_call ( pastix_data->isched, z_init_smp, &args );
    }
}

int
z_bvec_gemv_check( int check, int m, int n,
                         pastix_int_t    *iparm,
                         pastix_fixdbl_t *dparm )
{
    pastix_data_t      *pastix_data = NULL; /*< Pointer to the storage structure required by pastix */
    int                 i;
    int                 rc = 0;
    pastix_complex64_t *A, *x, *y, *y0;
    pastix_complex64_t  alpha = 3.5;
    pastix_complex64_t  beta = -4.8;
    struct z_solver     solver;

    /**
     * Startup pastix to start the scheduler
     */
    pastixInit( &pastix_data, MPI_COMM_WORLD, iparm, dparm );

    z_refine_init( &solver, pastix_data );

    A = malloc( sizeof(pastix_complex64_t) * m * n );
    x = malloc( sizeof(pastix_complex64_t) * n );
    y = malloc( sizeof(pastix_complex64_t) * m );

    /**
     * generate matrice of size 'n * nrhs'
     */
    z_init( pastix_data, m * n, A );

    /**
     * generate vectors of size 'nrhs'
     */
    z_init( pastix_data, m, y );
    z_init( pastix_data, n, x );

    if ( check ) {
        double normA, normX, normY, normR, result;
        double eps = LAPACKE_dlamch_work( 'e' );

        normA = LAPACKE_zlange( LAPACK_COL_MAJOR, 'I', m, n, A, m );
        normX = LAPACKE_zlange( LAPACK_COL_MAJOR, '1', n, 1, x, n );
        normY = LAPACKE_zlange( LAPACK_COL_MAJOR, 'I', m, 1, y, m );

        y0 = malloc( sizeof(pastix_complex64_t) * m );
        solver.copy( pastix_data, m, y, y0 );

        solver.gemv( pastix_data, m, n, 3.5, A, m, x, -4.8, y );

        cblas_zgemv( CblasColMajor, CblasNoTrans, m, n,
                     CBLAS_SADDR(alpha), A, m, x, 1,
                     CBLAS_SADDR(beta), y0, 1 );

        core_zgeadd( PastixNoTrans, m, 1, -1., y, m, 1., y0, m );

        normR = LAPACKE_zlange( LAPACK_COL_MAJOR, 'I', m, 1, y0, m );

        result = normR / ((normA + normX + normY) * n * eps);
        if (  isnan(result) || isinf(result) || (result > 10.0) ) {
            rc = 1;
        }

        free( y0 );
    }
    else {
        double t, flops;
        Clock timer;

        clockStart(timer);
        for ( i = 0; i < 50 ; ++i) {
            solver.gemv( pastix_data, m, n, alpha, A, m, x, beta, y );
        }
        clockStop(timer);

        t = clockVal(timer) / 50.;
        flops = FLOPS_ZGEMV( m, n ) / t;
        printf("    Time for zgemv ( %ld x %ld ) : %e s ( %8.2g %cFlop/s)\n",
               (long)m, (long)n, t,
               pastix_print_value( flops ), pastix_print_unit( flops ) );
    }
    free( A );
    free( x );
    free( y );

    pastixFinalize( &pastix_data );

    return rc;
}

int
z_bvec_check( pastix_data_t *pastix_data,
              pastix_int_t m )
{
    Clock timer;
    int                 i;
    int                 rc = 0;
    pastix_complex64_t *x, *y;
    pastix_complex64_t  alpha = 3.5;
    struct z_solver     solver;

    z_refine_init( &solver, pastix_data );

    x = malloc( sizeof(pastix_complex64_t) * m );
    y = malloc( sizeof(pastix_complex64_t) * m );

    /**
     * generate vectors of size 'nrhs'
     */
    z_init( pastix_data, m, y );
    z_init( pastix_data, m, x );

    printf("========== copy function time ==========\n");
    clockStart(timer);
    for( i=0; i<50; ++i ) {
        solver.copy( pastix_data, m, x, y );
    }
    clockStop(timer);
    printf("    Time for copy (N= %ld)         : %e \n", (long)m, clockVal(timer)/50.);

    printf("========== axpy function time ==========\n");
    clockStart(timer);
    for( i=0; i<50; ++i ) {
        solver.axpy( pastix_data, m, alpha, x, y );
    }
    clockStop(timer);
    printf("    Time for axpy (N= %ld)         : %e \n", (long)m, clockVal(timer)/50.);

    printf("========== dot function time ==========\n");
    clockStart(timer);
    for( i=0; i<50; ++i ) {
        solver.dot( pastix_data, m, x, y );
    }
    clockStop(timer);
    printf("    Time for dot  (N= %ld)         : %e \n", (long)m, clockVal(timer)/50.);

    printf("========== norm function time ==========\n");
    clockStart(timer);
    for( i=0; i<50; ++i ) {
        solver.norm( pastix_data, m, x );
    }
    clockStop(timer);
    printf("    Time for norm (N= %ld)         : %e \n", (long)m, clockVal(timer)/50.);

    printf("========== scal function time ==========\n");
    clockStart(timer);
    for( i=0; i<50; ++i ) {
        solver.scal( pastix_data, m, alpha, x );
    }
    clockStop(timer);
    printf("    Time for scal (N= %ld)         : %e \n", (long)m, clockVal(timer)/50.);

    free( x );
    free( y );

    return rc;
}

int
z_bcsc_spmv_time( pastix_data_t *pastix_data,
                  spmatrix_t    *spm,
                  pastix_int_t   nrhs )
{
    int                 i;
    int                 rc = 0;
    pastix_complex64_t *x, *y;
    pastix_complex64_t  alpha = 3.5;
    pastix_complex64_t  beta = -4.8;
    struct z_solver     solver;
    double t, flops;
    Clock timer;

    z_refine_init( &solver, pastix_data );

    x = malloc( sizeof(pastix_complex64_t) * spm->nexp * nrhs );
    y = malloc( sizeof(pastix_complex64_t) * spm->nexp * nrhs );

    /**
     * generate matrices
     */
    z_init( pastix_data, spm->nexp * nrhs, x );
    z_init( pastix_data, spm->nexp * nrhs, y );

    clockStart(timer);
    for ( i = 0; i < 50 ; ++i) {
        solver.spmv( pastix_data, PastixNoTrans, alpha, x, beta, y );
    }
    clockStop(timer);

    t = clockVal(timer) / 50.;
    flops = FLOPS_ZGEMM( spm->nnzexp, 1, 1 ) / t;
    printf("    Time for zspmv ( n=%ld; nnz=%ld ) : %e s ( %8.2g %cFlop/s)\n",
           (long)spm->nexp, (long)spm->nnzexp, t,
           pastix_print_value( flops ), pastix_print_unit( flops ) );

    free( x );
    free( y );

    return rc;
}
