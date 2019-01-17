/**
 *
 * @file z_refine_grad.c
 *
 * PaStiX refinement functions implementations.
 *
 * @copyright 2015-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.1
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @author Theophile Terraz
 * @date 2018-07-16
 * @precisions normal z -> c d s
 *
 **/
#include "common.h"
#include "bcsc.h"
#include "z_refine_functions.h"

/**
 *******************************************************************************
 *
 * @ingroup pastix_refine
 *
 * z_grad_smp - Refine the solution using conjugate gradian method.
 *
 *******************************************************************************
 *
 * @param[in] pastix_data
 *          The PaStiX data structure that describes the solver instance.
 *
 * @param[out] x
 *          The solution vector.
 *
 * @param[in] b
 *          The right hand side member (only one).
 *
 *******************************************************************************
 *
 * @return Number of iterations
 *
 *******************************************************************************/
pastix_int_t z_grad_smp(pastix_data_t *pastix_data, void *x, void *b)
{
    struct z_solver     solver;
    pastix_int_t        n;
    Clock               refine_clk;
    pastix_fixdbl_t     t0        = 0;
    pastix_fixdbl_t     t3        = 0;
    int                 itermax;
    int                 nb_iter   = 0;
    int                 precond   = 1;
    pastix_complex64_t *gradr;
    pastix_complex64_t *gradp;
    pastix_complex64_t *gradz;
    pastix_complex64_t *grad2;
    double normb, normx, normr, alpha, beta;
    double resid_b, eps;

    memset( &solver, 0, sizeof(struct z_solver) );
    z_refine_init( &solver, pastix_data );

    if ( !(pastix_data->steps & STEP_NUMFACT) ) {
        precond = 0;
    }

    n       = pastix_data->bcsc->n;
    itermax = pastix_data->iparm[IPARM_ITERMAX];
    eps     = pastix_data->dparm[DPARM_EPSILON_REFINEMENT];

    /* Initialize vectors */
    gradr = (pastix_complex64_t *)solver.malloc(n * sizeof(pastix_complex64_t));
    gradp = (pastix_complex64_t *)solver.malloc(n * sizeof(pastix_complex64_t));
    gradz = (pastix_complex64_t *)solver.malloc(n * sizeof(pastix_complex64_t));
    grad2 = (pastix_complex64_t *)solver.malloc(n * sizeof(pastix_complex64_t));

    clockInit(refine_clk);
    clockStart(refine_clk);

    normb = solver.norm( pastix_data, n, b );
    normx = solver.norm( pastix_data, n, x );

    /* Compute r0 = b - A * x */
    solver.copy( pastix_data, n, b, gradr );
    if ( normx > 0. ) {
        solver.spmv( pastix_data, PastixNoTrans, -1., x, 1., gradr );
    }
    normr = solver.norm( pastix_data, n, gradr );
    resid_b = normr / normb;

    /* z = M^{-1} r */
    solver.copy( pastix_data, n, gradr, gradz );
    if ( precond ) {
        solver.spsv( pastix_data, gradz );
    }

    /* p = z */
    solver.copy( pastix_data, n, gradz, gradp );

    while ((resid_b > eps) && (nb_iter < itermax))
    {
        clockStop((refine_clk));
        t0 = clockGet();
        nb_iter++;

        /* grad2 = A * p */
        solver.spmv( pastix_data, PastixNoTrans, 1.0, gradp, 0., grad2 );

        /* alpha = <r, z> / <Ap, p> */
        beta  = solver.dot( pastix_data, n, gradr, gradz );
        alpha = solver.dot( pastix_data, n, grad2, gradp );
        alpha = beta / alpha;

        /* x = x + alpha * p */
        solver.axpy( pastix_data, n, alpha, gradp, x );

        /* r = r - alpha * A * p */
        solver.axpy( pastix_data, n, -alpha, grad2, gradr );

        /* z = M-1 * r */
        solver.copy( pastix_data, n, gradr, gradz );
        if ( precond ) {
            solver.spsv( pastix_data, gradz );
        }

        /* beta = <r', z> / <r, z> */
        alpha = solver.dot( pastix_data, n, gradr, gradz );
        beta  = alpha / beta;

        /* p = z + beta * p */
        solver.scal( pastix_data, n, beta, gradp );
        solver.axpy( pastix_data, n, 1., gradz, gradp );

        normr = solver.norm( pastix_data, n, gradr );
        resid_b = normr / normb;

        clockStop((refine_clk));
        t3 = clockGet();
        if ( pastix_data->iparm[IPARM_VERBOSE] > PastixVerboseNot ) {
            solver.output_oneiter( t0, t3, resid_b, nb_iter );
        }
        t0 = t3;
    }

    solver.output_final(pastix_data, resid_b, nb_iter, t3, x, x);

    solver.free((void*) gradr);
    solver.free((void*) gradp);
    solver.free((void*) gradz);
    solver.free((void*) grad2);

    return nb_iter;
}
