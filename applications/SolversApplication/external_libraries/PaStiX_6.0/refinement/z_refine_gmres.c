/**
 *
 * @file z_refine_gmres.c
 *
 * PaStiX refinement functions implementations.
 *
 * @copyright 2015-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.1
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Theophile Terraz
 * @author Xavier Lacoste
 * @date 2018-07-16
 * @precisions normal z -> c d s
 *
 **/
#include "common.h"
#include "cblas.h"
#include "bcsc.h"
#include "z_refine_functions.h"

/**
 *******************************************************************************
 *
 * @ingroup pastix_refine
 *
 * z_gmres_smp - Function computing GMRES iterative refinement.
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
 *******************************************************************************
 *
 * @return Number of iterations
 *
 *******************************************************************************/
pastix_int_t z_gmres_smp(pastix_data_t *pastix_data, void *x, void *b)
{
    struct z_solver     solver;
    Clock               refine_clk;
    pastix_complex64_t *gmHi, *gmH;
    pastix_complex64_t *gmVi, *gmV;
    pastix_complex64_t *gmWi, *gmW;
    pastix_complex64_t *gmcos, *gmsin;
    pastix_complex64_t *gmG;
#if defined(PASTIX_DEBUG_GMRES)
    pastix_complex64_t *dbg_x, *dbg_r, *dbg_G;
#endif
    pastix_complex64_t  tmp;
    pastix_fixdbl_t     t0, t3;
    double              eps, resid, resid_b;
    double              norm, normb, normx;
    pastix_int_t        n, im, im1, itermax;
    pastix_int_t        i, j,  ldw, iters;
    int                 outflag, inflag;
    int                 savemem = 0;
    int                 precond = 1;

    memset( &solver, 0, sizeof(struct z_solver) );
    z_refine_init( &solver, pastix_data );

    /* if ( pastix_data->bcsc->mtxtype == PastixHermitian ) { */
    /*     /\* Check if we need dotu for non hermitian matrices (CEA patch) *\/ */
    /*     solver.dot = &z_Pastix_Dotc; */
    /* } */

    /* Get the parameters */
    n       = pastix_data->bcsc->n;
    im      = pastix_data->iparm[IPARM_GMRES_IM];
    im1     = im + 1;
    itermax = pastix_data->iparm[IPARM_ITERMAX];
    eps     = pastix_data->dparm[DPARM_EPSILON_REFINEMENT];
    ldw     = n;

    if ( !(pastix_data->steps & STEP_NUMFACT) ) {
        precond = 0;
    }

    if ((!precond) || savemem ) {
        ldw = 0;
    }

    gmcos = (pastix_complex64_t *)solver.malloc(im  * sizeof(pastix_complex64_t));
    gmsin = (pastix_complex64_t *)solver.malloc(im  * sizeof(pastix_complex64_t));
    gmG   = (pastix_complex64_t *)solver.malloc(im1 * sizeof(pastix_complex64_t));

    /**
     * H stores the h_{i,j} elements ot the upper hessenberg matrix H (See Alg. 9.5 p 270)
     * V stores the v_{i} vectors
     * W stores the M^{-1} v_{i} vectors to avoid the application of the
     *          preconditioner on the output result (See line 11 of Alg 9.5)
     *
     * If no preconditioner is applied, or the user wants to save memory, W
     * stores only temporarily one vector for the Ax product (ldw is set to 0 to
     * reuse the same vector at each iteration)
     */
    gmH = (pastix_complex64_t *)solver.malloc(im * im1 * sizeof(pastix_complex64_t));
    gmV = (pastix_complex64_t *)solver.malloc(n  * im1 * sizeof(pastix_complex64_t));
    if (precond && (!savemem) ) {
        gmW = (pastix_complex64_t *)solver.malloc(n * im  * sizeof(pastix_complex64_t));
    }
    else {
        gmW = (pastix_complex64_t *)solver.malloc(n       * sizeof(pastix_complex64_t));
    }
    memset( gmH, 0, im * im1 * sizeof(pastix_complex64_t) );

#if defined(PASTIX_DEBUG_GMRES)
    dbg_x = (pastix_complex64_t *)solver.malloc(n   * sizeof(pastix_complex64_t));
    dbg_r = (pastix_complex64_t *)solver.malloc(n   * sizeof(pastix_complex64_t));
    dbg_G = (pastix_complex64_t *)solver.malloc(im1 * sizeof(pastix_complex64_t));
    solver.copy( pastix_data, n, x, dbg_x );
#endif

    normb = solver.norm( pastix_data, n, b );
    normx = solver.norm( pastix_data, n, x );

    clockInit(refine_clk);
    clockStart(refine_clk);

    /**
     * Algorithm from Iterative Methods for Sparse Linear systems, Y. Saad, Second Ed. p267-273
     *
     * The version implemented is the Right preconditioned algorithm.
     */
    outflag = 1;
    iters = 0;
    while (outflag)
    {
        /* Initialize v_{0} and w_{0} */
        gmVi = gmV;

        /* Compute r0 = b - A * x */
        solver.copy( pastix_data, n, b, gmVi );
        if ( normx > 0. ) {
            solver.spmv( pastix_data, PastixNoTrans, -1., x, 1., gmVi );
        }

        /* Compute resid = ||r0||_f */
        resid = solver.norm( pastix_data, n, gmVi );
        resid_b = resid / normb;

        /* If residual is small enough, exit */
        if ( resid_b <= eps )
        {
            outflag = 0;
            break;
        }

        /* Compute v0 = r0 / resid */
        tmp = (pastix_complex64_t)( 1.0 / resid );
        solver.scal( pastix_data, n, tmp, gmVi );

        gmG[0] = (pastix_complex64_t)resid;
        inflag = 1;
        i = -1;
        gmHi = gmH - im1;
        gmWi = gmW - ldw;

        while( inflag )
        {
            clockStop( refine_clk );
            t0 = clockGet();

            i++;

            /* Set H and W pointers to the beginning of columns i */
            gmHi = gmHi + im1;
            gmWi = gmWi + ldw;

            /* Backup v_{i} into w_{i} for the end */
            solver.copy( pastix_data, n, gmVi, gmWi );

            /* Compute w_{i} = M^{-1} v_{i} */
            if ( precond ) {
                solver.spsv( pastix_data, gmWi );
            }

            /* v_{i+1} = A (M^{-1} v_{i}) = A w_{i} */
            gmVi += n;
            solver.spmv( pastix_data, PastixNoTrans, 1.0, gmWi, 0., gmVi );

            /* Classical Gram-Schmidt */
            for (j=0; j<=i; j++)
            {
                /* Compute h_{j,i} = < v_{i+1}, v_{j} > */
                gmHi[j] = solver.dot( pastix_data, n, gmVi, gmV + j * n );

                /* Compute v_{i+1} = v_{i+1} - h_{j,i} v_{j} */
                solver.axpy( pastix_data, n, -1. * gmHi[j],  gmV + j * n, gmVi );
            }

            /* Compute || v_{i+1} ||_f */
            norm = solver.norm( pastix_data, n, gmVi );
            gmHi[i+1] = norm;

            /* Compute v_{i+1} = v_{i+1} / h_{i+1,i} iff h_{i+1,i} is not too small */
            if ( norm > 1e-50 )
            {
                tmp = (pastix_complex64_t)(1.0 / norm);
                solver.scal( pastix_data, n, tmp, gmVi );
            }

            /* Apply the previous Givens rotation to the new column (should call LAPACKE_zrot_work())*/
            for (j=0; j<i;j++)
            {
                /*
                 * h_{j,  i} = cos_j * h_{j,  i} +      sin_{j}  * h_{j+1, i}
                 * h_{j+1,i} = cos_j * h_{j+1,i} - conj(sin_{j}) * h_{j,   i}
                 */
                tmp = gmHi[j];
                gmHi[j]   = gmcos[j] * tmp       +      gmsin[j]  * gmHi[j+1];
                gmHi[j+1] = gmcos[j] * gmHi[j+1] - conj(gmsin[j]) * tmp;
            }

            /*
             * Compute the new Givens rotation (zrotg)
             *
             * t = sqrt( h_{i,i}^2 + h_{i+1,i}^2 )
             * cos = h_{i,i}   / t
             * sin = h_{i+1,i} / t
             */
            {
                tmp = csqrt( gmHi[i]   * gmHi[i] +
                             gmHi[i+1] * gmHi[i+1] );

                if ( cabs(tmp) <= eps ) {
                    tmp = (pastix_complex64_t)eps;
                }
                gmcos[i] = gmHi[i]   / tmp;
                gmsin[i] = gmHi[i+1] / tmp;
            }

            /* Update the residuals (See p. 168, eq 6.35) */
            gmG[i+1] = -gmsin[i] * gmG[i];
            gmG[i]   =  gmcos[i] * gmG[i];

            /* Apply the last Givens rotation */
            gmHi[i] = gmcos[i] * gmHi[i] + gmsin[i] * gmHi[i+1];

            /* (See p. 169, eq 6.42) */
            resid = cabs( gmG[i+1] );

            resid_b = resid / normb;
            iters++;
            if ( (i+1 >= im) ||
                 (resid_b <= eps) ||
                 (iters >= itermax) )
            {
                inflag = 0;
            }

            clockStop((refine_clk));
            t3 = clockGet();
            if ( pastix_data->iparm[IPARM_VERBOSE] > PastixVerboseNot ) {
                solver.output_oneiter( t0, t3, resid_b, iters );

#if defined(PASTIX_DEBUG_GMRES)
                {
                    double normr2;

                    /* Compute y_m = H_m^{-1} g_m (See p. 169) */
                    memcpy( dbg_G, gmG, im1 * sizeof(pastix_complex64_t) );
                    cblas_ztrsv( CblasColMajor, CblasUpper, CblasNoTrans, CblasNonUnit,
                                 i+1, gmH, im1, dbg_G, 1 );

                    solver.copy( pastix_data, n, b, dbg_r );
                    solver.copy( pastix_data, n, x, dbg_x );

                    /* Accumulate the current v_m */
                    solver.gemv( pastix_data, n, i+1, 1.0, (precond ? gmW : gmV), n, dbg_G, 1.0, dbg_x );

                    /* Compute b - Ax */
                    solver.spmv( pastix_data, PastixNoTrans, -1., dbg_x, 1., dbg_r );

                    normr2 = solver.norm( pastix_data, n, dbg_r );
                    fprintf(stdout, OUT_ITERREFINE_ERR, normr2 / normb );
                }
#endif
            }
        }

        /* Compute y_m = H_m^{-1} g_m (See p. 169) */
        cblas_ztrsv( CblasColMajor, CblasUpper, CblasNoTrans, CblasNonUnit,
                     i+1, gmH, im1, gmG, 1 );

        /**
         * Compute x_m = x_0 + M^{-1} V_m y_m
         *             = x_0 +        W_m y_m
         */
        if (precond && savemem) {
            /**
             * Since we saved memory, we do not have (M^{-1} V_m) stored,
             * thus we compute:
             *     w = V_m y_m
             *     w = M^{-1} (V_m y_m)
             *     x = x0 + (M^{-1} (V_m y_m))
             */
            solver.gemv( pastix_data, n, i+1, 1.0, gmV, n, gmG, 0., gmW );
            solver.spsv( pastix_data, gmW );
            solver.axpy( pastix_data, n, 1.,  gmW, x );
        }
        else {
            /**
             * Since we did not saved memory, we do have (M^{-1} V_m) stored in
             * W_m if precond is true, thus we compute:
             *     x = x0 + W_m y_m, if precond
             *     x = x0 + V_m y_m, if not precond
             */
            gmWi = precond ? gmW : gmV;
            solver.gemv( pastix_data, n, i+1, 1.0, gmWi, n, gmG, 1.0, x );
        }

        if ((resid_b <= eps) || (iters >= itermax))
        {
            outflag = 0;
        }
    }

    clockStop( refine_clk );
    t3 = clockGet();

    solver.output_final( pastix_data, resid_b, iters, t3, x, x );

    solver.free(gmcos);
    solver.free(gmsin);
    solver.free(gmG);
    solver.free(gmH);
    solver.free(gmV);
    solver.free(gmW);
#if defined(PASTIX_DEBUG_GMRES)
    solver.free(dbg_x);
    solver.free(dbg_r);
    solver.free(dbg_G);
#endif

    return iters;
}
