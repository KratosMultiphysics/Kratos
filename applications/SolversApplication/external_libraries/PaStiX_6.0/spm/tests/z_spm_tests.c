/**
 *
 * @file z_spm_tests.c
 *
 * Tests and validate the spm_convert routines.
 *
 * @copyright 2015-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.0
 * @author Mathieu Faverge
 * @author Theophile Terraz
 * @date 2015-01-01
 *
 * @precisions normal z -> c d s
 *
 **/
#ifndef _GNU_SOURCE
#define _GNU_SOURCE 1
#endif
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <spm_tests.h>
#include "cblas.h"
#include "lapacke.h"
#include <z_spm.h>

void
core_zplrnt( int m, int n, spm_complex64_t *A, int lda,
             int gM, int m0, int n0, unsigned long long int seed );
int
core_zgeadd( spm_trans_t            trans,
             spm_int_t              M,
             spm_int_t              N,
             spm_complex64_t        alpha,
             const spm_complex64_t *A,
             spm_int_t              LDA,
             spm_complex64_t        beta,
             spm_complex64_t       *B,
             spm_int_t              LDB);

/*------------------------------------------------------------------------
 *  Check the accuracy of the solution
 */
void
z_spm_print_check( char *filename, const spmatrix_t *spm )
{
    spm_complex64_t *A;
    char *file;
    FILE *f;
    int rc;

    rc = asprintf( &file, "expand_%s_sparse_cp.dat", filename );
    if ( (f = fopen( file, "w" )) == NULL ) {
        perror("z_spm_print_check:sparse_cp");
        return;
    }
    z_spmPrint( f, spm );
    fclose(f);
    free(file);

    A = z_spm2dense( spm );
    rc = asprintf( &file, "expand_%s_dense_cp.dat", filename );
    if ( (f = fopen( file, "w" )) == NULL ) {
        perror("z_spm_print_check:dense_cp");
        return;
    }
    z_spmDensePrint( f, spm->nexp, spm->nexp, A, spm->nexp );
    fclose(f);
    free(file);
    free(A);

    if ( spm->dof != 1 ) {
        spmatrix_t *espm = z_spmExpand( spm );

        rc = asprintf( &file, "expand_%s_sparse_ucp.dat", filename );
        if ( (f = fopen( file, "w" )) == NULL ) {
            perror("z_spm_print_check:sparse_ucp");
            return;
        }
        z_spmPrint( f, espm );
        fclose(f);
        free(file);

        A = z_spm2dense( espm );
        rc = asprintf( &file, "expand_%s_dense_ucp.dat", filename );
        if ( (f = fopen( file, "w" )) == NULL ) {
            perror("z_spm_print_check:dense_ucp");
            return;
        }
        z_spmDensePrint( f, espm->nexp, espm->nexp, A, espm->nexp );
        fclose(f);
        free(file);
        free(A);

        spmExit( espm );
        free(espm);
    }

    (void)rc; (void)A;
    return;
}

/*------------------------------------------------------------------------
 *  Check the accuracy of the solution
 */
int
z_spm_matvec_check( spm_trans_t trans, const spmatrix_t *spm )
{
    unsigned long long int seed = 35469;
    spm_complex64_t *A, *x, *y0, *ys, *yd;
    /*
     * Alpha and beta are complex for cblas, but only the real part is used for
     * matvec/matmat subroutines
     */
    spm_complex64_t zalpha = 0.;
    spm_complex64_t zbeta  = 0.;
    double dalpha = 0.;
    double dbeta  = 0.;

    double Anorm, Xnorm, Y0norm, Ysnorm, Ydnorm, Rnorm;
    double eps, result;
    int rc, info_solution, start = 1;

    eps = LAPACKE_dlamch_work('e');

    core_dplrnt( 1, 1, &dalpha, 1, 1, start, 0, seed ); start++;
    core_dplrnt( 1, 1, &dbeta,  1, 1, start, 0, seed ); start++;

    /* Make sure alpha/beta are doubles */
    zalpha = dalpha;
    zbeta = dbeta;

    x = (spm_complex64_t*)malloc(spm->gNexp * sizeof(spm_complex64_t));
    core_zplrnt( spm->gNexp, 1, x, spm->gNexp, 1, start, 0, seed ); start += spm->gNexp;

    y0 = (spm_complex64_t*)malloc(spm->gNexp * sizeof(spm_complex64_t));
    core_zplrnt( spm->gNexp, 1, y0, spm->gNexp, 1, start, 0, seed ); start += spm->gNexp;

    /* Create a dense backup of spm */
    A = z_spm2dense( spm );

    /* Allocate cs/cd */
    ys = (spm_complex64_t*)malloc(spm->gNexp * sizeof(spm_complex64_t));
    yd = (spm_complex64_t*)malloc(spm->gNexp * sizeof(spm_complex64_t));

    /* Initialize cs/cd */
    memcpy( ys, y0, spm->gNexp * sizeof(spm_complex64_t) );
    memcpy( yd, y0, spm->gNexp * sizeof(spm_complex64_t) );

    /* Compute the sparse matrix-vector product */
    //rc = spmMatVec( trans, dalpha, spm, x, dbeta, ys );
    rc = spmMatMat( trans, 1, dalpha, spm, x, spm->nexp, dbeta, ys, spm->nexp );
    if ( rc != SPM_SUCCESS ) {
        info_solution = 1;
        goto end;
    }

    /* Compute the dense matrix-vector product */
    cblas_zgemm( CblasColMajor, (CBLAS_TRANSPOSE)trans, CblasNoTrans,
                 spm->gNexp, 1, spm->gNexp,
                 CBLAS_SADDR(zalpha), A, spm->gNexp,
                                      x, spm->gNexp,
                 CBLAS_SADDR(zbeta), yd, spm->gNexp );

    Anorm  = LAPACKE_zlange( LAPACK_COL_MAJOR, 'I', spm->gNexp, spm->gNexp,  A, spm->gNexp );
    Xnorm  = LAPACKE_zlange( LAPACK_COL_MAJOR, 'I', spm->gNexp, 1,           x, spm->gNexp );
    Y0norm = LAPACKE_zlange( LAPACK_COL_MAJOR, 'I', spm->gNexp, 1,          y0, spm->gNexp );
    Ysnorm = LAPACKE_zlange( LAPACK_COL_MAJOR, 'I', spm->gNexp, 1,          ys, spm->gNexp );
    Ydnorm = LAPACKE_zlange( LAPACK_COL_MAJOR, 'I', spm->gNexp, 1,          yd, spm->gNexp );

    core_zgeadd(SpmNoTrans, spm->gNexp, 1,
                -1., ys, spm->gNexp,
                 1., yd, spm->gNexp);
    Rnorm = LAPACKE_zlange( LAPACK_COL_MAJOR, 'M', spm->gNexp, 1, yd, spm->gNexp );

    if ( 1 ) {
        printf("  ||A||_inf = %e, ||x||_inf = %e, ||y||_inf = %e\n"
               "  ||dense(a*A*x+b*y)||_inf = %e, ||sparse(a*A*x+b*y)||_inf = %e, ||R||_m = %e\n",
               Anorm, Xnorm, Y0norm, Ydnorm, Ysnorm, Rnorm);
    }

    result = Rnorm / ((Anorm + Xnorm + Y0norm) * spm->gNexp* eps);
    if (  isinf(Ydnorm) || isinf(Ysnorm) ||
          isnan(result) || isinf(result) || (result > 10.0) ) {
        info_solution = 1;
    }
    else {
        info_solution = 0;
    }

  end:
    free(A); free(x); free(y0); free(ys); free(yd);

    return info_solution;
}

/*------------------------------------------------------------------------
 *  Check the accuracy of the solution
 */
int
z_spm_norm_check( const spmatrix_t *spm )
{
    spm_complex64_t *A;
    double norms, normd;
    double eps, result;
    int ret = 0;

    eps = LAPACKE_dlamch_work('e');

    /* Create a dense backup of spm */
    A = z_spm2dense( spm );

    /**
     * Test Norm Max
     */
    printf(" -- Test norm Max :");
    norms = spmNorm( SpmMaxNorm, spm );
    normd = LAPACKE_zlange( LAPACK_COL_MAJOR, 'M', spm->gNexp, spm->gNexp, A, spm->gNexp );
    result = fabs(norms - normd) / (normd * eps);
    ret += spm_norm_print_result( norms, normd, result );

    /**
     * Test Norm Inf
     */
    printf(" -- Test norm Inf :");
    norms = spmNorm( SpmInfNorm, spm );
    normd = LAPACKE_zlange( LAPACK_COL_MAJOR, 'I', spm->gNexp, spm->gNexp, A, spm->gNexp );
    result = fabs(norms - normd) / (normd * eps);
    result = result * ((double)(spm->gNexp)) / ((double)(spm->gnnzexp));
    ret += spm_norm_print_result( norms, normd, result );

    /**
     * Test Norm One
     */
    printf(" -- Test norm One :");
    norms = spmNorm( SpmOneNorm, spm );
    normd = LAPACKE_zlange( LAPACK_COL_MAJOR, 'O', spm->gNexp, spm->gNexp, A, spm->gNexp );
    result = fabs(norms - normd) / (normd * eps);
    result = result * ((double)(spm->gNexp)) / ((double)(spm->gnnzexp));
    ret += spm_norm_print_result( norms, normd, result );

    /**
     * Test Norm Frobenius
     */
    printf(" -- Test norm Frb :");
    norms = spmNorm( SpmFrobeniusNorm, spm );
    normd = LAPACKE_zlange( LAPACK_COL_MAJOR, 'F', spm->gNexp, spm->gNexp, A, spm->gNexp );
    result = fabs(norms - normd) / (normd * eps);
    result = result / ((double)spm->gnnzexp);
    ret += spm_norm_print_result( norms, normd, result );

    free(A);
    return ret;
}
