/**
 *
 * @file spm_tests.h
 *
 * SParse Matrix package testings header.
 *
 * @copyright 2016-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @date 2013-06-24
 *
 **/
#ifndef _spm_tests_h_
#define _spm_tests_h_

#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <spm.h>

void spmGetOptions( int argc, char **argv,
                    spm_driver_t *driver, char **filename );

void core_zplrnt( int m, int n, spm_complex64_t *A, int lda,
                  int gM, int m0, int n0, unsigned long long int seed );
int  core_zgeadd( spm_trans_t            trans,
                  spm_int_t              M,
                  spm_int_t              N,
                  spm_complex64_t        alpha,
                  const spm_complex64_t *A,
                  spm_int_t              LDA,
                  spm_complex64_t        beta,
                  spm_complex64_t       *B,
                  spm_int_t              LDB );

void core_cplrnt( int m, int n, spm_complex32_t *A, int lda,
                  int gM, int m0, int n0, unsigned long long int seed );
int  core_cgeadd( spm_trans_t            trans,
                  spm_int_t              M,
                  spm_int_t              N,
                  spm_complex32_t        alpha,
                  const spm_complex32_t *A,
                  spm_int_t              LDA,
                  spm_complex32_t        beta,
                  spm_complex32_t       *B,
                  spm_int_t              LDB );

void core_dplrnt( int m, int n, double *A, int lda,
                  int gM, int m0, int n0, unsigned long long int seed );
int  core_dgeadd( spm_trans_t   trans,
                  spm_int_t     M,
                  spm_int_t     N,
                  double        alpha,
                  const double *A,
                  spm_int_t     LDA,
                  double        beta,
                  double       *B,
                  spm_int_t     LDB );

void core_splrnt( int m, int n, float *A, int lda,
                  int gM, int m0, int n0, unsigned long long int seed );
int  core_sgeadd( spm_trans_t  trans,
                  spm_int_t    M,
                  spm_int_t    N,
                  float        alpha,
                  const float *A,
                  spm_int_t    LDA,
                  float        beta,
                  float       *B,
                  spm_int_t    LDB );

void z_spm_print_check( char *filename, const spmatrix_t *spm );
int  z_spm_matvec_check( spm_trans_t trans, const spmatrix_t *spm );
int  z_spm_norm_check( const spmatrix_t *spm );

void c_spm_print_check( char *filename, const spmatrix_t *spm );
int  c_spm_matvec_check( spm_trans_t trans, const spmatrix_t *spm );
int  c_spm_norm_check( const spmatrix_t *spm );

void d_spm_print_check( char *filename, const spmatrix_t *spm );
int  d_spm_matvec_check( spm_trans_t trans, const spmatrix_t *spm );
int  d_spm_norm_check( const spmatrix_t *spm );

void s_spm_print_check( char *filename, const spmatrix_t *spm );
int  s_spm_matvec_check( spm_trans_t trans, const spmatrix_t *spm );
int  s_spm_norm_check( const spmatrix_t *spm );

static inline int
spm_norm_print_result( double norms, double normd, double result )
{
    int rc = 0;
    if ( (result >= 0.) && (result < 1.) ) {
        printf("SUCCESS !\n");
    } else {
        printf("FAILED !\n");
        rc=1;
    }

    printf("   Nsparse = %e, Ndense = %e\n", norms, normd );
    printf("  | Nsparse - Ndense | / Ndense = %e\n", result);

    return rc;
}

#endif /* _spm_tests_h_ */
