/**
 *
 * @file z_lrmm_tests.c
 *
 * Tests and validate the core_zlrmm() routine.
 *
 * @copyright 2015-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.1
 * @author Gregoire Pichon
 * @date 2018-07-16
 *
 * @precisions normal z -> c d s
 *
 **/
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include <pastix.h>
#include "common/common.h"
#include <lapacke.h>
#include <cblas.h>
#include "blend/solver.h"
#include "kernels/pastix_zcores.h"
#include "kernels/pastix_zlrcores.h"

#define PRINT_RES(_ret_)                        \
    if(_ret_ == -1) {                           \
        printf("UNDEFINED\n");                  \
    }                                           \
    else if(_ret_ > 0) {                        \
        printf("FAILED(%d)\n", _ret_);          \
        err++;                                  \
    }                                           \
    else {                                      \
        printf("SUCCESS\n");                    \
    }

int
z_lowrank_genmat( int mode, double tolerance, pastix_int_t rank,
                  pastix_int_t m, pastix_int_t n, pastix_int_t lda,
                  pastix_complex64_t **Aptr,
                  pastix_lrblock_t    *lrA_svd,
                  pastix_lrblock_t    *lrA_rrqr,
                  double              *normA );

int
z_lowrank_lrmm( pastix_int_t m, pastix_int_t n, pastix_int_t k,
                const pastix_lrblock_t *lrA, double normA,
                const pastix_lrblock_t *lrB, double normB,
                pastix_int_t Cm, pastix_int_t Cn, pastix_int_t offx, pastix_int_t offy,
                const pastix_complex64_t *C, pastix_lrblock_t *lrC, double normC,
                pastix_lr_t *lowrank );

extern pastix_lr_t z_lr_svd;
extern pastix_lr_t z_lr_rrqr;

int main (int argc, char **argv)
{
    (void) argc;
    (void) argv;
    int err = 0;
    int ret, rc, mode = 0;
    pastix_int_t m, n, k, Cm, Cn, offx, offy;
    pastix_int_t compress_type;
    double       eps = LAPACKE_dlamch_work('e');
    double       tolerance = sqrt(eps);
    pastix_complex64_t *A, *B, *C;
    pastix_lrblock_t    lrA_svd, lrA_rrqr;
    pastix_lrblock_t    lrB_svd, lrB_rrqr;
    pastix_lrblock_t    lrC_svd, lrC_rrqr;
    double              norm_dense_A, norm_dense_B, norm_dense_C;
    int ranks[3], rA, rB, rC, s, i, j, l;

    pastix_complex64_t mzone = -1.0;
    pastix_complex64_t zone  = 1.0;

    compress_type = 3;
    if ( argc > 1 ) {
        compress_type = atoi( argv[1] );
    }

    /* Initialize tolerance of low-rank structure */
    z_lr_rrqr.tolerance = tolerance;
    z_lr_svd.tolerance  = tolerance;

    for (s=100; s<=200; s = 2*s) {
        ranks[0] = s + 1;
        ranks[1] = 16;
        ranks[2] = 2;

        m = s / 2;
        n = s / 4;
        k = s / 6;

        offx = 1;
        offy = 2;

        Cm = s;
        Cn = s;

        /* Matrix A */
        for (i=0; i<3; i++) {
            rA = pastix_imin( m, k ) / ranks[i];

            z_lowrank_genmat( mode, tolerance, rA, m, k, m,
                              &A, &lrA_svd, &lrA_rrqr, &norm_dense_A );

            /* Matrix B */
            for (j=0; j<3; j++) {
                rB = pastix_imin( n, k ) / ranks[j];

                z_lowrank_genmat( mode, tolerance, rB, n, k, n,
                                  &B, &lrB_svd, &lrB_rrqr, &norm_dense_B );

                /* Matrix C */
                for (l=0; l<3; l++) {
                    rC = pastix_imin( Cm, Cn ) / ranks[l];

                    z_lowrank_genmat( mode, tolerance, rC, Cm, Cn, Cm,
                                      &C, &lrC_svd, &lrC_rrqr, &norm_dense_C );


                    printf( "  -- Test LRMM Cm=%ld, Cn=%ld, m=%ld, n=%ld, k=%ld, rA=%ld, rB=%ld, rC=%ld\n",
                            (long)Cm, (long)Cn, (long)m, (long)n, (long)k, (long)rA, (long)rB, (long)rC );

                    /* Compute the full rank GEMM */
                    cblas_zgemm( CblasColMajor, CblasNoTrans, CblasConjTrans,
                                 m, n, k,
                                 CBLAS_SADDR( mzone ), A, m,
                                                       B, n,
                                 CBLAS_SADDR( zone ),  C + Cm * offy + offx, Cm );

                    /*
                     * Generate matrices of rankA and rankB and their compressed SVD/RRQR versions
                     */
                    ret = 0;
                    if (compress_type & 1)
                    {
                        printf("     RRQR: A %2d B %2d C %2d ", lrA_rrqr.rk, lrB_rrqr.rk, lrC_rrqr.rk );
                        rc = z_lowrank_lrmm( m, n, k,
                                             &lrA_rrqr, norm_dense_A,
                                             &lrB_rrqr, norm_dense_B,
                                             Cm, Cn, offx, offy,
                                             C, &lrC_rrqr, norm_dense_C,
                                             &z_lr_rrqr );
                        assert( rc >= 0 );
                        ret += rc;
                    }
                    core_zlrfree( &lrC_rrqr );

                    if (compress_type & 2)
                    {
                        printf("     SVD:  A %2d B %2d C %2d ", lrA_svd.rk, lrB_svd.rk, lrC_svd.rk );
                        rc = z_lowrank_lrmm( m, n, k,
                                             &lrA_svd, norm_dense_A,
                                             &lrB_svd, norm_dense_B,
                                             Cm, Cn, offx, offy,
                                             C, &lrC_svd, norm_dense_C,
                                             &z_lr_svd );
                        assert( rc >= 0 );
                        ret += rc;
                    }
                    core_zlrfree( &lrC_svd );
                    free( C );

                    PRINT_RES(ret);
                }

                free( B );
                core_zlrfree( &lrB_rrqr );
                core_zlrfree( &lrB_svd );
            }
            free( A );
            core_zlrfree( &lrA_rrqr );
            core_zlrfree( &lrA_svd );
        }
    }

    if( err == 0 ) {
        printf(" -- All tests PASSED --\n");
        return EXIT_SUCCESS;
    }
    else
    {
        printf(" -- %d tests FAILED --\n", err);
        return EXIT_FAILURE;
    }
}
