/**
 *
 * @file z_ge2lr_tests.c
 *
 * Tests and validate the Xge2lr routine.
 *
 * @copyright 2015-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.1
 * @author Gregoire Pichon
 * @date 2018-07-16
 *
 * @precisions normal z -> z c d s
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
z_ge2lr_test( int mode, double tolerance, pastix_int_t rank,
              pastix_int_t m, pastix_int_t n, pastix_int_t lda )
{

    pastix_complex64_t *A, *A2;
    pastix_lrblock_t    lrA_rrqr, lrA_svd;

    double norm_dense;
    double norm_diff_RRQR, norm_diff_SVD;
    double res_SVD, res_RRQR;

    pastix_int_t rankmax  = core_get_rklimit(m, n);
    int          rc = 0;

    /*
     * Generate a matrix of rank rank and its compressed SVD/RRQR versions
     */
    z_lowrank_genmat( mode, tolerance, rank, m, n, lda,
                      &A, &lrA_svd, &lrA_rrqr, &norm_dense );

    printf(" The rank of A is: RRQR %d SVD %d rkmax %d\n", lrA_rrqr.rk, lrA_svd.rk, (int)rankmax);

    /*
     * Check || A - c(A) || < tol * || A ||
     */
    A2 = malloc(n * lda * sizeof(pastix_complex64_t));

    /* SVD */
    {
        core_zlr2ge( PastixNoTrans, m, n,
                     &lrA_svd, A2, lda );

        core_zgeadd( PastixNoTrans, m, n,
                     -1., A,  lda,
                      1., A2, lda );

        norm_diff_SVD = LAPACKE_zlange_work( LAPACK_COL_MAJOR, 'f', m, n,
                                             A2, lda, NULL );
    }

    /* RRQR */
    {
        core_zlr2ge( PastixNoTrans, m, n,
                     &lrA_svd, A2, lda );

        core_zgeadd( PastixNoTrans, m, n,
                     -1., A,  lda,
                      1., A2, lda );

        norm_diff_RRQR = LAPACKE_zlange_work( LAPACK_COL_MAJOR, 'f', m, n,
                                              A2, lda, NULL );
    }

    if (rank != 0){
        res_RRQR = norm_diff_RRQR / ( tolerance * norm_dense );
        res_SVD  = norm_diff_SVD  / ( tolerance * norm_dense );
    }
    else{
        res_RRQR = norm_diff_RRQR;
        res_SVD  = norm_diff_SVD;
    }

    free(A);
    free(A2);

    /* Check the correctness of the compression */
    if (res_RRQR > 10.0) {
        rc += 1;
    }
    if (res_SVD > 10.0) {
        rc += 2;
    }

    /* Check that SVD rank is equal to the desired rank */
    if ( ((rank >  rankmax) && (lrA_svd.rk != -1  )) ||
         ((rank <= rankmax) && ((lrA_svd.rk < (rank-2)) || (lrA_svd.rk > (rank+2)))) )
    {
        rc += 4;
    }

    /* Check that RRQR rank is larger or equal to SVD rank */
    if (lrA_svd.rk == -1) {
        if (lrA_rrqr.rk != -1) {
            rc += 8;
        }
    }
    else {
        if ( (lrA_rrqr.rk != -1) &&
             ((lrA_rrqr.rk < lrA_svd.rk) || (lrA_rrqr.rk > (lrA_svd.rk + 1.25 * rank ))) )
        {
            rc += 16;
        }
    }

    core_zlrfree( &lrA_rrqr );
    core_zlrfree( &lrA_svd );
    return rc;
}

int main (int argc, char **argv)
{
    (void) argc;
    (void) argv;
    int err = 0;
    int ret;
    pastix_int_t m, r;
    double eps = LAPACKE_dlamch_work('e');
    double tolerance = sqrt(eps);

    for (m=100; m<300; m+=100){
        for (r=0; r <= (m/2); r += ( r + 1 ) ) {
            printf("   -- Test GE2LR M=N=LDA=%ld R=%ld\n", (long)m, (long)r);

            ret = z_ge2lr_test(0, tolerance, r, m, m, m);
            PRINT_RES(ret);
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
