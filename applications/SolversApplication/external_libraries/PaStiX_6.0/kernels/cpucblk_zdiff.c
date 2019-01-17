/**
 *
 * @file cpucblk_zdiff.c
 *
 * Precision dependent routines to differentiate two solver matrix structures
 * when debuging.
 *
 * @copyright 2015-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.1
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @author Gregoire Pichon
 * @author Mathieu Faverge
 * @date 2018-07-16
 *
 * @precisions normal z -> s d c
 *
 **/
#include "common.h"
#include "solver.h"
#include "lapacke.h"
#include "pastix_zcores.h"
#include "pastix_zlrcores.h"

/**
 *******************************************************************************
 *
 * @brief Compare two column blocks in full-rank format.
 *
 * The second cblk is overwritten by the difference of the two column blocks.
 * The frobenius norm of the difference is computed and the functions returns 0
 * if the result:
 *      || B - A || / ( || A || * eps )
 *
 * is below 10. Otherwise, an error message is printed and 1 is returned.
 *
 *******************************************************************************
 *
 * @param[in] side
 *          Define which side of the cblk must be tested.
 *          @arg PastixLCoef if lower part only
 *          @arg PastixUCoef if upper part only
 *          @arg PastixLUCoef if both sides.
 *
 * @param[in] cblkA
 *          The column block of the A matrix.
 *
 * @param[inout] cblkB
 *          The column block of the B matrix that matches the A matrix in
 *          stucture.
 *          On exit, cblkB coefficient arrays are overwritten by the result of
 *          (B-A).
 *
 *******************************************************************************
 *
 * @return 0 if the test is passed, >= 0 otherwise.
 *
 *******************************************************************************/
int
cpucblk_zdiff( pastix_coefside_t side,
               const SolverCblk *cblkA,
               SolverCblk       *cblkB )
{
    pastix_complex64_t *coefA;
    pastix_complex64_t *coefB;
    pastix_int_t        ncols  = cblk_colnbr( cblkA );
    pastix_int_t        stride = cblkA->stride;
    double normdiff, normfull, normlowr, res, eps;
    int rc = 0;

    assert( ncols  == cblk_colnbr( cblkB ) );
    assert( stride == cblkB->stride );
    assert( !(cblkA->cblktype & CBLK_LAYOUT_2D) ); /* Not yet implemented */

    eps = LAPACKE_dlamch_work( 'e' );

    if ( side != PastixUCoef ) {
        coefA = cblkA->lcoeftab;
        coefB = cblkB->lcoeftab;

        assert( (coefA != NULL) && (coefB != NULL) );

        normfull = LAPACKE_zlange_work( LAPACK_COL_MAJOR, 'f', stride, ncols,
                                        coefA, stride, NULL );
        normlowr = LAPACKE_zlange_work( LAPACK_COL_MAJOR, 'f', stride, ncols,
                                        coefB, stride, NULL );
        core_zgeadd( PastixNoTrans, stride, ncols,
                     -1., coefA, stride,
                      1., coefB, stride );

        normdiff = LAPACKE_zlange( LAPACK_COL_MAJOR, 'M', stride, ncols,
                                   coefB, stride );
        res = (normfull == 0.) ? 0. : (normdiff / (normfull * eps));

        if ( res > 10 ) {
            fprintf(stderr, "KO on L: ||full(A)||_f=%e, ||comp(A)||_f=%e, ||comp(A)-full(A)||_0=%e, ||comp(A)-full(A)||_0 / (||full(A)||_2 * eps)=%e\n",
                    normfull, normlowr, normdiff, res );
            rc++;
        }
    }

    if ( side != PastixLCoef ) {
        coefA = cblkA->ucoeftab;
        coefB = cblkB->ucoeftab;

        assert( (coefA != NULL) && (coefB != NULL) );

        normfull = LAPACKE_zlange_work( LAPACK_COL_MAJOR, 'f', stride, ncols,
                                        coefA, stride, NULL );
        normlowr = LAPACKE_zlange_work( LAPACK_COL_MAJOR, 'f', stride, ncols,
                                        coefB, stride, NULL );
        core_zgeadd( PastixNoTrans, stride, ncols,
                     -1., coefA, stride,
                      1., coefB, stride );

        normdiff = LAPACKE_zlange_work( LAPACK_COL_MAJOR, 'M', stride, ncols,
                                        coefB, stride, NULL );
        res = (normfull == 0.) ? 0. : (normdiff / (normfull * eps));

        if ( res > 10 ) {
            fprintf(stderr, "KO on L: ||full(A)||_f=%e, ||comp(A)||_f=%e, ||comp(A)-full(A)||_0=%e, ||comp(A)-full(A)||_0 / (||full(A)||_2 * eps)=%e\n",
                    normfull, normlowr, normdiff, res );
            rc++;
        }
    }

    return rc;
}
