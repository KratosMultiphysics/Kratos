/**
 *
 * @file core_zgetro.c
 *
 * PaStiX kernel routines
 *
 * @copyright 2010-2015 Univ. of Tennessee, Univ. of California Berkeley and
 *                      Univ. of Colorado Denver. All rights reserved.
 * @copyright 2015-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.1
 * @author Mathieu Faverge
 * @date 2018-07-16
 * @precisions normal z -> c d s
 *
 **/
#include "common.h"

/**
 ******************************************************************************
 *
 * @brief Transposes a m-by-n matrix out of place using an extra workspace of size
 * m-by-n.
 *
 *******************************************************************************
 *
 * @param[in] m
 *         Number of rows of A.
 *
 * @param[in] n
 *         Number of columns of A.
 *
 * @param[in] A
 *         Matrix to be transposed.
 *
 * @param[in] lda
 *         Leading dimension of matrix A.
 *
 * @param[inout] B
 *         On exit B = trans(A).
 *
 * @param[in] ldb
 *         Leading dimension of matrix B.
 *
 ******************************************************************************/
void
core_zgetro( int m, int n,
             const pastix_complex64_t *A, int lda,
             pastix_complex64_t *B, int ldb )
{
    int i, j;

    /* rectangular transposition (use workspace) */
    for (i=0; i<m; i++) {
        for (j=0; j<n; j++) {
            B[j+i*ldb] = A[i+j*lda];
        }
    }
}
