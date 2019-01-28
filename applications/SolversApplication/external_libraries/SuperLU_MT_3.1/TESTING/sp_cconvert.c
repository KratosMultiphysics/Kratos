/*! \file
Copyright (c) 2003, The Regents of the University of California, through
Lawrence Berkeley National Laboratory (subject to receipt of any required 
approvals from U.S. Dept. of Energy) 

All rights reserved. 

The source code is distributed under BSD license, see the file License.txt
at the top-level directory.
*/

#include "slu_mt_cdefs.h"

int_t
sp_cconvert(int_t m, int_t n, complex *A, int_t lda, int_t kl, int_t ku,
	   complex *a, int_t *asub, int_t *xa, int_t *nnz)
{
/*
 * -- SuperLU MT routine (version 1.0) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * August 15, 1997
 *
 * Convert a full matrix into a sparse matrix format. 
 */
    int_t     lasta = 0;
    int_t     i, j, ilow, ihigh;
    int_t     *row;
    complex  *val;

    for (j = 0; j < n; ++j) {
	xa[j] = lasta;
	val = &a[xa[j]];
	row = &asub[xa[j]];

	ilow = SUPERLU_MAX(0, j - ku);
	ihigh = SUPERLU_MIN(n-1, j + kl);
	for (i = ilow; i <= ihigh; ++i) {
	    val[i-ilow] = A[i + j*lda];
	    row[i-ilow] = i;
	}
	lasta += ihigh - ilow + 1;
    }

    xa[n] = *nnz = lasta;
    return 0;
}


