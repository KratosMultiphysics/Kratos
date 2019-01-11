/*! \file
Copyright (c) 2003, The Regents of the University of California, through
Lawrence Berkeley National Laboratory (subject to receipt of any required 
approvals from U.S. Dept. of Energy) 

All rights reserved. 

The source code is distributed under BSD license, see the file License.txt
at the top-level directory.
*/

#include "slu_mt_cdefs.h"

int_t pcgst02(trans_t trans, int_t m, int_t n, int_t nrhs, SuperMatrix *A,
	    complex *x, int_t ldx, complex *b, int_t ldb, float *resid)
{
/*
 * -- SuperLU MT routine (version 2.0) --
 * Lawrence Berkeley National Lab, Univ. of California Berkeley,
 * and Xerox Palo Alto Research Center.
 * September 10, 2007
 *
 *  Purpose   
 *  =======   
 *
 *  pcgst02() computes the residual for a solution of a system of linear   
 *  equations  A*x = b  or  A'*x = b:   
 *      RESID = norm(B - A*X) / ( norm(A) * norm(X) * EPS ),   
 *  where EPS is the machine epsilon.   
 *
 *  Arguments   
 *  =========   
 *
 *  TRANS   (input) trans_t
 *          Specifies the form of the system of equations:   
 *          = NOTRANS: A *x = b   
 *          = TRANS:   A'*x = b, where A' is the transpose of A   
 *          = CONJ:    A'*x = b, where A' is the conjugate transpose of A   
 *
 *  M       (input) INTEGER   
 *
 *  N       (input) INTEGER   
 *          The number of columns of the matrix A.  N >= 0.   
 *
 *  NRHS    (input) INTEGER   
 *          The number of columns of B, the matrix of right hand sides.   
 *          NRHS >= 0.
 *
 *  A       (input) SuperMatrix*, dimension (LDA,N)   
 *          The original M x N sparse matrix A.   
 *
 *  X       (input) DOUBLE PRECISION array, dimension (LDX,NRHS)   
 *          The computed solution vectors for the system of linear   
 *          equations.   
 *
 *  LDX     (input) INTEGER   
 *          The leading dimension of the array X.  If TRANS = NOTRANS,   
 *          LDX >= max(1,N); if TRANS = TRANS or CONJ, LDX >= max(1,M).   
 *
 *  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)   
 *          On entry, the right hand side vectors for the system of   
 *          linear equations.   
 *          On exit, B is overwritten with the difference B - A*X.   
 *
 *  LDB     (input) INTEGER   
 *          The leading dimension of the array B.  IF TRANS = NOTRANS,
 *          LDB >= max(1,M); if TRANS = TRANS or CONJ, LDB >= max(1,N).
 *
 *  RESID   (output) DOUBLE PRECISION   
 *          The maximum over the number of right hand sides of   
 *
 *  =====================================================================
*/

    /* Table of constant values */
    complex alpha = {-1., 0.0};
    complex beta  = {1., 0.0};
    int    c__1  = 1;
    
    /* System generated locals */
    float d__1, d__2;

    /* Local variables */
    int_t j;
    int n1, n2;
    float anorm, bnorm;
    float xnorm;
    float eps;
    char transc[1];

    /* Function prototypes */
    extern int lsame_(char *, char *);
    extern float clangs(char *, SuperMatrix *);
    extern float scasum_(int *, complex *, int *);
    extern double slamch_(char *);
    
    /* Function Body */
    if ( m <= 0 || n <= 0 || nrhs == 0) {
	*resid = 0.;
	return 0;
    }

    if ( trans == TRANS || trans == CONJ ) {
	n1 = n;
	n2 = m;
    } else {
	n1 = m;
	n2 = n;
    }

    /* Exit with RESID = 1/EPS if ANORM = 0. */

    eps = slamch_("Epsilon");
    anorm = clangs("1", A);
    if (anorm <= 0.) {
	*resid = 1. / eps;
	return 0;
    }

    /* Compute  B - A*X  (or  B - A'*X ) and store in B. */

    if ( trans == NOTRANS ) *transc = 'N';
    else if ( trans == TRANS ) *transc = 'T';
    else if ( trans == CONJ ) *transc = 'C';

    sp_cgemm(transc, n1, nrhs, n2, alpha, A, x, ldx, beta, b, ldb);

    /*for (j = 0; j < m; ++j)
      if ( b[j] > 0.001 ) { printf("b-Ax: %d, %f\n", j, b[j]); b[j] = 1.; }*/
 
    /* Compute the maximum over the number of right hand sides of   
       norm(B - A*X) / ( norm(A) * norm(X) * EPS ) . */

    *resid = 0.;
    for (j = 0; j < nrhs; ++j) {
        bnorm = scasum_(&n1, &b[j*ldb], &c__1);
        xnorm = scasum_(&n2, &x[j*ldx], &c__1);
	if (xnorm <= 0.) {
	    *resid = 1. / eps;
	} else {
	    /* Computing MAX */
	    d__1 = *resid, d__2 = bnorm / anorm / xnorm / eps;
	    *resid = SUPERLU_MAX(d__1, d__2);
	}
    }

    return 0;

} /* pcgst02 */

