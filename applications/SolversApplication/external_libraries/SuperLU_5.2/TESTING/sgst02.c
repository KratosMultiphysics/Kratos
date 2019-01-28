/*! \file
Copyright (c) 2003, The Regents of the University of California, through
Lawrence Berkeley National Laboratory (subject to receipt of any required 
approvals from U.S. Dept. of Energy) 

All rights reserved. 

The source code is distributed under BSD license, see the file License.txt
at the top-level directory.
*/

/*
 * -- SuperLU routine (version 3.0) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * October 15, 2003
 *
 */
#include "slu_sdefs.h"

int sgst02(trans_t trans, int m, int n, int nrhs, SuperMatrix *A,
	      float *x, int ldx, float *b, int ldb, float *resid)
{
/*  
    Purpose   
    =======   

    SGST02 computes the residual for a solution of a system of linear   
    equations  A*x = b  or  A'*x = b:   
       RESID = norm(B - A*X) / ( norm(A) * norm(X) * EPS ),   
    where EPS is the machine epsilon.   

    Arguments   
    =========   

    TRANS   (input) trans_t
            Specifies the form of the system of equations:   
            = NOTRANS:  A *x = b   
            = TRANS  :  A'*x = b, where A' is the transpose of A   
            = CONJ   :  A'*x = b, where A' is the transpose of A   

    M       (input) INTEGER   
            The number of rows of the matrix A.  M >= 0.   

    N       (input) INTEGER   
            The number of columns of the matrix A.  N >= 0.   

    NRHS    (input) INTEGER   
            The number of columns of B, the matrix of right hand sides.   
            NRHS >= 0.
	    
    A       (input) SuperMatrix*, dimension (LDA,N)   
            The original M x N sparse matrix A.   

    X       (input) FLOAT PRECISION array, dimension (LDX,NRHS)   
            The computed solution vectors for the system of linear   
            equations.   

    LDX     (input) INTEGER   
            The leading dimension of the array X.  If TRANS = NOTRANS,   
            LDX >= max(1,N); if TRANS = TRANS or CONJ, LDX >= max(1,M).   

    B       (input/output) FLOAT PRECISION array, dimension (LDB,NRHS)   
            On entry, the right hand side vectors for the system of   
            linear equations.   
            On exit, B is overwritten with the difference B - A*X.   

    LDB     (input) INTEGER   
            The leading dimension of the array B.  IF TRANS = NOTRANS,
            LDB >= max(1,M); if TRANS = TRANS or CONJ, LDB >= max(1,N).
	    
    RESID   (output) FLOAT PRECISION   
            The maximum over the number of right hand sides of   
            norm(B - A*X) / ( norm(A) * norm(X) * EPS ).   

    =====================================================================
*/

    /* Table of constant values */
    float alpha = -1.;
    float beta  = 1.;
    int    c__1  = 1;
    
    /* System generated locals */
    float d__1, d__2;

    /* Local variables */
    int j;
    int n1, n2;
    float anorm, bnorm;
    float xnorm;
    float eps;
    char transc[1];

    /* Function prototypes */
    extern float slangs(char *, SuperMatrix *);
    extern float sasum_(int *, float *, int *);
    
    /* Function Body */
    if ( m <= 0 || n <= 0 || nrhs == 0) {
	*resid = 0.;
	return 0;
    }

    if ( (trans == TRANS) || (trans == CONJ) ) {
	n1 = n;
	n2 = m;
        *transc = 'T';
    } else {
	n1 = m;
	n2 = n;
	*transc = 'N';
    }

    /* Exit with RESID = 1/EPS if ANORM = 0. */
    eps = smach("Epsilon");
    anorm = slangs("1", A);
    if (anorm <= 0.) {
	*resid = 1. / eps;
	return 0;
    }

    /* Compute  B - A*X  (or  B - A'*X ) and store in B. */

    sp_sgemm(transc, "N", n1, nrhs, n2, alpha, A, x, ldx, beta, b, ldb);

    /* Compute the maximum over the number of right hand sides of   
       norm(B - A*X) / ( norm(A) * norm(X) * EPS ) . */

    *resid = 0.;
    for (j = 0; j < nrhs; ++j) {
	bnorm = sasum_(&n1, &b[j*ldb], &c__1);
	xnorm = sasum_(&n2, &x[j*ldx], &c__1);
	if (xnorm <= 0.) {
	    *resid = 1. / eps;
	} else {
	    /* Computing MAX */
	    d__1 = *resid, d__2 = bnorm / anorm / xnorm / eps;
	    *resid = SUPERLU_MAX(d__1, d__2);
	}
    }

    return 0;

} /* sgst02 */

