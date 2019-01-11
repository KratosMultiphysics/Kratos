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
 * Lawrence Berkeley National Lab, Univ. of California Berkeley,
 * and Xerox Palo Alto Research Center.
 * September 10, 2007
 *
 */
#include <math.h>
#include "slu_mt_Cnames.h"
#include "slu_mt_ddefs.h"
#include "slu_mt_util.h"

int_t
dlacon_(int_t *n, double *v, double *x, int_t *isgn, double *est, int_t *kase)

{
/*
    Purpose   
    =======   

    DLACON estimates the 1-norm of a square matrix A.   
    Reverse communication is used for evaluating matrix-vector products. 
  

    Arguments   
    =========   

    N      (input) INT_T
           The order of the matrix.  N >= 1.   

    V      (workspace) DOUBLE PRECISION array, dimension (N)   
           On the final return, V = A*W,  where  EST = norm(V)/norm(W)   
           (W is not returned).   

    X      (input/output) DOUBLE PRECISION array, dimension (N)   
           On an intermediate return, X should be overwritten by   
                 A * X,   if KASE=1,   
                 A' * X,  if KASE=2,
           and DLACON must be re-called with all the other parameters   
           unchanged.   

    ISGN   (workspace) INT_T array, dimension (N)

    EST    (output) DOUBLE PRECISION   
           An estimate (a lower bound) for norm(A).   

    KASE   (input/output) INT_T
           On the initial call to DLACON, KASE should be 0.   
           On an intermediate return, KASE will be 1 or 2, indicating   
           whether X should be overwritten by A * X  or A' * X.   
           On the final return from DLACON, KASE will again be 0.   

    Further Details   
    ======= =======   

    Contributed by Nick Higham, University of Manchester.   
    Originally named CONEST, dated March 16, 1988.   

    Reference: N.J. Higham, "FORTRAN codes for estimating the one-norm of 
    a real or complex matrix, with applications to condition estimation", 
    ACM Trans. Math. Soft., vol. 14, no. 4, pp. 381-396, December 1988.   
    ===================================================================== 
*/

    /* Table of constant values */
    int c__1 = 1;
    int ni = *n;
    double      zero = 0.0;
    double      one = 1.0;
    
    /* Local variables */
    static int_t iter;
    static int_t jump, jlast;
    static double altsgn, estold;
    static int_t i, j;
    double temp;
    extern int idamax_(int *, double *, int *);
    extern double dasum_(int *, double *, int *);
    extern int dcopy_(int *, double *, int *, double *, int *);
#define d_sign(a, b) (b >= 0 ? fabs(a) : -fabs(a))    /* Copy sign */
#define i_dnnt(a) \
	( a>=0 ? floor(a+.5) : -floor(.5-a) ) /* Round to nearest integer */

    if ( *kase == 0 ) {
	for (i = 0; i < *n; ++i) {
	    x[i] = 1. / (double) (*n);
	}
	*kase = 1;
	jump = 1;
	return 0;
    }

    switch (jump) {
	case 1:  goto L20;
	case 2:  goto L40;
	case 3:  goto L70;
	case 4:  goto L110;
	case 5:  goto L140;
    }

    /*     ................ ENTRY   (JUMP = 1)   
	   FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY A*X. */
  L20:
    if (*n == 1) {
	v[0] = x[0];
	*est = fabs(v[0]);
	/*        ... QUIT */
	goto L150;
    }
    *est = dasum_(&ni, x, &c__1);

    for (i = 0; i < *n; ++i) {
	x[i] = d_sign(one, x[i]);
	isgn[i] = i_dnnt(x[i]);
    }
    *kase = 2;
    jump = 2;
    return 0;

    /*     ................ ENTRY   (JUMP = 2)   
	   FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY TRANSPOSE(A)*X. */
L40:
    j = idamax_(&ni, &x[0], &c__1);
    --j;
    iter = 2;

    /*     MAIN LOOP - ITERATIONS 2,3,...,ITMAX. */
L50:
    for (i = 0; i < *n; ++i) x[i] = zero;
    x[j] = one;
    *kase = 1;
    jump = 3;
    return 0;

    /*     ................ ENTRY   (JUMP = 3)   
	   X HAS BEEN OVERWRITTEN BY A*X. */
L70:
    dcopy_(&ni, &x[0], &c__1, &v[0], &c__1);
    estold = *est;
    *est = dasum_(&ni, v, &c__1);

    for (i = 0; i < *n; ++i)
	if (i_dnnt(d_sign(one, x[i])) != isgn[i])
	    goto L90;

    /*     REPEATED SIGN VECTOR DETECTED, HENCE ALGORITHM HAS CONVERGED. */
    goto L120;

L90:
    /*     TEST FOR CYCLING. */
    if (*est <= estold) goto L120;

    for (i = 0; i < *n; ++i) {
	x[i] = d_sign(one, x[i]);
	isgn[i] = i_dnnt(x[i]);
    }
    *kase = 2;
    jump = 4;
    return 0;

    /*     ................ ENTRY   (JUMP = 4)   
	   X HAS BEEN OVERWRITTEN BY TRANDPOSE(A)*X. */
L110:
    jlast = j;
    j = idamax_(&ni, &x[0], &c__1);
    --j;
    if (x[jlast] != fabs(x[j]) && iter < 5) {
	++iter;
	goto L50;
    }

    /*     ITERATION COMPLETE.  FINAL STAGE. */
L120:
    altsgn = 1.;
    for (i = 1; i <= *n; ++i) {
	x[i-1] = altsgn * ((double) (i - 1) / (double) (*n - 1) + 1.);
	altsgn = -altsgn;
    }
    *kase = 1;
    jump = 5;
    return 0;
    
    /*     ................ ENTRY   (JUMP = 5)   
	   X HAS BEEN OVERWRITTEN BY A*X. */
L140:
    temp = dasum_(&ni, x, &c__1) / (double) (*n * 3) * 2.;
    if (temp > *est) {
	dcopy_(&ni, &x[0], &c__1, &v[0], &c__1);
	*est = temp;
    }

L150:
    *kase = 0;
    return 0;

} /* dlacon_ */
