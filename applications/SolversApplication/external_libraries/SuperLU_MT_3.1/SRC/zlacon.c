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
#include "slu_mt_zdefs.h"
#include "slu_mt_util.h"
#include "slu_dcomplex.h"

int_t
zlacon_(int_t *n, doublecomplex *v, doublecomplex *x, double *est, int_t *kase)

{
/*
    Purpose   
    =======   

    ZLACON estimates the 1-norm of a square matrix A.   
    Reverse communication is used for evaluating matrix-vector products. 
  

    Arguments   
    =========   

    N      (input) INT_T
           The order of the matrix.  N >= 1.   

    V      (workspace) DOUBLE COMPLEX PRECISION array, dimension (N)   
           On the final return, V = A*W,  where  EST = norm(V)/norm(W)   
           (W is not returned).   

    X      (input/output) DOUBLE COMPLEX PRECISION array, dimension (N)   
           On an intermediate return, X should be overwritten by   
                 A * X,   if KASE=1,   
                 A' * X,  if KASE=2,
           where A' is the conjugate transpose of A,
           and ZLACON must be re-called with all the other parameters   
           unchanged.   


    EST    (output) DOUBLE PRECISION   
           An estimate (a lower bound) for norm(A).   

    KASE   (input/output) INT_T
           On the initial call to ZLACON, KASE should be 0.   
           On an intermediate return, KASE will be 1 or 2, indicating   
           whether X should be overwritten by A * X  or A' * X.   
           On the final return from ZLACON, KASE will again be 0.   

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
    doublecomplex      zero = {0.0, 0.0};
    doublecomplex      one = {1.0, 0.0};

    /* System generated locals */
    double d__1;
    
    /* Local variables */
    static int_t iter;
    static int_t jump, jlast;
    static double altsgn, estold;
    static int_t i, j;
    double temp;
    double safmin;
    extern double dlamch_(char *);
    extern int izmax1_(int *, doublecomplex *, int *);
    extern double dzsum1_(int *, doublecomplex *, int *);

    safmin = dlamch_("Safe minimum");
    if ( *kase == 0 ) {
	for (i = 0; i < *n; ++i) {
	    x[i].r = 1. / (double) (*n);
	    x[i].i = 0.;
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
	*est = z_abs(&v[0]);
	/*        ... QUIT */
	goto L150;
    }
    *est = dzsum1_(&ni, x, &c__1);

    for (i = 0; i < *n; ++i) {
	d__1 = z_abs(&x[i]);
	if (d__1 > safmin) {
	    d__1 = 1 / d__1;
	    x[i].r *= d__1;
	    x[i].i *= d__1;
	} else {
	    x[i] = one;
	}
    }
    *kase = 2;
    jump = 2;
    return 0;

    /*     ................ ENTRY   (JUMP = 2)   
	   FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY TRANSPOSE(A)*X. */
L40:
    j = izmax1_(&ni, &x[0], &c__1);
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
    zcopy_(&ni, &x[0], &c__1, &v[0], &c__1);
    estold = *est;
    *est = dzsum1_(&ni, v, &c__1);


L90:
    /*     TEST FOR CYCLING. */
    if (*est <= estold) goto L120;

    for (i = 0; i < *n; ++i) {
	d__1 = z_abs(&x[i]);
	if (d__1 > safmin) {
	    d__1 = 1 / d__1;
	    x[i].r *= d__1;
	    x[i].i *= d__1;
	} else {
	    x[i] = one;
	}
    }
    *kase = 2;
    jump = 4;
    return 0;

    /*     ................ ENTRY   (JUMP = 4)   
	   X HAS BEEN OVERWRITTEN BY TRANDPOSE(A)*X. */
L110:
    jlast = j;
    j = izmax1_(&ni, &x[0], &c__1);
    --j;
    if (x[jlast].r != (d__1 = x[j].r, fabs(d__1)) && iter < 5) {
	++iter;
	goto L50;
    }

    /*     ITERATION COMPLETE.  FINAL STAGE. */
L120:
    altsgn = 1.;
    for (i = 1; i <= *n; ++i) {
	x[i-1].r = altsgn * ((double) (i - 1) / (double) (*n - 1) + 1.);
	x[i-1].i = 0.;
	altsgn = -altsgn;
    }
    *kase = 1;
    jump = 5;
    return 0;
    
    /*     ................ ENTRY   (JUMP = 5)   
	   X HAS BEEN OVERWRITTEN BY A*X. */
L140:
    temp = dzsum1_(&ni, x, &c__1) / (double) (*n * 3) * 2.;
    if (temp > *est) {
	zcopy_(&ni, &x[0], &c__1, &v[0], &c__1);
	*est = temp;
    }

L150:
    *kase = 0;
    return 0;

} /* zlacon_ */
