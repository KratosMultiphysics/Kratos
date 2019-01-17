/*  -- translated by f2c (version 19940927).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__3 = 3;
static integer c__1 = 1;
static doublereal c_b11 = 1.;
static doublereal c_b13 = 0.;

/* Subroutine */ int dlagge_slu(integer *m, integer *n, integer *kl, integer *ku,
	 doublereal *d, doublereal *a, integer *lda, integer *iseed, 
	doublereal *work, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    doublereal d__1;

    /* Builtin functions */
    double d_sign(doublereal *, doublereal *);

    /* Local variables */
    extern /* Subroutine */ int dger_(integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    static integer i, j;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *), dgemv_(char *, integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *);
    static doublereal wa, wb, wn;
    extern /* Subroutine */ int dlarnv_slu(integer *, integer *, integer *, doublereal *);
    extern int input_error(char *, int *);
    static doublereal tau;


/*  -- LAPACK auxiliary test routine (version 2.0)   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       February 29, 1992   


    Purpose   
    =======   

    DLAGGE generates a real general m by n matrix A, by pre- and post-   
    multiplying a real diagonal matrix D with random orthogonal matrices: 
  
    A = U*D*V. The lower and upper bandwidths may then be reduced to   
    kl and ku by additional orthogonal transformations.   

    Arguments   
    =========   

    M       (input) INTEGER   
            The number of rows of the matrix A.  M >= 0.   

    N       (input) INTEGER   
            The number of columns of the matrix A.  N >= 0.   

    KL      (input) INTEGER   
            The number of nonzero subdiagonals within the band of A.   
            0 <= KL <= M-1.   

    KU      (input) INTEGER   
            The number of nonzero superdiagonals within the band of A.   
            0 <= KU <= N-1.   

    D       (input) DOUBLE PRECISION array, dimension (min(M,N))   
            The diagonal elements of the diagonal matrix D.   

    A       (output) DOUBLE PRECISION array, dimension (LDA,N)   
            The generated m by n matrix A.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= M.   

    ISEED   (input/output) INTEGER array, dimension (4)   
            On entry, the seed of the random number generator; the array 
  
            elements must be between 0 and 4095, and ISEED(4) must be   
            odd.   
            On exit, the seed is updated.   

    WORK    (workspace) DOUBLE PRECISION array, dimension (M+N)   

    INFO    (output) INTEGER   
            = 0: successful exit   
            < 0: if INFO = -i, the i-th argument had an illegal value   

    ===================================================================== 
  


       Test the input arguments   

       Parameter adjustments */
    --d;
    a_dim1 = *lda;
    a_offset = a_dim1 + 1;
    a -= a_offset;
    --iseed;
    --work;

    /* Function Body */
    *info = 0;
    if (*m < 0) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*kl < 0 || *kl > *m - 1) {
	*info = -3;
    } else if (*ku < 0 || *ku > *n - 1) {
	*info = -4;
    } else if (*lda < max(1,*m)) {
	*info = -7;
    }
    if (*info < 0) {
	i__1 = -(*info);
	input_error("DLAGGE", &i__1);
	return 0;
    }

/*     initialize A to diagonal matrix */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *m;
	for (i = 1; i <= i__2; ++i) {
	    a[i + j * a_dim1] = 0.;
/* L10: */
	}
/* L20: */
    }
    i__1 = min(*m,*n);
    for (i = 1; i <= i__1; ++i) {
	a[i + i * a_dim1] = d[i];
/* L30: */
    }

/*     pre- and post-multiply A by random orthogonal matrices */

    for (i = min(*m,*n); i >= 1; --i) {
	if (i < *m) {

/*           generate random reflection */

	    i__1 = *m - i + 1;
	    dlarnv_slu(&c__3, &iseed[1], &i__1, &work[1]);
	    i__1 = *m - i + 1;
	    wn = dnrm2_(&i__1, &work[1], &c__1);
	    wa = d_sign(&wn, &work[1]);
	    if (wn == 0.) {
		tau = 0.;
	    } else {
		wb = work[1] + wa;
		i__1 = *m - i;
		d__1 = 1. / wb;
		dscal_(&i__1, &d__1, &work[2], &c__1);
		work[1] = 1.;
		tau = wb / wa;
	    }

/*           multiply A(i:m,i:n) by random reflection from the lef
t */

	    i__1 = *m - i + 1;
	    i__2 = *n - i + 1;
	    dgemv_("Transpose", &i__1, &i__2, &c_b11, &a[i + i * a_dim1], lda,
		     &work[1], &c__1, &c_b13, &work[*m + 1], &c__1);
	    i__1 = *m - i + 1;
	    i__2 = *n - i + 1;
	    d__1 = -tau;
	    dger_(&i__1, &i__2, &d__1, &work[1], &c__1, &work[*m + 1], &c__1, 
		    &a[i + i * a_dim1], lda);
	}
	if (i < *n) {

/*           generate random reflection */

	    i__1 = *n - i + 1;
	    dlarnv_slu(&c__3, &iseed[1], &i__1, &work[1]);
	    i__1 = *n - i + 1;
	    wn = dnrm2_(&i__1, &work[1], &c__1);
	    wa = d_sign(&wn, &work[1]);
	    if (wn == 0.) {
		tau = 0.;
	    } else {
		wb = work[1] + wa;
		i__1 = *n - i;
		d__1 = 1. / wb;
		dscal_(&i__1, &d__1, &work[2], &c__1);
		work[1] = 1.;
		tau = wb / wa;
	    }

/*           multiply A(i:m,i:n) by random reflection from the rig
ht */

	    i__1 = *m - i + 1;
	    i__2 = *n - i + 1;
	    dgemv_("No transpose", &i__1, &i__2, &c_b11, &a[i + i * a_dim1], 
		    lda, &work[1], &c__1, &c_b13, &work[*n + 1], &c__1);
	    i__1 = *m - i + 1;
	    i__2 = *n - i + 1;
	    d__1 = -tau;
	    dger_(&i__1, &i__2, &d__1, &work[*n + 1], &c__1, &work[1], &c__1, 
		    &a[i + i * a_dim1], lda);
	}
/* L40: */
    }

/*     Reduce number of subdiagonals to KL and number of superdiagonals   
       to KU   

   Computing MAX */
    i__2 = *m - 1 - *kl, i__3 = *n - 1 - *ku;
    i__1 = max(i__2,i__3);
    for (i = 1; i <= i__1; ++i) {
	if (*kl <= *ku) {

/*           annihilate subdiagonal elements first (necessary if K
L = 0)   

   Computing MIN */
	    i__2 = *m - 1 - *kl;
	    if (i <= min(i__2,*n)) {

/*              generate reflection to annihilate A(kl+i+1:m,i
) */

		i__2 = *m - *kl - i + 1;
		wn = dnrm2_(&i__2, &a[*kl + i + i * a_dim1], &c__1);
		wa = d_sign(&wn, &a[*kl + i + i * a_dim1]);
		if (wn == 0.) {
		    tau = 0.;
		} else {
		    wb = a[*kl + i + i * a_dim1] + wa;
		    i__2 = *m - *kl - i;
		    d__1 = 1. / wb;
		    dscal_(&i__2, &d__1, &a[*kl + i + 1 + i * a_dim1], &c__1);
		    a[*kl + i + i * a_dim1] = 1.;
		    tau = wb / wa;
		}

/*              apply reflection to A(kl+i:m,i+1:n) from the l
eft */

		i__2 = *m - *kl - i + 1;
		i__3 = *n - i;
		dgemv_("Transpose", &i__2, &i__3, &c_b11, &a[*kl + i + (i + 1)
			 * a_dim1], lda, &a[*kl + i + i * a_dim1], &c__1, &
			c_b13, &work[1], &c__1);
		i__2 = *m - *kl - i + 1;
		i__3 = *n - i;
		d__1 = -tau;
		dger_(&i__2, &i__3, &d__1, &a[*kl + i + i * a_dim1], &c__1, &
			work[1], &c__1, &a[*kl + i + (i + 1) * a_dim1], lda);
		a[*kl + i + i * a_dim1] = -wa;
	    }

/* Computing MIN */
	    i__2 = *n - 1 - *ku;
	    if (i <= min(i__2,*m)) {

/*              generate reflection to annihilate A(i,ku+i+1:n
) */

		i__2 = *n - *ku - i + 1;
		wn = dnrm2_(&i__2, &a[i + (*ku + i) * a_dim1], lda);
		wa = d_sign(&wn, &a[i + (*ku + i) * a_dim1]);
		if (wn == 0.) {
		    tau = 0.;
		} else {
		    wb = a[i + (*ku + i) * a_dim1] + wa;
		    i__2 = *n - *ku - i;
		    d__1 = 1. / wb;
		    dscal_(&i__2, &d__1, &a[i + (*ku + i + 1) * a_dim1], lda);
		    a[i + (*ku + i) * a_dim1] = 1.;
		    tau = wb / wa;
		}

/*              apply reflection to A(i+1:m,ku+i:n) from the r
ight */

		i__2 = *m - i;
		i__3 = *n - *ku - i + 1;
		dgemv_("No transpose", &i__2, &i__3, &c_b11, &a[i + 1 + (*ku 
			+ i) * a_dim1], lda, &a[i + (*ku + i) * a_dim1], lda, 
			&c_b13, &work[1], &c__1);
		i__2 = *m - i;
		i__3 = *n - *ku - i + 1;
		d__1 = -tau;
		dger_(&i__2, &i__3, &d__1, &work[1], &c__1, &a[i + (*ku + i) *
			 a_dim1], lda, &a[i + 1 + (*ku + i) * a_dim1], lda);
		a[i + (*ku + i) * a_dim1] = -wa;
	    }
	} else {

/*           annihilate superdiagonal elements first (necessary if
   
             KU = 0)   

   Computing MIN */
	    i__2 = *n - 1 - *ku;
	    if (i <= min(i__2,*m)) {

/*              generate reflection to annihilate A(i,ku+i+1:n
) */

		i__2 = *n - *ku - i + 1;
		wn = dnrm2_(&i__2, &a[i + (*ku + i) * a_dim1], lda);
		wa = d_sign(&wn, &a[i + (*ku + i) * a_dim1]);
		if (wn == 0.) {
		    tau = 0.;
		} else {
		    wb = a[i + (*ku + i) * a_dim1] + wa;
		    i__2 = *n - *ku - i;
		    d__1 = 1. / wb;
		    dscal_(&i__2, &d__1, &a[i + (*ku + i + 1) * a_dim1], lda);
		    a[i + (*ku + i) * a_dim1] = 1.;
		    tau = wb / wa;
		}

/*              apply reflection to A(i+1:m,ku+i:n) from the r
ight */

		i__2 = *m - i;
		i__3 = *n - *ku - i + 1;
		dgemv_("No transpose", &i__2, &i__3, &c_b11, &a[i + 1 + (*ku 
			+ i) * a_dim1], lda, &a[i + (*ku + i) * a_dim1], lda, 
			&c_b13, &work[1], &c__1);
		i__2 = *m - i;
		i__3 = *n - *ku - i + 1;
		d__1 = -tau;
		dger_(&i__2, &i__3, &d__1, &work[1], &c__1, &a[i + (*ku + i) *
			 a_dim1], lda, &a[i + 1 + (*ku + i) * a_dim1], lda);
		a[i + (*ku + i) * a_dim1] = -wa;
	    }

/* Computing MIN */
	    i__2 = *m - 1 - *kl;
	    if (i <= min(i__2,*n)) {

/*              generate reflection to annihilate A(kl+i+1:m,i
) */

		i__2 = *m - *kl - i + 1;
		wn = dnrm2_(&i__2, &a[*kl + i + i * a_dim1], &c__1);
		wa = d_sign(&wn, &a[*kl + i + i * a_dim1]);
		if (wn == 0.) {
		    tau = 0.;
		} else {
		    wb = a[*kl + i + i * a_dim1] + wa;
		    i__2 = *m - *kl - i;
		    d__1 = 1. / wb;
		    dscal_(&i__2, &d__1, &a[*kl + i + 1 + i * a_dim1], &c__1);
		    a[*kl + i + i * a_dim1] = 1.;
		    tau = wb / wa;
		}

/*              apply reflection to A(kl+i:m,i+1:n) from the l
eft */

		i__2 = *m - *kl - i + 1;
		i__3 = *n - i;
		dgemv_("Transpose", &i__2, &i__3, &c_b11, &a[*kl + i + (i + 1)
			 * a_dim1], lda, &a[*kl + i + i * a_dim1], &c__1, &
			c_b13, &work[1], &c__1);
		i__2 = *m - *kl - i + 1;
		i__3 = *n - i;
		d__1 = -tau;
		dger_(&i__2, &i__3, &d__1, &a[*kl + i + i * a_dim1], &c__1, &
			work[1], &c__1, &a[*kl + i + (i + 1) * a_dim1], lda);
		a[*kl + i + i * a_dim1] = -wa;
	    }
	}

	i__2 = *m;
	for (j = *kl + i + 1; j <= i__2; ++j) {
	    a[j + i * a_dim1] = 0.;
/* L50: */
	}

	i__2 = *n;
	for (j = *ku + i + 1; j <= i__2; ++j) {
	    a[i + j * a_dim1] = 0.;
/* L60: */
	}
/* L70: */
    }
    return 0;

/*     End of DLAGGE */

} /* dlagge_slu */

