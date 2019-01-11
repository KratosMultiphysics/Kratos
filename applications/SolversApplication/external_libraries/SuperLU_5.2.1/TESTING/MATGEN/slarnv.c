#include "f2c.h"

/* Subroutine */ int slarnv_slu(integer *idist, integer *iseed, integer *n, real 
	*x)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    SLARNV returns a vector of n random real numbers from a uniform or   
    normal distribution.   

    Arguments   
    =========   

    IDIST   (input) INTEGER   
            Specifies the distribution of the random numbers:   
            = 1:  uniform (0,1)   
            = 2:  uniform (-1,1)   
            = 3:  normal (0,1)   

    ISEED   (input/output) INTEGER array, dimension (4)   
            On entry, the seed of the random number generator; the array 
  
            elements must be between 0 and 4095, and ISEED(4) must be   
            odd.   
            On exit, the seed is updated.   

    N       (input) INTEGER   
            The number of random numbers to be generated.   

    X       (output) REAL array, dimension (N)   
            The generated random numbers.   

    Further Details   
    ===============   

    This routine calls the auxiliary routine SLARUV to generate random   
    real numbers from a uniform (0,1) distribution, in batches of up to   
    128 using vectorisable code. The Box-Muller method is used to   
    transform numbers from a uniform to a normal distribution.   

    ===================================================================== 
  


    
   Parameter adjustments   
       Function Body */
    /* System generated locals */
    integer i__1, i__2, i__3;
    /* Builtin functions */
    double log(doublereal), sqrt(doublereal), cos(doublereal);
    /* Local variables */
    static integer i;
    static real u[128];
    static integer il, iv, il2;
    extern /* Subroutine */ int slaruv_slu(integer *, integer *, real *);


#define X(I) x[(I)-1]
#define ISEED(I) iseed[(I)-1]


    i__1 = *n;
    for (iv = 1; iv <= *n; iv += 64) {
/* Computing MIN */
	i__2 = 64, i__3 = *n - iv + 1;
	il = min(i__2,i__3);
	if (*idist == 3) {
	    il2 = il << 1;
	} else {
	    il2 = il;
	}

/*        Call SLARUV to generate IL2 numbers from a uniform (0,1)   
          distribution (IL2 <= LV) */

	slaruv_slu(&ISEED(1), &il2, u);

	if (*idist == 1) {

/*           Copy generated numbers */

	    i__2 = il;
	    for (i = 1; i <= il; ++i) {
		X(iv + i - 1) = u[i - 1];
/* L10: */
	    }
	} else if (*idist == 2) {

/*           Convert generated numbers to uniform (-1,1) distribut
ion */

	    i__2 = il;
	    for (i = 1; i <= il; ++i) {
		X(iv + i - 1) = u[i - 1] * 2.f - 1.f;
/* L20: */
	    }
	} else if (*idist == 3) {

/*           Convert generated numbers to normal (0,1) distributio
n */

	    i__2 = il;
	    for (i = 1; i <= il; ++i) {
		X(iv + i - 1) = sqrt(log(u[(i << 1) - 2]) * -2.f) * cos(u[(i 
			<< 1) - 1] * 6.2831853071795864769252867663f);
/* L30: */
	    }
	}
/* L40: */
    }
    return 0;

/*     End of SLARNV */

} /* slarnv_slu */

