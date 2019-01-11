/*  -- translated by f2c (version 19940927).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include <string.h>
#include "f2c.h"

/* Table of constant values */

static integer c__2 = 2;

/* Subroutine */ int clatb4_slu(char *path, integer *imat, integer *m, integer *
	n, char *type, integer *kl, integer *ku, real *anorm, integer *mode, 
	real *cndnum, char *dist)
{
    /* Initialized data */

    static logical first = TRUE_;

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static real badc1, badc2, large, small;
    static char c2[2];
    extern /* Subroutine */ int slabad_slu(real *, real *);
    extern float smach(char *);
    static integer mat;
    static real eps;


/*  -- LAPACK test routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       February 29, 1992   


    Purpose   
    =======   

    CLATB4 sets parameters for the matrix generator based on the type of 
  
    matrix to be generated.   

    Arguments   
    =========   

    PATH    (input) CHARACTER*3   
            The LAPACK path name.   

    IMAT    (input) INTEGER   
            An integer key describing which matrix to generate for this   
            path.   

    M       (input) INTEGER   
            The number of rows in the matrix to be generated.   

    N       (input) INTEGER   
            The number of columns in the matrix to be generated.   

    TYPE    (output) CHARACTER*1   
            The type of the matrix to be generated:   
            = 'S':  symmetric matrix   
            = 'P':  symmetric positive (semi)definite matrix   
            = 'N':  nonsymmetric matrix   

    KL      (output) INTEGER   
            The lower band width of the matrix to be generated.   

    KU      (output) INTEGER   
            The upper band width of the matrix to be generated.   

    ANORM   (output) REAL   
            The desired norm of the matrix to be generated.  The diagonal 
  
            matrix of singular values or eigenvalues is scaled by this   
            value.   

    MODE    (output) INTEGER   
            A key indicating how to choose the vector of eigenvalues.   

    CNDNUM  (output) REAL   
            The desired condition number.   

    DIST    (output) CHARACTER*1   
            The type of distribution to be used by the random number   
            generator.   

    ===================================================================== 
  


       Set some constants for use in the subroutine. */

    if (first) {
	first = FALSE_;
	eps = smach("Precision");
	badc2 = .1f / eps;
	badc1 = sqrt(badc2);
	small = smach("Safe minimum");
	large = 1.f / small;

/*        If it looks like we're on a Cray, take the square root of   
          SMALL and LARGE to avoid overflow and underflow problems. */

	slabad_slu(&small, &large);
	small = small / eps * .25f;
	large = 1.f / small;
    }

    /*    s_copy(c2, path + 1, 2L, 2L);*/
    strncpy(c2, path + 1, 2);

/*     Set some parameters we don't plan to change. */

    *(unsigned char *)dist = 'S';
    *mode = 3;

/*     xQR, xLQ, xQL, xRQ:  Set parameters to generate a general   
                            M x N matrix. */

    if (strncmp(c2, "QR", 2)==0 || strncmp(c2, "LQ", 2)==0  
	|| strncmp(c2, "QL", 2)==0 || strncmp(c2, "RQ", 2)==0) {

/*        Set TYPE, the type of matrix to be generated. */

	*(unsigned char *)type = 'N';

/*        Set the lower and upper bandwidths. */

	if (*imat == 1) {
	    *kl = 0;
	    *ku = 0;
	} else if (*imat == 2) {
	    *kl = 0;
/* Computing MAX */
	    i__1 = *n - 1;
	    *ku = max(i__1,0);
	} else if (*imat == 3) {
/* Computing MAX */
	    i__1 = *m - 1;
	    *kl = max(i__1,0);
	    *ku = 0;
	} else {
/* Computing MAX */
	    i__1 = *m - 1;
	    *kl = max(i__1,0);
/* Computing MAX */
	    i__1 = *n - 1;
	    *ku = max(i__1,0);
	}

/*        Set the condition number and norm. */

	if (*imat == 5) {
	    *cndnum = badc1;
	} else if (*imat == 6) {
	    *cndnum = badc2;
	} else {
	    *cndnum = 2.f;
	}

	if (*imat == 7) {
	    *anorm = small;
	} else if (*imat == 8) {
	    *anorm = large;
	} else {
	    *anorm = 1.f;
	}

    } else if (strncmp(c2, "GE", 2)==0) {

/*        xGE:  Set parameters to generate a general M x N matrix.   

          Set TYPE, the type of matrix to be generated. */

	*(unsigned char *)type = 'N';

/*        Set the lower and upper bandwidths. */

	if (*imat == 1) {
	    *kl = 0;
	    *ku = 0;
	} else if (*imat == 2) {
	    *kl = 0;
/* Computing MAX */
	    i__1 = *n - 1;
	    *ku = max(i__1,0);
	} else if (*imat == 3) {
/* Computing MAX */
	    i__1 = *m - 1;
	    *kl = max(i__1,0);
	    *ku = 0;
	} else {
/* Computing MAX */
	    i__1 = *m - 1;
	    *kl = max(i__1,0);
/* Computing MAX */
	    i__1 = *n - 1;
	    *ku = max(i__1,0);
	}

/*        Set the condition number and norm. */

	if (*imat == 8) {
	    *cndnum = badc1;
	} else if (*imat == 9) {
	    *cndnum = badc2;
	} else {
	    *cndnum = 2.f;
	}

	if (*imat == 10) {
	    *anorm = small;
	} else if (*imat == 11) {
	    *anorm = large;
	} else {
	    *anorm = 1.f;
	}

    } else if (strncmp(c2, "GB", 2)==0) {

/*        xGB:  Set parameters to generate a general banded matrix.   

          Set TYPE, the type of matrix to be generated. */

	*(unsigned char *)type = 'N';

/*        Set the condition number and norm. */

	if (*imat == 5) {
	    *cndnum = badc1;
	} else if (*imat == 6) {
	    *cndnum = badc2 * .1f;
	} else {
	    *cndnum = 2.f;
	}

	if (*imat == 7) {
	    *anorm = small;
	} else if (*imat == 8) {
	    *anorm = large;
	} else {
	    *anorm = 1.f;
	}

    } else if (strncmp(c2, "GT", 2)==0) {

/*        xGT:  Set parameters to generate a general tridiagonal matri
x.   

          Set TYPE, the type of matrix to be generated. */

	*(unsigned char *)type = 'N';

/*        Set the lower and upper bandwidths. */

	if (*imat == 1) {
	    *kl = 0;
	} else {
	    *kl = 1;
	}
	*ku = *kl;

/*        Set the condition number and norm. */

	if (*imat == 3) {
	    *cndnum = badc1;
	} else if (*imat == 4) {
	    *cndnum = badc2;
	} else {
	    *cndnum = 2.f;
	}

	if (*imat == 5 || *imat == 11) {
	    *anorm = small;
	} else if (*imat == 6 || *imat == 12) {
	    *anorm = large;
	} else {
	    *anorm = 1.f;
	}

    } else if (strncmp(c2, "PO", 2)==0 || strncmp(c2, "PP", 2)==0 || 
	       strncmp(c2, "HE", 2)==0 || strncmp(c2, "HP", 2)==0 ||
	       strncmp(c2, "SY", 2)==0 || strncmp(c2, "SP", 2)==0) {

/*        xPO, xPP, xHE, xHP, xSY, xSP: Set parameters to generate a 
  
          symmetric or Hermitian matrix.   

          Set TYPE, the type of matrix to be generated. */

	*(unsigned char *)type = *(unsigned char *)c2;

/*        Set the lower and upper bandwidths. */

	if (*imat == 1) {
	    *kl = 0;
	} else {
/* Computing MAX */
	    i__1 = *n - 1;
	    *kl = max(i__1,0);
	}
	*ku = *kl;

/*        Set the condition number and norm. */

	if (*imat == 6) {
	    *cndnum = badc1;
	} else if (*imat == 7) {
	    *cndnum = badc2;
	} else {
	    *cndnum = 2.f;
	}

	if (*imat == 8) {
	    *anorm = small;
	} else if (*imat == 9) {
	    *anorm = large;
	} else {
	    *anorm = 1.f;
	}

    } else if (strncmp(c2, "PB", 2)==0) {

/*        xPB:  Set parameters to generate a symmetric band matrix.   

          Set TYPE, the type of matrix to be generated. */

	*(unsigned char *)type = 'P';

/*        Set the norm and condition number. */

	if (*imat == 5) {
	    *cndnum = badc1;
	} else if (*imat == 6) {
	    *cndnum = badc2;
	} else {
	    *cndnum = 2.f;
	}

	if (*imat == 7) {
	    *anorm = small;
	} else if (*imat == 8) {
	    *anorm = large;
	} else {
	    *anorm = 1.f;
	}

    } else if (strncmp(c2, "PT", 2)==0) {

/*        xPT:  Set parameters to generate a symmetric positive defini
te   
          tridiagonal matrix. */

	*(unsigned char *)type = 'P';
	if (*imat == 1) {
	    *kl = 0;
	} else {
	    *kl = 1;
	}
	*ku = *kl;

/*        Set the condition number and norm. */

	if (*imat == 3) {
	    *cndnum = badc1;
	} else if (*imat == 4) {
	    *cndnum = badc2;
	} else {
	    *cndnum = 2.f;
	}

	if (*imat == 5 || *imat == 11) {
	    *anorm = small;
	} else if (*imat == 6 || *imat == 12) {
	    *anorm = large;
	} else {
	    *anorm = 1.f;
	}

    } else if (strncmp(c2, "TR", 2)==0 || strncmp(c2, "TP", 2)==0) {

/*        xTR, xTP:  Set parameters to generate a triangular matrix   

          Set TYPE, the type of matrix to be generated. */

	*(unsigned char *)type = 'N';

/*        Set the lower and upper bandwidths. */

	mat = abs(*imat);
	if (mat == 1 || mat == 7) {
	    *kl = 0;
	    *ku = 0;
	} else if (*imat < 0) {
/* Computing MAX */
	    i__1 = *n - 1;
	    *kl = max(i__1,0);
	    *ku = 0;
	} else {
	    *kl = 0;
/* Computing MAX */
	    i__1 = *n - 1;
	    *ku = max(i__1,0);
	}

/*        Set the condition number and norm. */

	if (mat == 3 || mat == 9) {
	    *cndnum = badc1;
	} else if (mat == 4 || mat == 10) {
	    *cndnum = badc2;
	} else {
	    *cndnum = 2.f;
	}

	if (mat == 5) {
	    *anorm = small;
	} else if (mat == 6) {
	    *anorm = large;
	} else {
	    *anorm = 1.f;
	}

    } else if (strncmp(c2, "TB", 2)==0) {

/*        xTB:  Set parameters to generate a triangular band matrix. 
  

          Set TYPE, the type of matrix to be generated. */

	*(unsigned char *)type = 'N';

/*        Set the norm and condition number. */

	if (*imat == 2 || *imat == 8) {
	    *cndnum = badc1;
	} else if (*imat == 3 || *imat == 9) {
	    *cndnum = badc2;
	} else {
	    *cndnum = 2.f;
	}

	if (*imat == 4) {
	    *anorm = small;
	} else if (*imat == 5) {
	    *anorm = large;
	} else {
	    *anorm = 1.f;
	}
    }
    if (*n <= 1) {
	*cndnum = 1.f;
    }

    return 0;

/*     End of CLATB4 */

} /* clatb4_slu */

