
/*  -- translated by f2c (version 19940927).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/
#include <string.h>
#include "f2c.h"

/* Subroutine */ int ctrsv_(char *uplo, char *trans, char *diag, integer *n, 
	complex *a, integer *lda, complex *x, integer *incx)
{


    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5;
    complex q__1, q__2, q__3;

    /* Builtin functions */
    void c_div(complex *, complex *, complex *), r_cnjg(complex *, complex *);

    /* Local variables */
    static integer info;
    static complex temp;
    static integer i, j;
    static integer ix, jx, kx;
    static logical noconj, nounit;

    extern int input_error(char *, int *);

/*  Purpose   
    =======   

    CTRSV  solves one of the systems of equations   

       A*x = b,   or   A'*x = b,   or   conjg( A' )*x = b,   

    where b and x are n element vectors and A is an n by n unit, or   
    non-unit, upper or lower triangular matrix.   

    No test for singularity or near-singularity is included in this   
    routine. Such tests must be performed before calling this routine.   

    Parameters   
    ==========   

    UPLO   - CHARACTER*1.   
             On entry, UPLO specifies whether the matrix is an upper or   
             lower triangular matrix as follows:   

                UPLO = 'U' or 'u'   A is an upper triangular matrix.   

                UPLO = 'L' or 'l'   A is a lower triangular matrix.   

             Unchanged on exit.   

    TRANS  - CHARACTER*1.   
             On entry, TRANS specifies the equations to be solved as   
             follows:   

                TRANS = 'N' or 'n'   A*x = b.   

                TRANS = 'T' or 't'   A'*x = b.   

                TRANS = 'C' or 'c'   conjg( A' )*x = b.   

             Unchanged on exit.   

    DIAG   - CHARACTER*1.   
             On entry, DIAG specifies whether or not A is unit   
             triangular as follows:   

                DIAG = 'U' or 'u'   A is assumed to be unit triangular.   

                DIAG = 'N' or 'n'   A is not assumed to be unit   
                                    triangular.   

             Unchanged on exit.   

    N      - INTEGER.   
             On entry, N specifies the order of the matrix A.   
             N must be at least zero.   
             Unchanged on exit.   

    A      - COMPLEX          array of DIMENSION ( LDA, n ).   
             Before entry with  UPLO = 'U' or 'u', the leading n by n   
             upper triangular part of the array A must contain the upper 
  
             triangular matrix and the strictly lower triangular part of 
  
             A is not referenced.   
             Before entry with UPLO = 'L' or 'l', the leading n by n   
             lower triangular part of the array A must contain the lower 
  
             triangular matrix and the strictly upper triangular part of 
  
             A is not referenced.   
             Note that when  DIAG = 'U' or 'u', the diagonal elements of 
  
             A are not referenced either, but are assumed to be unity.   
             Unchanged on exit.   

    LDA    - INTEGER.   
             On entry, LDA specifies the first dimension of A as declared 
  
             in the calling (sub) program. LDA must be at least   
             max( 1, n ).   
             Unchanged on exit.   

    X      - COMPLEX          array of dimension at least   
             ( 1 + ( n - 1 )*abs( INCX ) ).   
             Before entry, the incremented array X must contain the n   
             element right-hand side vector b. On exit, X is overwritten 
  
             with the solution vector x.   

    INCX   - INTEGER.   
             On entry, INCX specifies the increment for the elements of   
             X. INCX must not be zero.   
             Unchanged on exit.   


    Level 2 Blas routine.   

    -- Written on 22-October-1986.   
       Jack Dongarra, Argonne National Lab.   
       Jeremy Du Croz, Nag Central Office.   
       Sven Hammarling, Nag Central Office.   
       Richard Hanson, Sandia National Labs.   



       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
#define X(I) x[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    info = 0;
    if ( strncmp(uplo, "U", 1)!=0 && strncmp(uplo, "L", 1)!=0 ) {
	info = 1;
    } else if ( strncmp(trans, "N", 1)!=0 && strncmp(trans, "T", 1)!=0 &&
		strncmp(trans, "C", 1)!=0 ) {
	info = 2;
    } else if ( strncmp(diag, "U", 1)!=0 && strncmp(diag, "N", 1)!=0 ) {
	info = 3;
    } else if (*n < 0) {
	info = 4;
    } else if (*lda < max(1,*n)) {
	info = 6;
    } else if (*incx == 0) {
	info = 8;
    }
    if (info != 0) {
	input_error("CTRSV ", &info);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0) {
	return 0;
    }

    noconj = (strncmp(trans, "T", 1)==0);
    nounit = (strncmp(diag, "N", 1)==0);

/*     Set up the start point in X if the increment is not unity. This   
       will be  ( N - 1 )*INCX  too small for descending loops. */

    if (*incx <= 0) {
	kx = 1 - (*n - 1) * *incx;
    } else if (*incx != 1) {
	kx = 1;
    }

/*     Start the operations. In this version the elements of A are   
       accessed sequentially with one pass through A. */

    if (strncmp(trans, "N", 1)==0) {

/*        Form  x := inv( A )*x. */

	if (strncmp(uplo, "U", 1)==0) {
	    if (*incx == 1) {
		for (j = *n; j >= 1; --j) {
		    i__1 = j;
		    if (X(j).r != 0.f || X(j).i != 0.f) {
			if (nounit) {
			    i__1 = j;
			    c_div(&q__1, &X(j), &A(j,j));
			    X(j).r = q__1.r, X(j).i = q__1.i;
			}
			i__1 = j;
			temp.r = X(j).r, temp.i = X(j).i;
			for (i = j - 1; i >= 1; --i) {
			    i__1 = i;
			    i__2 = i;
			    i__3 = i + j * a_dim1;
			    q__2.r = temp.r * A(i,j).r - temp.i * A(i,j).i, 
				    q__2.i = temp.r * A(i,j).i + temp.i * A(i,j).r;
			    q__1.r = X(i).r - q__2.r, q__1.i = X(i).i - 
				    q__2.i;
			    X(i).r = q__1.r, X(i).i = q__1.i;
/* L10: */
			}
		    }
/* L20: */
		}
	    } else {
		jx = kx + (*n - 1) * *incx;
		for (j = *n; j >= 1; --j) {
		    i__1 = jx;
		    if (X(jx).r != 0.f || X(jx).i != 0.f) {
			if (nounit) {
			    i__1 = jx;
			    c_div(&q__1, &X(jx), &A(j,j));
			    X(jx).r = q__1.r, X(jx).i = q__1.i;
			}
			i__1 = jx;
			temp.r = X(jx).r, temp.i = X(jx).i;
			ix = jx;
			for (i = j - 1; i >= 1; --i) {
			    ix -= *incx;
			    i__1 = ix;
			    i__2 = ix;
			    i__3 = i + j * a_dim1;
			    q__2.r = temp.r * A(i,j).r - temp.i * A(i,j).i, 
				    q__2.i = temp.r * A(i,j).i + temp.i * A(i,j).r;
			    q__1.r = X(ix).r - q__2.r, q__1.i = X(ix).i - 
				    q__2.i;
			    X(ix).r = q__1.r, X(ix).i = q__1.i;
/* L30: */
			}
		    }
		    jx -= *incx;
/* L40: */
		}
	    }
	} else {
	    if (*incx == 1) {
		i__1 = *n;
		for (j = 1; j <= *n; ++j) {
		    i__2 = j;
		    if (X(j).r != 0.f || X(j).i != 0.f) {
			if (nounit) {
			    i__2 = j;
			    c_div(&q__1, &X(j), &A(j,j));
			    X(j).r = q__1.r, X(j).i = q__1.i;
			}
			i__2 = j;
			temp.r = X(j).r, temp.i = X(j).i;
			i__2 = *n;
			for (i = j + 1; i <= *n; ++i) {
			    i__3 = i;
			    i__4 = i;
			    i__5 = i + j * a_dim1;
			    q__2.r = temp.r * A(i,j).r - temp.i * A(i,j).i, 
				    q__2.i = temp.r * A(i,j).i + temp.i * A(i,j).r;
			    q__1.r = X(i).r - q__2.r, q__1.i = X(i).i - 
				    q__2.i;
			    X(i).r = q__1.r, X(i).i = q__1.i;
/* L50: */
			}
		    }
/* L60: */
		}
	    } else {
		jx = kx;
		i__1 = *n;
		for (j = 1; j <= *n; ++j) {
		    i__2 = jx;
		    if (X(jx).r != 0.f || X(jx).i != 0.f) {
			if (nounit) {
			    i__2 = jx;
			    c_div(&q__1, &X(jx), &A(j,j));
			    X(jx).r = q__1.r, X(jx).i = q__1.i;
			}
			i__2 = jx;
			temp.r = X(jx).r, temp.i = X(jx).i;
			ix = jx;
			i__2 = *n;
			for (i = j + 1; i <= *n; ++i) {
			    ix += *incx;
			    i__3 = ix;
			    i__4 = ix;
			    i__5 = i + j * a_dim1;
			    q__2.r = temp.r * A(i,j).r - temp.i * A(i,j).i, 
				    q__2.i = temp.r * A(i,j).i + temp.i * A(i,j).r;
			    q__1.r = X(ix).r - q__2.r, q__1.i = X(ix).i - 
				    q__2.i;
			    X(ix).r = q__1.r, X(ix).i = q__1.i;
/* L70: */
			}
		    }
		    jx += *incx;
/* L80: */
		}
	    }
	}
    } else {

/*        Form  x := inv( A' )*x  or  x := inv( conjg( A' ) )*x. */

	if (strncmp(uplo, "U", 1)==0) {
	    if (*incx == 1) {
		i__1 = *n;
		for (j = 1; j <= *n; ++j) {
		    i__2 = j;
		    temp.r = X(j).r, temp.i = X(j).i;
		    if (noconj) {
			i__2 = j - 1;
			for (i = 1; i <= j-1; ++i) {
			    i__3 = i + j * a_dim1;
			    i__4 = i;
			    q__2.r = A(i,j).r * X(i).r - A(i,j).i * X(
				    i).i, q__2.i = A(i,j).r * X(i).i + 
				    A(i,j).i * X(i).r;
			    q__1.r = temp.r - q__2.r, q__1.i = temp.i - 
				    q__2.i;
			    temp.r = q__1.r, temp.i = q__1.i;
/* L90: */
			}
			if (nounit) {
			    c_div(&q__1, &temp, &A(j,j));
			    temp.r = q__1.r, temp.i = q__1.i;
			}
		    } else {
			i__2 = j - 1;
			for (i = 1; i <= j-1; ++i) {
			    r_cnjg(&q__3, &A(i,j));
			    i__3 = i;
			    q__2.r = q__3.r * X(i).r - q__3.i * X(i).i, 
				    q__2.i = q__3.r * X(i).i + q__3.i * X(
				    i).r;
			    q__1.r = temp.r - q__2.r, q__1.i = temp.i - 
				    q__2.i;
			    temp.r = q__1.r, temp.i = q__1.i;
/* L100: */
			}
			if (nounit) {
			    r_cnjg(&q__2, &A(j,j));
			    c_div(&q__1, &temp, &q__2);
			    temp.r = q__1.r, temp.i = q__1.i;
			}
		    }
		    i__2 = j;
		    X(j).r = temp.r, X(j).i = temp.i;
/* L110: */
		}
	    } else {
		jx = kx;
		i__1 = *n;
		for (j = 1; j <= *n; ++j) {
		    ix = kx;
		    i__2 = jx;
		    temp.r = X(jx).r, temp.i = X(jx).i;
		    if (noconj) {
			i__2 = j - 1;
			for (i = 1; i <= j-1; ++i) {
			    i__3 = i + j * a_dim1;
			    i__4 = ix;
			    q__2.r = A(i,j).r * X(ix).r - A(i,j).i * X(
				    ix).i, q__2.i = A(i,j).r * X(ix).i + 
				    A(i,j).i * X(ix).r;
			    q__1.r = temp.r - q__2.r, q__1.i = temp.i - 
				    q__2.i;
			    temp.r = q__1.r, temp.i = q__1.i;
			    ix += *incx;
/* L120: */
			}
			if (nounit) {
			    c_div(&q__1, &temp, &A(j,j));
			    temp.r = q__1.r, temp.i = q__1.i;
			}
		    } else {
			i__2 = j - 1;
			for (i = 1; i <= j-1; ++i) {
			    r_cnjg(&q__3, &A(i,j));
			    i__3 = ix;
			    q__2.r = q__3.r * X(ix).r - q__3.i * X(ix).i, 
				    q__2.i = q__3.r * X(ix).i + q__3.i * X(
				    ix).r;
			    q__1.r = temp.r - q__2.r, q__1.i = temp.i - 
				    q__2.i;
			    temp.r = q__1.r, temp.i = q__1.i;
			    ix += *incx;
/* L130: */
			}
			if (nounit) {
			    r_cnjg(&q__2, &A(j,j));
			    c_div(&q__1, &temp, &q__2);
			    temp.r = q__1.r, temp.i = q__1.i;
			}
		    }
		    i__2 = jx;
		    X(jx).r = temp.r, X(jx).i = temp.i;
		    jx += *incx;
/* L140: */
		}
	    }
	} else {
	    if (*incx == 1) {
		for (j = *n; j >= 1; --j) {
		    i__1 = j;
		    temp.r = X(j).r, temp.i = X(j).i;
		    if (noconj) {
			i__1 = j + 1;
			for (i = *n; i >= j+1; --i) {
			    i__2 = i + j * a_dim1;
			    i__3 = i;
			    q__2.r = A(i,j).r * X(i).r - A(i,j).i * X(
				    i).i, q__2.i = A(i,j).r * X(i).i + 
				    A(i,j).i * X(i).r;
			    q__1.r = temp.r - q__2.r, q__1.i = temp.i - 
				    q__2.i;
			    temp.r = q__1.r, temp.i = q__1.i;
/* L150: */
			}
			if (nounit) {
			    c_div(&q__1, &temp, &A(j,j));
			    temp.r = q__1.r, temp.i = q__1.i;
			}
		    } else {
			i__1 = j + 1;
			for (i = *n; i >= j+1; --i) {
			    r_cnjg(&q__3, &A(i,j));
			    i__2 = i;
			    q__2.r = q__3.r * X(i).r - q__3.i * X(i).i, 
				    q__2.i = q__3.r * X(i).i + q__3.i * X(
				    i).r;
			    q__1.r = temp.r - q__2.r, q__1.i = temp.i - 
				    q__2.i;
			    temp.r = q__1.r, temp.i = q__1.i;
/* L160: */
			}
			if (nounit) {
			    r_cnjg(&q__2, &A(j,j));
			    c_div(&q__1, &temp, &q__2);
			    temp.r = q__1.r, temp.i = q__1.i;
			}
		    }
		    i__1 = j;
		    X(j).r = temp.r, X(j).i = temp.i;
/* L170: */
		}
	    } else {
		kx += (*n - 1) * *incx;
		jx = kx;
		for (j = *n; j >= 1; --j) {
		    ix = kx;
		    i__1 = jx;
		    temp.r = X(jx).r, temp.i = X(jx).i;
		    if (noconj) {
			i__1 = j + 1;
			for (i = *n; i >= j+1; --i) {
			    i__2 = i + j * a_dim1;
			    i__3 = ix;
			    q__2.r = A(i,j).r * X(ix).r - A(i,j).i * X(
				    ix).i, q__2.i = A(i,j).r * X(ix).i + 
				    A(i,j).i * X(ix).r;
			    q__1.r = temp.r - q__2.r, q__1.i = temp.i - 
				    q__2.i;
			    temp.r = q__1.r, temp.i = q__1.i;
			    ix -= *incx;
/* L180: */
			}
			if (nounit) {
			    c_div(&q__1, &temp, &A(j,j));
			    temp.r = q__1.r, temp.i = q__1.i;
			}
		    } else {
			i__1 = j + 1;
			for (i = *n; i >= j+1; --i) {
			    r_cnjg(&q__3, &A(i,j));
			    i__2 = ix;
			    q__2.r = q__3.r * X(ix).r - q__3.i * X(ix).i, 
				    q__2.i = q__3.r * X(ix).i + q__3.i * X(
				    ix).r;
			    q__1.r = temp.r - q__2.r, q__1.i = temp.i - 
				    q__2.i;
			    temp.r = q__1.r, temp.i = q__1.i;
			    ix -= *incx;
/* L190: */
			}
			if (nounit) {
			    r_cnjg(&q__2, &A(j,j));
			    c_div(&q__1, &temp, &q__2);
			    temp.r = q__1.r, temp.i = q__1.i;
			}
		    }
		    i__1 = jx;
		    X(jx).r = temp.r, X(jx).i = temp.i;
		    jx -= *incx;
/* L200: */
		}
	    }
	}
    }

    return 0;

/*     End of CTRSV . */

} /* ctrsv_ */

