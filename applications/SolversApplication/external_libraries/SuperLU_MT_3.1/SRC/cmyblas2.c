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
/*
 * File name:		cmyblas2.c
 * Purpose:
 *     Level 2 BLAS operations: solves and matvec, written in C.
 * Note:
 *     This is only used when the system lacks an efficient BLAS library.
 */
#include "slu_scomplex.h"
#include "slu_mt_cdefs.h"


/*
 * Solves a dense UNIT lower triangular system. The unit lower 
 * triangular matrix is stored in a 2D array M(1:nrow,1:ncol). 
 * The solution will be returned in the rhs vector.
 */
void clsolve ( int_t ldm, int_t ncol, complex *M, complex *rhs )
{
    int_t k;
    complex x0, x1, x2, x3, temp;
    complex *M0;
    complex *Mki0, *Mki1, *Mki2, *Mki3;
    register int_t firstcol = 0;

    M0 = &M[0];


    while ( firstcol < ncol - 3 ) { /* Do 4 columns */
      	Mki0 = M0 + 1;
      	Mki1 = Mki0 + ldm + 1;
      	Mki2 = Mki1 + ldm + 1;
      	Mki3 = Mki2 + ldm + 1;

      	x0 = rhs[firstcol];
      	cc_mult(&temp, &x0, Mki0); Mki0++;
      	c_sub(&x1, &rhs[firstcol+1], &temp);
      	cc_mult(&temp, &x0, Mki0); Mki0++;
	c_sub(&x2, &rhs[firstcol+2], &temp);
	cc_mult(&temp, &x1, Mki1); Mki1++;
	c_sub(&x2, &x2, &temp);
      	cc_mult(&temp, &x0, Mki0); Mki0++;
	c_sub(&x3, &rhs[firstcol+3], &temp);
	cc_mult(&temp, &x1, Mki1); Mki1++;
	c_sub(&x3, &x3, &temp);
	cc_mult(&temp, &x2, Mki2); Mki2++;
	c_sub(&x3, &x3, &temp);

 	rhs[++firstcol] = x1;
      	rhs[++firstcol] = x2;
      	rhs[++firstcol] = x3;
      	++firstcol;
    
      	for (k = firstcol; k < ncol; k++) {
	    cc_mult(&temp, &x0, Mki0); Mki0++;
	    c_sub(&rhs[k], &rhs[k], &temp);
	    cc_mult(&temp, &x1, Mki1); Mki1++;
	    c_sub(&rhs[k], &rhs[k], &temp);
	    cc_mult(&temp, &x2, Mki2); Mki2++;
	    c_sub(&rhs[k], &rhs[k], &temp);
	    cc_mult(&temp, &x3, Mki3); Mki3++;
	    c_sub(&rhs[k], &rhs[k], &temp);
	}

        M0 += 4 * ldm + 4;
    }

    if ( firstcol < ncol - 1 ) { /* Do 2 columns */
        Mki0 = M0 + 1;
        Mki1 = Mki0 + ldm + 1;

        x0 = rhs[firstcol];
	cc_mult(&temp, &x0, Mki0); Mki0++;
	c_sub(&x1, &rhs[firstcol+1], &temp);

      	rhs[++firstcol] = x1;
      	++firstcol;
    
      	for (k = firstcol; k < ncol; k++) {
	    cc_mult(&temp, &x0, Mki0); Mki0++;
	    c_sub(&rhs[k], &rhs[k], &temp);
	    cc_mult(&temp, &x1, Mki1); Mki1++;
	    c_sub(&rhs[k], &rhs[k], &temp);
	} 
    }
    
}

/*
 * Solves a dense upper triangular system. The upper triangular matrix is
 * stored in a 2-dim array M(1:ldm,1:ncol). The solution will be returned
 * in the rhs vector.
 */
void
cusolve (
int_t ldm,	/* in */
int_t ncol,	/* in */
complex *M,	/* in */
complex *rhs	/* modified */
)
{
    complex xj, temp;
    int_t jcol, j, irow;

    jcol = ncol - 1;

    for (j = 0; j < ncol; j++) {

	c_div(&xj, &rhs[jcol], &M[jcol + jcol*ldm]); /* M(jcol, jcol) */
	rhs[jcol] = xj;
	
	for (irow = 0; irow < jcol; irow++) {
	    cc_mult(&temp, &xj, &M[irow+jcol*ldm]); /* M(irow, jcol) */
	    c_sub(&rhs[irow], &rhs[irow], &temp);
	}

	jcol--;

    }
}


/*
 * Performs a dense matrix-vector multiply: Mxvec = Mxvec + M * vec.
 * The input matrix is M(1:nrow,1:ncol); The product is returned in Mxvec[].
 */
void cmatvec (
int_t ldm,	/* in -- leading dimension of M */
int_t nrow,	/* in */ 
int_t ncol,	/* in */
complex *M,	/* in */
complex *vec,	/* in */
complex *Mxvec	/* in/out */
)
{
    complex vi0, vi1, vi2, vi3;
    complex *M0, temp;
    complex *Mki0, *Mki1, *Mki2, *Mki3;
    register int_t firstcol = 0;
    int_t k;

    M0 = &M[0];

    while ( firstcol < ncol - 3 ) {	/* Do 4 columns */
	Mki0 = M0;
	Mki1 = Mki0 + ldm;
	Mki2 = Mki1 + ldm;
	Mki3 = Mki2 + ldm;

	vi0 = vec[firstcol++];
	vi1 = vec[firstcol++];
	vi2 = vec[firstcol++];
	vi3 = vec[firstcol++];	
	for (k = 0; k < nrow; k++) {
	    cc_mult(&temp, &vi0, Mki0); Mki0++;
	    c_add(&Mxvec[k], &Mxvec[k], &temp);
	    cc_mult(&temp, &vi1, Mki1); Mki1++;
	    c_add(&Mxvec[k], &Mxvec[k], &temp);
	    cc_mult(&temp, &vi2, Mki2); Mki2++;
	    c_add(&Mxvec[k], &Mxvec[k], &temp);
	    cc_mult(&temp, &vi3, Mki3); Mki3++;
	    c_add(&Mxvec[k], &Mxvec[k], &temp);
	}

	M0 += 4 * ldm;
    }

    while ( firstcol < ncol ) {		/* Do 1 column */
 	Mki0 = M0;
	vi0 = vec[firstcol++];
	for (k = 0; k < nrow; k++) {
	    cc_mult(&temp, &vi0, Mki0); Mki0++;
	    c_add(&Mxvec[k], &Mxvec[k], &temp);
	}
	M0 += ldm;
    }
	
}

/*
 * Performs dense matrix-vector multiply with 2 vectors:
 *        y0 = y0 + A * x0
 *        y1 = y1 + A * x1
 */
void cmatvec2 (
               int_t lda,     /* leading dimension of A */
               int_t m,
               int_t n,
               complex *A,   /* in - size m-by-n */
               complex *x0,  /* in - size n-by-1 */
               complex *x1,  /* in - size n-by-1 */
               complex *y0,  /* out - size n-by-1 */
               complex *y1   /* out - size n-by-1 */
               )
{
    complex v00, v10, v20, v30, v01, v11, v21, v31;
    complex *M0, temp, t0, t1, t2, t3;
    complex *Mki0, *Mki1, *Mki2, *Mki3;
    register int_t firstcol = 0;
    int_t k;

    M0 = &A[0];

    while ( firstcol < n - 3 ) {	/* Do 4 columns */
	Mki0 = M0;
	Mki1 = Mki0 + lda;
	Mki2 = Mki1 + lda;
	Mki3 = Mki2 + lda;

        v00 = x0[firstcol];   v01 = x1[firstcol++];
        v10 = x0[firstcol];   v11 = x1[firstcol++];
        v20 = x0[firstcol];   v21 = x1[firstcol++];
        v30 = x0[firstcol];   v31 = x1[firstcol++];

	for (k = 0; k < m; k++) {
	    t0 = Mki0[k];
	    t1 = Mki1[k];
	    t2 = Mki2[k];
	    t3 = Mki3[k];
	    
	    cc_mult(&temp, &v00, &t0);
	    c_add(&y0[k], &y0[k], &temp);
	    cc_mult(&temp, &v01, &t0);
	    c_add(&y1[k], &y1[k], &temp);

	    cc_mult(&temp, &v10, &t1);
	    c_add(&y0[k], &y0[k], &temp);
	    cc_mult(&temp, &v11, &t1);
	    c_add(&y1[k], &y1[k], &temp);

	    cc_mult(&temp, &v20, &t2);
	    c_add(&y0[k], &y0[k], &temp);
	    cc_mult(&temp, &v21, &t2);
	    c_add(&y1[k], &y1[k], &temp);

	    cc_mult(&temp, &v30, &t3);
	    c_add(&y0[k], &y0[k], &temp);
	    cc_mult(&temp, &v31, &t3);
	    c_add(&y1[k], &y1[k], &temp);
	}

	M0 += 4 * lda;
    }

    while ( firstcol < n ) {		/* Do 1 column */
 	Mki0 = M0;
        v00 = x0[firstcol];   v01 = x1[firstcol++];

	for (k = 0; k < m; k++) {
	    t0 = Mki0[k];
	    
	    cc_mult(&temp, &v00, &t0);
	    c_add(&y0[k], &y0[k], &temp);
	    cc_mult(&temp, &v01, &t0);
	    c_add(&y1[k], &y1[k], &temp);
	}
	M0 += lda;
    }
	
}

