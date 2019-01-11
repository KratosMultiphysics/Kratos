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
#ifndef __SUPERLU_DCOMPLEX /* allow multiple inclusions */
#define __SUPERLU_DCOMPLEX

/* 
 * This header file is to be included in source files z*.c
 */
#ifndef DCOMPLEX_INCLUDE
#define DCOMPLEX_INCLUDE

typedef struct { double r, i; } doublecomplex;


/* Macro definitions */

/* Complex Addition c = a + b */
#define z_add(c, a, b) { (c)->r = (a)->r + (b)->r; \
			 (c)->i = (a)->i + (b)->i; }

/* Complex Subtraction c = a - b */
#define z_sub(c, a, b) { (c)->r = (a)->r - (b)->r; \
			 (c)->i = (a)->i - (b)->i; }

/* Complex-Double Multiplication */
#define zd_mult(c, a, b) { (c)->r = (a)->r * (b); \
                           (c)->i = (a)->i * (b); }

/* Complex-Complex Multiplication */
#define zz_mult(c, a, b) { \
	double cr, ci; \
    	cr = (a)->r * (b)->r - (a)->i * (b)->i; \
    	ci = (a)->i * (b)->r + (a)->r * (b)->i; \
    	(c)->r = cr; \
    	(c)->i = ci; \
    }

#define zz_conj(a, b) { \
        (a)->r = (b)->r; \
        (a)->i = -((b)->i); \
    }

/* Complex equality testing */
#define z_eq(a, b)  ( (a)->r == (b)->r && (a)->i == (b)->i )


#ifdef __cplusplus
extern "C" {
#endif

/* Prototypes for functions in dcomplex.c */
void z_div(doublecomplex *, doublecomplex *, doublecomplex *);
double z_abs(doublecomplex *);     /* exact */
double z_abs1(doublecomplex *);    /* approximate */
void z_exp(doublecomplex *, doublecomplex *);
void d_cnjg(doublecomplex *r, doublecomplex *z);
double d_imag(doublecomplex *);


#ifdef __cplusplus
  }
#endif

#endif

#endif  /* __SUPERLU_DCOMPLEX */
