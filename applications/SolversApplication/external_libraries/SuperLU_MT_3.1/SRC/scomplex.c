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
 * This file defines common arithmetic operations for complex type.
 */
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "slu_scomplex.h"


/* Complex Division c = a/b */
void c_div(complex *c, complex *a, complex *b)
{
    float ratio, den;
    float abr, abi, cr, ci;
  
    if( (abr = b->r) < 0.)
	abr = - abr;
    if( (abi = b->i) < 0.)
	abi = - abi;
    if( abr <= abi ) {
	if (abi == 0) {
	    fprintf(stderr, "c_div: division by zero\n");
            exit(-1);
	}	  
	ratio = b->r / b->i ;
	den = b->i * (1 + ratio*ratio);
	cr = (a->r*ratio + a->i) / den;
	ci = (a->i*ratio - a->r) / den;
    } else {
	ratio = b->i / b->r ;
	den = b->r * (1 + ratio*ratio);
	cr = (a->r + a->i*ratio) / den;
	ci = (a->i - a->r*ratio) / den;
    }
    c->r = cr;
    c->i = ci;
}


/* Returns sqrt(z.r^2 + z.i^2) */
double c_abs(complex *z)
{
    float temp;
    float real = z->r;
    float imag = z->i;

    if (real < 0) real = -real;
    if (imag < 0) imag = -imag;
    if (imag > real) {
	temp = real;
	real = imag;
	imag = temp;
    }
    if ((real+imag) == real) return(real);
  
    temp = imag/real;
    temp = real*sqrt(1.0 + temp*temp);  /*overflow!!*/
    return (temp);
}


/* Approximates the abs */
/* Returns abs(z.r) + abs(z.i) */
double c_abs1(complex *z)
{
    float real = z->r;
    float imag = z->i;
  
    if (real < 0) real = -real;
    if (imag < 0) imag = -imag;

    return (real + imag);
}

/* Return the exponentiation */
void c_exp(complex *r, complex *z)
{
    float expx;

    expx = exp(z->r);
    r->r = expx * cos(z->i);
    r->i = expx * sin(z->i);
}

/* Return the complex conjugate */
void r_cnjg(complex *r, complex *z)
{
    r->r = z->r;
    r->i = -z->i;
}

/* Return the imaginary part */
double r_imag(complex *z)
{
    return (z->i);
}


