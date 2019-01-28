/*! \file
Copyright (c) 2003, The Regents of the University of California, through
Lawrence Berkeley National Laboratory (subject to receipt of any required 
approvals from U.S. Dept. of Energy) 

All rights reserved. 

The source code is distributed under BSD license, see the file License.txt
at the top-level directory.
*/

#include <math.h>
#include <stdlib.h>
#include "slu_mt_sdefs.h"

int_t
psgstrf_pivotL(
	       const int_t  pnum,  /* process number */
	       const int_t  jcol,  /* current column */
	       const float u,   /* diagonal pivoting threshold */
	       yes_no_t *usepr,  /* re-use the pivot sequence given by
				    perm_r[]/inv_perm_r[] */
	       int_t   *perm_r,    /* modified - row pivotings */
	       int_t   *inv_perm_r,/* modified - inverse of perm_r */
	       int_t   *inv_perm_c,/* in - used to find diagonal of Pc*A*Pc' */
	       int_t   *pivrow,    /* the pivot row for this column */
	       GlobalLU_t *Glu,  /* modified - global LU data structures */
	       Gstat_t *Gstat    /* modified */
       )
{
/*
 * -- SuperLU MT routine (version 2.0) --
 * Lawrence Berkeley National Lab, Univ. of California Berkeley,
 * and Xerox Palo Alto Research Center.
 * September 10, 2007
 *
 * Purpose
 * =======
 *   Performs the numerical pivoting on the current column of L,
 *   and the CDIV operation.
 *
 * Pivot policy
 * ============
 *   (1) Compute thresh = u * max_(i>=j) abs(A_ij);
 *   (2) IF user specifies pivot row k and abs(A_kj) >= thresh THEN
 *           pivot row = k;
 *       ELSE IF abs(A_jj) >= thresh THEN
 *           pivot row = j;
 *       ELSE
 *           pivot row = m;
 * 
 *   Note: If you absolutely want to use a given pivot order, then set u=0.0.
 *
 * Return value
 * ============
 *   0      success;
 *   i > 0  U(i,i) is exactly zero.
 *
 */
    register int_t fsupc; /* first column in the supernode */
    register int_t nsupc; /* no of columns in the supernode */
    register int_t nsupr; /* no of rows in the supernode */
    register int_t lptr;  /* the starting subscript of the supernode */
    register int_t pivptr, old_pivptr, diag, diagind;
    register int_t isub, icol, k, itemp;
    register float pivmax, rtemp, thresh;
    float       *lu_sup_ptr; 
    float       *lu_col_ptr;
    float	temp;
    int_t          *lsub_ptr;
    int_t          *lsub;
    float       *lusup;
    int_t          *xlusup;

    /* Initialize pointers */
    lsub       = Glu->lsub;
    lusup      = Glu->lusup;
    xlusup     = Glu->xlusup;
    fsupc      = Glu->xsup[Glu->supno[jcol]];
    nsupc      = jcol - fsupc;	        /* excluding jcol; nsupc >= 0 */
    lptr       = Glu->xlsub[fsupc];
    nsupr      = Glu->xlsub_end[fsupc] - lptr;
    lu_sup_ptr = &lusup[xlusup[fsupc]];	/* start of the current supernode */
    lu_col_ptr = &lusup[xlusup[jcol]];	/* start of jcol in the supernode */
    lsub_ptr   = &lsub[lptr];	/* start of row indices of the supernode */

#ifdef CHK_PIVOT
    printf("Before cdiv: col %d\n", jcol);
    for (k = nsupc; k < nsupr; k++) 
	printf("  lu[%d] %f\n", lsub_ptr[k], lu_col_ptr[k]);
#endif
    
    /* Determine the largest abs numerical value for partial pivoting;
       Also search for user-specified pivot, and diagonal element. */
    if ( *usepr == YES ) *pivrow = inv_perm_r[jcol];
    diagind = inv_perm_c[jcol];
    pivmax = 0.0;
    pivptr = nsupc;
    diag = EMPTY;
    old_pivptr = nsupc;
    for (isub = nsupc; isub < nsupr; ++isub) {
        rtemp = fabs (lu_col_ptr[isub]);
	if ( rtemp > pivmax ) {
	    pivmax = rtemp;
	    pivptr = isub;
	}
	if ( *usepr == YES && lsub_ptr[isub] == *pivrow ) old_pivptr = isub;
	if ( lsub_ptr[isub] == diagind ) diag = isub;
    }

    /* Test for singularity */
    if ( pivmax == 0.0 ) {
	*pivrow = lsub_ptr[pivptr];
	perm_r[*pivrow] = jcol;
	inv_perm_r[jcol] = *pivrow;
	*usepr = NO;
	return (jcol+1);
    }

    thresh = u * pivmax;
    
    /* Choose appropriate pivotal element by our policy. */
    if ( *usepr == YES ) {
        rtemp = fabs (lu_col_ptr[old_pivptr]);
	if ( rtemp != 0.0 && rtemp >= thresh )
	    pivptr = old_pivptr;
	else
	    *usepr = NO;
    }
    if ( *usepr == NO ) {
	/* Can we use diagonal as pivot? */
	if ( diag >= 0 ) { /* diagonal exists */
            rtemp = fabs (lu_col_ptr[diag]);
            if ( rtemp != 0.0 && rtemp >= thresh ) pivptr = diag;
	}
	*pivrow = lsub_ptr[pivptr];
    }
    
    /* Record pivot row */
    perm_r[*pivrow] = jcol;
    inv_perm_r[jcol] = *pivrow;
    
    /* Interchange row subscripts */
    if ( pivptr != nsupc ) {
	itemp = lsub_ptr[pivptr];
	lsub_ptr[pivptr] = lsub_ptr[nsupc];
	lsub_ptr[nsupc] = itemp;

	/* Interchange numerical values as well, for the whole supernode,
	 * such that L is indexed the same way as A.
 	 */
	k = 0;
	for (icol = 0; icol <= nsupc; ++icol, k += nsupr) {
	    itemp = pivptr + k;
	    temp = lu_sup_ptr[itemp];
	    lu_sup_ptr[itemp] = lu_sup_ptr[nsupc + k];
	    lu_sup_ptr[nsupc + k] = temp;
	}
    } /* if */

    
    /* CDIV operation */
/*    ops[FACT] += nsupr - nsupc;*/
    Gstat->procstat[pnum].fcops += nsupr - nsupc;

    temp = 1.0 / lu_col_ptr[nsupc];
    for (k = nsupc+1; k < nsupr; k++)
        lu_col_ptr[k] *= temp;

#ifdef CHK_PIVOT
    printf("After cdiv: col %d\n", jcol);
    for (k = nsupc; k < nsupr; k++) 
	printf("  lu[%d] %f\n", lsub_ptr[k], lu_col_ptr[k]);
#endif

    return 0;
}

