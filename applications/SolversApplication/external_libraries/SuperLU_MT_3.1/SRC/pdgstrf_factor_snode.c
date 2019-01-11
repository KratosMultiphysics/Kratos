/*! \file
Copyright (c) 2003, The Regents of the University of California, through
Lawrence Berkeley National Laboratory (subject to receipt of any required 
approvals from U.S. Dept. of Energy) 

All rights reserved. 

The source code is distributed under BSD license, see the file License.txt
at the top-level directory.
*/

#include "slu_mt_ddefs.h"

int_t
pdgstrf_factor_snode(
		     const int_t pnum,  /* process number */
		     const int_t jcol,
		     SuperMatrix *A,
		     const double diag_pivot_thresh,
		     yes_no_t *usepr,
		     int_t    *perm_r,
		     int_t    *inv_perm_r, /* modified */
		     int_t    *inv_perm_c, /* in - used to find diagonal of Pc*A*Pc' */
		     int_t    *xprune,
		     int_t    *marker,
		     int_t    *col_lsub, /* values are irrevelant on entry 
					  and on return */
		     double *dense,
		     double *tempv,
		     pxgstrf_shared_t *pxgstrf_shared,
		     int_t    *info
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
 *
 *   Factorize the artificial supernodes grouped at the bottom
 *   of the etree.
 *
 */
    GlobalLU_t   *Glu = pxgstrf_shared->Glu;
    int_t          singular;
    NCPformat    *Astore;
    register int_t kcol, icol, k, jsupno, fsupc, nsupr;
    register int_t ifrom, ito;
    int_t          nextu, nextlu;
    int_t          pivrow;
    double       *a;
    int_t          *asub, *xa_begin, *xa_end, *xusub, *xusub_end,
                 *xsup, *supno, *xlusup, *lsub, *xlsub, *xlsub_end;

    lsub      = Glu->lsub;
    xlsub     = Glu->xlsub;
    xlsub_end = Glu->xlsub_end;
    xusub     = Glu->xusub;
    xusub_end = Glu->xusub_end;
    xsup      = Glu->xsup;
    supno     = Glu->supno;
    xlusup    = Glu->xlusup;
    
    singular = 0;
    Astore   = A->Store;
    a        = Astore->nzval;
    asub     = Astore->rowind;
    xa_begin = Astore->colbeg;
    xa_end   = Astore->colend;
    
    kcol = jcol + pxgstrf_shared->pan_status[jcol].size;
	
    /* Determine the union of the row structure of the supernode */
    if ( (*info = pdgstrf_snode_dfs(pnum, jcol, kcol-1, asub, xa_begin, xa_end,
				   xprune, marker, col_lsub, pxgstrf_shared)) )
	return 0;
    
    /*
     * Factorize the relaxed supernode (jcol:kcol-1)
     */
    nextu        = Glu->nextu; /* xiaoye - race condition (no problem!) */
    jsupno       = supno[jcol];
    fsupc        = xsup[jsupno];
    nsupr        = xlsub_end[fsupc] - xlsub[fsupc];
    if ( (*info = Glu_alloc(pnum, jcol, nsupr*(kcol-jcol), LUSUP, &nextlu,
			  pxgstrf_shared)) )
	return 0;
    
    for (icol = jcol; icol < kcol; icol++) {
	xusub[icol] = xusub_end[icol] = nextu;
	xlusup[icol] = nextlu;
	
	/* Scatter into SPA dense[*] */
	for (k = xa_begin[icol]; k < xa_end[icol]; k++)
	    dense[asub[k]] = a[k];
	
	/* Numeric update within the supernode */
	pdgstrf_snode_bmod(pnum, icol, jsupno, fsupc, dense, tempv, 
			   Glu, pxgstrf_shared->Gstat);
	
	if ( (*info = pdgstrf_pivotL
	                 (pnum, icol, diag_pivot_thresh, usepr, perm_r,
			  inv_perm_r, inv_perm_c, &pivrow, 
			  Glu, pxgstrf_shared->Gstat)) )
	    if ( singular == 0 ) singular = *info;
	
	nextlu += nsupr;

#if ( DEBUGlevel>= 2 )
  if ( icol>=LOCOL && icol<=HICOL )
    dprint_lu_col(pnum,"relax:",jcol,icol,kcol-jcol,pivrow,xprune,Glu);
#endif
	
    }

    /* Store the row subscripts of kcol-1 for pruned graph */
    k = ito = xlsub_end[jcol];
    for (ifrom = xlsub[jcol]+kcol-jcol-1; ifrom < k; ++ifrom)
	lsub[ito++] = lsub[ifrom];
    k = ito;
    xprune[kcol-1] = k;
    if (jcol < kcol-1) {    /* not a singleton */
	for (icol = jcol+1; icol < kcol; ++icol) xlsub_end[icol] = k;
	k = xlsub_end[jcol];
	xprune[jcol] = k;
	for (icol = jcol+1; icol < kcol; ++icol) xlsub[icol] = k;
    }
    
    *info = singular;
    return 0;
}
