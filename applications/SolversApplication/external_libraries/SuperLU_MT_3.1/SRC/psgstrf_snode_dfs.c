/*! \file
Copyright (c) 2003, The Regents of the University of California, through
Lawrence Berkeley National Laboratory (subject to receipt of any required 
approvals from U.S. Dept. of Energy) 

All rights reserved. 

The source code is distributed under BSD license, see the file License.txt
at the top-level directory.
*/

#include "slu_mt_sdefs.h"

int_t
psgstrf_snode_dfs(
		  const int_t  pnum,      /* process number */
		  const int_t  jcol,	  /* in - start of the supernode */
		  const int_t  kcol, 	  /* in - end of the supernode */
		  const int_t  *asub,     /* in */
		  const int_t  *xa_begin, /* in */
		  const int_t  *xa_end,   /* in */
		  int_t        *xprune,   /* out */
		  int_t        *marker,   /* modified */
		  int_t        *col_lsub, /* values are irrelevant on entry 
					   and on return */
		  pxgstrf_shared_t *pxgstrf_shared /* modified */
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
 *    psgstrf_snode_dfs() determines the union of the row structures of 
 *    those columns within the relaxed snode.
 *    Note: The relaxed snodes are leaves of the supernodal etree, 
 *    therefore, the portion outside the rectangular supernode must be zero.
 *
 * Return value
 * ============
 *     0   success;
 *    >0   number of bytes allocated when run out of memory.
 *
 */
    GlobalLU_t *Glu = pxgstrf_shared->Glu;
    register int_t i, k, ifrom, nextl, nsuper;
    int_t          ito;
    int_t          krow, kmark, mem_error;
    int_t          *supno, *lsub, *xlsub, *xlsub_end;
    
    supno                 = Glu->supno;
    xlsub                 = Glu->xlsub;
    xlsub_end             = Glu->xlsub_end;
    nsuper = NewNsuper(pnum, pxgstrf_shared, &Glu->nsuper);
    Glu->xsup[nsuper]     = jcol;
    Glu->xsup_end[nsuper] = kcol + 1;
    
    nextl = 0;
    for (i = jcol; i <= kcol; i++) {
	/* for each nonzero in A[*,i] */
	for (k = xa_begin[i]; k < xa_end[i]; k++) {	
	    krow = asub[k];
	    kmark = marker[krow];
	    if ( kmark != kcol ) { /* First time visit krow */
		marker[krow] = kcol;
		col_lsub[nextl++] = krow;
	    }
    	}
	supno[i] = nsuper;
    }

    if ( (mem_error = Glu_alloc(pnum, jcol, 2*nextl, LSUB, &ito, 
				pxgstrf_shared)) )
	return mem_error;
    
    xlsub[jcol] = ito;
    lsub        = Glu->lsub;
    for (ifrom = 0; ifrom < nextl; ++ifrom)
	lsub[ito++] = col_lsub[ifrom];
    xlsub_end[jcol] = ito;

    return 0;
}
