/*! \file
Copyright (c) 2003, The Regents of the University of California, through
Lawrence Berkeley National Laboratory (subject to receipt of any required 
approvals from U.S. Dept. of Energy) 

All rights reserved. 

The source code is distributed under BSD license, see the file License.txt
at the top-level directory.
*/
#include "slu_mt_ddefs.h"

void
pxgstrf_mark_busy_descends(int_t pnum, int_t jcol, int_t *etree, 
			   pxgstrf_shared_t *pxgstrf_shared,
			   int_t *bcol, int_t *lbusy)
{
/*
 * -- SuperLU MT routine (version 1.0) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * August 15, 1997
 *
 * Purpose
 * =======
 *
 *   Mark busy panels in local "lbusy" array, used for linear pipelining.
 *
 *   When jcol begins, its busy descendant panels (if any) are bcol and
 *   all the e-tree ancestors of bcol between bcol and jcol. This routine
 *   marks those columns in the array lbusy, which is local to this
 *   processor, to preserve a snapshot regardless of what the other
 *   processors are doing meanwhile.
 *
 * Arguments
 * =========
 *
 * jcol    (input) int_t
 *         Current panel, with leading column jcol.
 *
 * etree   (input) int_t*
 *         Elimination tree parent pointers.
 *
 * bcol    (input/output) int_t*
 *         Farthest busy descendant of jcol.
 *         On entry, it is the first column of the farthest busy panel.
 *         On exit, it may be adjusted to the first column of the
 *                  farthest busy supernode.
 *
 * lbusy   (input/output) int_t*
 *         Initially all -1, lbusy(r) = jcol means that r was busy
 *         at the beginning of computing jcol.
 *
 */
    GlobalLU_t *Glu = pxgstrf_shared->Glu;
    register int_t w,  kcol, fsupc, bcol_reg;
    int_t *xsup;

    bcol_reg = *bcol;
    if ( bcol_reg < jcol ) {
	
	/* -----------------------------------------------------------
	   Instead of waiting for the completion of "bcol", we can
	   pessimistically assume supno[bcol] == supno[bcol-1],
	   hence always mark as busy the supernode containing "bcol-1".
	   ----------------------------------------------------------- */
	if (pxgstrf_shared->pan_status[bcol_reg].type == RELAXED_SNODE) {
#if 0	    
	    if ( pxgstrf_shared->pan_status[bcol_reg].size < 0 )
	  	fsupc = bcol_reg + pxgstrf_shared->pan_status[bcol_reg].size;
	    else fsupc = bcol_reg;
#endif
	    fsupc = bcol_reg;
	    w = pxgstrf_shared->pan_status[fsupc].size;
	    bcol_reg += w;
	    for (kcol = fsupc; kcol < bcol_reg; ++kcol)
		lbusy[kcol] = jcol;
	} else {
	    /* Find leading column "fsupc" in the supernode that
	       contains column "bcol-1" */
#if 0
	    if ( pxgstrf_shared->spin_locks[bcol_reg] ) /* WORSE PERFORMANCE!! */
		await( &pxgstrf_shared->spin_locks[bcol_reg] );
#endif
	    xsup = Glu->xsup;
	    fsupc = SUPER_FSUPC ( Glu->supno[bcol_reg-1] );
	    for (kcol = fsupc; kcol < bcol_reg; ++kcol)	lbusy[kcol] = jcol;
	}
	
#if ( DEBUGlevel>=1 )
if (jcol >= LOCOL && jcol <= HICOL)
    printf("(%d) mark_busy_descends[1] jcol %d, bcol_reg %d, fsupc %d\n",
           pnum, jcol, bcol_reg, fsupc);
#endif
	
	/* Mark as busy all columns on the path between bcol_reg and jcol */
	for (kcol = bcol_reg; kcol < jcol; kcol = etree[kcol]) {
	    lbusy[kcol] = jcol;
	}

	/* INVARIANT: *bcol must be the first column of the farthest
	   busy supernode */
	*bcol = fsupc;
			 
    } /* if bcol_reg < jcol */
}
