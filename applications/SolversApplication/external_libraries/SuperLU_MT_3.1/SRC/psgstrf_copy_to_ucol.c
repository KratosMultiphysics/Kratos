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
psgstrf_copy_to_ucol(
		     const int_t  pnum,    /* process number */
		     const int_t  jcol,	 /* current column */
		     const int_t  nseg,	 /* number of U-segments */
		     const int_t  *segrep, /* in */
		     const int_t  *repfnz, /* in */
		     const int_t  *perm_r, /* in */
		     float	 *dense,  /* modified - reset to zero on exit */
		     pxgstrf_shared_t *pxgstrf_shared /* modified */
		     )
{
/*
 * -- SuperLU MT routine (version 2.0) --
 * Lawrence Berkeley National Lab, Univ. of California Berkeley,
 * and Xerox Palo Alto Research Center.
 * September 10, 2007
 *
 * Gather the nonzeros from SPA dense[*,jcol] int_to global ucol[*].
 */
    register int_t ksub, krep, ksupno, i, k, kfnz, segsze;
    register int_t fsupc, isub, irow, jsupno, colsize;
    int_t      nextu, mem_error;
    int_t      *xsup, *supno, *lsub, *xlsub, *usub;
    float   *ucol;
    GlobalLU_t *Glu = pxgstrf_shared->Glu; /* modified */

    float zero = 0.0;

    xsup    = Glu->xsup;
    supno   = Glu->supno;
    lsub    = Glu->lsub;
    xlsub   = Glu->xlsub;
    jsupno  = supno[jcol];

    /* find the size of column jcol */
    colsize = 0;
    k = nseg - 1;
    for (ksub = 0; ksub < nseg; ++ksub) {
	krep = segrep[k--];
	ksupno = supno[krep];

	if ( ksupno != jsupno ) { /* should go into ucol[] */
	    kfnz = repfnz[krep];
	    if ( kfnz != EMPTY )  /* nonzero U-segment */
		colsize += krep - kfnz + 1;;
	}
    } /* for each segment... */

    if ( (mem_error = Glu_alloc(pnum, jcol, colsize, UCOL, &nextu, 
				pxgstrf_shared)) )
	return mem_error;
    Glu->xusub[jcol] = nextu;
    ucol = Glu->ucol;
    usub = Glu->usub;

    /* Now, it does not have to be in topological order! */
    k = nseg - 1;
    for (ksub = 0; ksub < nseg; ++ksub) {
	
	krep = segrep[k--];
	ksupno = supno[krep];

	if ( ksupno != jsupno ) { /* should go into ucol[] */
	    kfnz = repfnz[krep];
	    if ( kfnz != EMPTY ) { /* nonzero U-segment */
	    	fsupc = xsup[ksupno];
	        isub = xlsub[fsupc] + kfnz - fsupc;
	        segsze = krep - kfnz + 1;
#pragma ivdep
		for (i = 0; i < segsze; i++) {
		    irow = lsub[isub];
		    usub[nextu] = perm_r[irow];
		    ucol[nextu] = dense[irow];
		    dense[irow] = zero;
#ifdef DEBUG
if (jcol == EMPTY)
    printf("(%d) pcopy_to_ucol[]: jcol %d, krep %d, irow %d, ucol %.10e\n",
	   ME, jcol, krep, irow, ucol[nextu]);
#endif		    
		    nextu++;
		    isub++;
		} 
	    }
	}

    } /* for each segment... */

    Glu->xusub_end[jcol] = nextu; /* close U[*,jcol] */
    return 0;
}
