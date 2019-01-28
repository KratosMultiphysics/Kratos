/*! \file
Copyright (c) 2003, The Regents of the University of California, through
Lawrence Berkeley National Laboratory (subject to receipt of any required 
approvals from U.S. Dept. of Energy) 

All rights reserved. 

The source code is distributed under BSD license, see the file License.txt
at the top-level directory.
*/

#include <stdio.h>
#include <stdlib.h>
#include "slu_mt_zdefs.h"

#define PRINT_SPIN_TIME(where)  { \
  if ( t2 > 0.001 ) { \
      printf("[" IFMT "] Panel " IFMT " on P" IFMT " waits s-node"IFMT " for %8.2f msec.\n",\
	     where, jcol, pnum, kcol, t2*1e3);\
      fflush(stdout); } }

static int_t
check_panel_dfs_list(int_t pnum, char *msg, int_t jcol, int_t nseg, int_t *segrep)
{
    register int_t k;

    printf("(" IFMT ") pzgstrf_panel_bmod(%s) jcol " IFMT ", nseg " IFMT " in top. order ->\n",
	   pnum, msg, jcol, nseg);
    for (k = nseg-1; k >= 0; --k)
	printf("(" IFMT ") segrep-" IFMT, pnum, segrep[k]);
    printf("\n");
    fflush(stdout);
    return 0;
}

#ifdef PREDICT_OPT

/* comparison function used by qsort() - in decreasing order */
int_t numcomp(desc_eft_t *a, desc_eft_t *b)
{
    if ( a->eft < b->eft )
	return -1;
    else if ( a->eft > b->eft )
	return 1;
    else
	return 0;
}
#endif

void
pzgstrf_panel_bmod(
		   const int_t  pnum, /* process number */
		   const int_t  m,    /* number of rows in the matrix */
		   const int_t  w,    /* current panel width */
		   const int_t  jcol, /* leading column of the current panel */
		   const int_t  bcol, /* first column of the farthest busy snode*/
		   int_t   *inv_perm_r,/* in; inverse of the row pivoting */
		   int_t   *etree,     /* in */
		   int_t   *nseg,      /* modified */
		   int_t   *segrep,    /* modified */
		   int_t   *repfnz,    /* modified, size n-by-w */
		   int_t   *panel_lsub,/* modified */
		   int_t   *w_lsub_end,/* modified */
		   int_t   *spa_marker,/* modified; size n-by-w */
		   doublecomplex *dense, /* modified, size n-by-w */
		   doublecomplex *tempv, /* working array - zeros on input/output */
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
 *
 *    Performs numeric block updates (sup-panel) in topological order.
 *    It features combined 1D and 2D blocking of the source updating s-node.
 *    It consists of two steps:
 *       (1) accumulates updates from "done" s-nodes.
 *       (2) accumulates updates from "busy" s-nodes.
 *
 *    Before entering this routine, the nonzeros of the original A in
 *    this panel were already copied into the SPA dense[n,w].
 *
 * Updated/Output arguments
 * ========================
 *    L[*,j:j+w-1] and U[*,j:j+w-1] are returned collectively in the
 *    m-by-w vector dense[*,w]. The locations of nonzeros in L[*,j:j+w-1]
 *    are given by lsub[*] and U[*,j:j+w-1] by (nseg,segrep,repfnz).
 *
 */
    GlobalLU_t *Glu = pxgstrf_shared->Glu;  /* modified */
    Gstat_t *Gstat = pxgstrf_shared->Gstat; /* modified */
    register int_t j, k, ksub;
    register int_t fsupc, nsupc, nsupr, nrow;
    register int_t kcol, krep, ksupno, dadsupno;
    register int_t jj;	      /* index through each column in the panel */
    int_t          *xsup, *xsup_end, *supno;
    int_t          *lsub, *xlsub, *xlsub_end;
    int_t          *repfnz_col; /* repfnz[] for a column in the panel */
    doublecomplex       *dense_col;  /* dense[] for a column in the panel */
    int_t          *col_marker; /* each column of the spa_marker[*,w] */
    int_t          *col_lsub;   /* each column of the panel_lsub[*,w] */
    /*    static   int first = 1, rowblk, colblk;*/
    register int_t rowblk = sp_ienv(4), colblk = sp_ienv(5);

#ifdef PROFILE
    double   t1, t2; /* temporary time */
#endif
    
#ifdef PREDICT_OPT    
    register float pmod, max_child_eft = 0, sum_pmod = 0, min_desc_eft = 0;
    register float pmod_eft;
    register int_t   kid, ndesc = 0;
#endif
    
#if ( DEBUGlevel>=2 )
    int_t dbg_addr = 0*m;
#endif

#if 0    
    if ( first ) {
	rowblk   = sp_ienv(4);
	colblk   = sp_ienv(5);
	first = 0;
    }
#endif
    
    xsup      = Glu->xsup;
    xsup_end  = Glu->xsup_end;
    supno     = Glu->supno;
    lsub      = Glu->lsub;
    xlsub     = Glu->xlsub;
    xlsub_end = Glu->xlsub_end;

#if ( DEBUGlevel>=2 )
    /*if (jcol >= LOCOL && jcol <= HICOL)
    check_panel_dfs_list(pnum, "begin", jcol, *nseg, segrep);*/
if (jcol == BADPAN)
    printf("(%d) Enter pzgstrf_panel_bmod() jcol %d,BADCOL %d,dense_col[%d] %.10f\n",
	   pnum, jcol, BADCOL, BADROW, dense[dbg_addr+BADROW]);
#endif    

    /* --------------------------------------------------------------------
       For each non-busy supernode segment of U[*,jcol] in topological order,
       perform sup-panel update.
       -------------------------------------------------------------------- */
    k = *nseg - 1;
    for (ksub = 0; ksub < *nseg; ++ksub) {
	/*
	 * krep = representative of current k-th supernode
	 * fsupc = first supernodal column
	 * nsupc = no of columns in a supernode
	 * nsupr = no of rows in a supernode
	 */
        krep = segrep[k--];
	fsupc = xsup[supno[krep]];
	nsupc = krep - fsupc + 1;
	nsupr = xlsub_end[fsupc] - xlsub[fsupc];
	nrow = nsupr - nsupc;

#ifdef PREDICT_OPT
	pmod = Gstat->procstat[pnum].fcops;
#endif
	    
	if ( nsupc >= colblk && nrow >= rowblk ) {
	    /* 2-D block update */
#ifdef GEMV2
	    pzgstrf_bmod2D_mv2(pnum, m, w, jcol, fsupc, krep, nsupc, nsupr, 
			       nrow, repfnz, panel_lsub, w_lsub_end, 
			       spa_marker, dense, tempv, Glu, Gstat);
#else
	    pzgstrf_bmod2D(pnum, m, w, jcol, fsupc, krep, nsupc, nsupr, nrow,
			   repfnz, panel_lsub, w_lsub_end, spa_marker,
			   dense, tempv, Glu, Gstat);
#endif
	} else {
	    /* 1-D block update */
#ifdef GEMV2
	    pzgstrf_bmod1D_mv2(pnum, m, w, jcol, fsupc, krep, nsupc, nsupr,
			       nrow, repfnz, panel_lsub, w_lsub_end, 
			       spa_marker, dense, tempv, Glu, Gstat);
#else
	    pzgstrf_bmod1D(pnum, m, w, jcol, fsupc, krep, nsupc, nsupr, nrow,
			   repfnz, panel_lsub, w_lsub_end, spa_marker,
			   dense, tempv, Glu, Gstat);
#endif
	}
	
#ifdef PREDICT_OPT
	pmod = Gstat->procstat[pnum].fcops - pmod;
	kid = (Glu->pan_status[krep].size > 0) ?
	    krep : (krep + Glu->pan_status[krep].size);
	desc_eft[ndesc].eft = cp_panel[kid].est + cp_panel[kid].pdiv;
	desc_eft[ndesc++].pmod = pmod;
#endif
	
#if ( DEBUGlevel>=2 )
if (jcol == BADPAN)
    printf("(%d) non-busy update: krep %d, repfnz %d, dense_col[%d] %.10e\n",
	   pnum, krep, repfnz[dbg_addr+krep], BADROW, dense[dbg_addr+BADROW]);
#endif

    } /* for each updating supernode ... */
    
#if ( DEBUGlevel>=2 )
if (jcol == BADPAN)
    printf("(%d) After non-busy update: dense_col[%d] %.10e\n",
	   pnum, BADROW, dense[dbg_addr+BADROW]);
#endif
    
    /* ---------------------------------------------------------------------
     * Now wait for the "busy" s-nodes to become "done" -- this amounts to
     * climbing up the e-tree along the path starting from "bcol".
     * Several points are worth noting:
     *
     *  (1) There are two possible relations between supernodes and panels
     *      along the path of the e-tree:
     *      o |s-node| < |panel|
     *        want to climb up the e-tree one column at a time in order
     *        to achieve more concurrency
     *      o |s-node| > |panel|
     *        want to climb up the e-tree one panel at a time; this
     *        processor is stalled anyway while waiting for the panel.
     *
     *  (2) Need to accommodate new fills, append them in panel_lsub[*,w].
     *      o use an n-by-w marker array, as part of the SPA (not scalable!)
     *
     *  (3) Symbolically, need to find out repfnz[S, w], for each (busy)
     *      supernode S.
     *      o use dense[inv_perm_r[kcol]], filter all zeros
     *      o detect the first nonzero in each segment
     *        (at this moment, the boundary of the busy supernode/segment
     *         S has already been identified)
     *
     * --------------------------------------------------------------------- */

    kcol = bcol;
    while ( kcol < jcol ) {
        /* Pointers to each column of the w-wide arrays. */
	repfnz_col = repfnz;
	dense_col = dense;
	col_marker = spa_marker;
	col_lsub = panel_lsub;

	/* Wait for the supernode, and collect wait-time statistics. */
	if ( pxgstrf_shared->spin_locks[kcol] ) {
#ifdef PROFILE
	    TIC(t1);
#endif
	    await( &pxgstrf_shared->spin_locks[kcol] );

#ifdef PROFILE
	    TOC(t2, t1);
	    Gstat->panstat[jcol].pipewaits++;
	    Gstat->panstat[jcol].spintime += t2;
	    Gstat->procstat[pnum].spintime += t2;
#ifdef DOPRINT
	    PRINT_SPIN_TIME(1);
#endif
#endif		
	}
	
        /* Find leading column "fsupc" in the supernode that
           contains column "kcol" */
	ksupno = supno[kcol];
	fsupc = kcol;

#if ( DEBUGlevel>=2 )
	/*if (jcol >= LOCOL && jcol <= HICOL)    */
  if ( jcol==BADCOL )
    printf("(%d) pzgstrf_panel_bmod[1] kcol %d, ksupno %d, fsupc %d\n",
	   pnum, kcol, ksupno, fsupc);
#endif
	
	/* Wait for the whole supernode to become "done" --
	   climb up e-tree one column at a time */
	do {
	    krep = SUPER_REP( ksupno );
	    kcol = etree[kcol];
	    if ( kcol >= jcol ) break;
	    if ( pxgstrf_shared->spin_locks[kcol] ) {
#ifdef PROFILE
		TIC(t1);
#endif
		await ( &pxgstrf_shared->spin_locks[kcol] );

#ifdef PROFILE
		TOC(t2, t1);
		Gstat->panstat[jcol].pipewaits++;
		Gstat->panstat[jcol].spintime += t2;
		Gstat->procstat[pnum].spintime += t2;
#ifdef DOPRINT
		PRINT_SPIN_TIME(2);
#endif
#endif		
	    }

	    dadsupno = supno[kcol];

#if ( DEBUGlevel>=2 )
	    /*if (jcol >= LOCOL && jcol <= HICOL)*/
if ( jcol==BADCOL )
    printf("(%d) pzgstrf_panel_bmod[2] krep %d, dad=kcol %d, dadsupno %d\n",
	   pnum, krep, kcol, dadsupno);
#endif	    

	} while ( dadsupno == ksupno );

	/* Append the new segment into segrep[*]. After column_bmod(),
	   copy_to_ucol() will use them. */
	segrep[*nseg] = krep;
        ++(*nseg);
        
	/* Determine repfnz[krep, w] for each column in the panel */
	for (jj = jcol; jj < jcol + w; ++jj, dense_col += m, 
	       repfnz_col += m, col_marker += m, col_lsub += m) {
	    /*
	     * Note: relaxed supernode may not form a path on the e-tree,
	     *       but its column numbers are contiguous.
	     */
#ifdef SCATTER_FOUND
 	    for (kcol = fsupc; kcol <= krep; ++kcol) {
		if ( col_marker[inv_perm_r[kcol]] == jj ) {
		    repfnz_col[krep] = kcol;

 		    /* Append new fills in panel_lsub[*,jj]. */
		    j = w_lsub_end[jj - jcol];
/*#pragma ivdep*/
		    for (k = xlsub[krep]; k < xlsub_end[krep]; ++k) {
			ksub = lsub[k];
			if ( col_marker[ksub] != jj ) {
			    col_marker[ksub] = jj;
			    col_lsub[j++] = ksub;
			}
		    }
		    w_lsub_end[jj - jcol] = j;

		    break; /* found the leading nonzero in the segment */
		}
	    }

#else
	    for (kcol = fsupc; kcol <= krep; ++kcol) {
		if ( (dense_col[inv_perm_r[kcol]].r != 0.0) || (dense_col[inv_perm_r[kcol]].i != 0.0) ) {
		    repfnz_col[krep] = kcol;
		    break; /* Found the leading nonzero in the U-segment */
		}
	    }

	    /* In this case, we always treat the L-subscripts of the 
	       busy s-node [kcol : krep] as the new fills, even if the
	       corresponding U-segment may be all zero. */

	    /* Append new fills in panel_lsub[*,jj]. */
	    j = w_lsub_end[jj - jcol];
/*#pragma ivdep*/
	    for (k = xlsub[krep]; k < xlsub_end[krep]; ++k) {
	        ksub = lsub[k];
		if ( col_marker[ksub] != jj ) {
		    col_marker[ksub] = jj;
		    col_lsub[j++] = ksub;
		}
	    }
	    w_lsub_end[jj - jcol] = j;
#endif

#if ( DEBUGlevel>=2 )
if (jj == BADCOL) {
printf("(%d) pzgstrf_panel_bmod[fills]: jj %d, repfnz_col[%d] %d, inv_pr[%d] %d\n",
	   pnum, jj, krep, repfnz_col[krep], fsupc, inv_perm_r[fsupc]);
printf("(%d) pzgstrf_panel_bmod[fills] xlsub %d, xlsub_end %d, #lsub[%d] %d\n",
       pnum,xlsub[krep],xlsub_end[krep],krep, xlsub_end[krep]-xlsub[krep]);
}
#endif	   
	} /* for jj ... */

#ifdef PREDICT_OPT
	pmod = Gstat->procstat[pnum].fcops;
#endif
	
	/* Perform sup-panel updates - use combined 1D + 2D updates. */
	nsupc = krep - fsupc + 1;
	nsupr = xlsub_end[fsupc] - xlsub[fsupc];
	nrow = nsupr - nsupc;
	if ( nsupc >= colblk && nrow >= rowblk ) {
	    /* 2-D block update */
#ifdef GEMV2
	    pzgstrf_bmod2D_mv2(pnum, m, w, jcol, fsupc, krep, nsupc, nsupr,
			       nrow, repfnz, panel_lsub, w_lsub_end, 
			       spa_marker, dense, tempv, Glu, Gstat);
#else
	    pzgstrf_bmod2D(pnum, m, w, jcol, fsupc, krep, nsupc, nsupr, nrow,
			   repfnz, panel_lsub, w_lsub_end, spa_marker,
			   dense, tempv, Glu, Gstat);
#endif
	} else {
	    /* 1-D block update */
#ifdef GEMV2
	    pzgstrf_bmod1D_mv2(pnum, m, w, jcol, fsupc, krep, nsupc, nsupr,
			       nrow, repfnz, panel_lsub, w_lsub_end, 
			       spa_marker, dense, tempv, Glu, Gstat);
#else
	    pzgstrf_bmod1D(pnum, m, w, jcol, fsupc, krep, nsupc, nsupr, nrow,
			   repfnz, panel_lsub, w_lsub_end, spa_marker,
			   dense, tempv, Glu, Gstat);
#endif
	}

#ifdef PREDICT_OPT
	pmod = Gstat->procstat[pnum].fcops - pmod;
	kid = (pxgstrf_shared->pan_status[krep].size > 0) ?
	       krep : (krep + pxgstrf_shared->pan_status[krep].size);
	desc_eft[ndesc].eft = cp_panel[kid].est + cp_panel[kid].pdiv;
	desc_eft[ndesc++].pmod = pmod;
#endif
	
#if ( DEBUGlevel>=2 )
if (jcol == BADPAN)
    printf("(%d) After busy update: dense_col[%d] %.10f\n",
	   pnum, BADROW, dense[dbg_addr+BADROW]);
#endif
	
	/* Go to the parent of "krep" */
	kcol = etree[krep];
	
    } /* while kcol < jcol ... */
    
#if ( DEBUGlevel>=2 )
    /*if (jcol >= LOCOL && jcol <= HICOL)*/
if ( jcol==BADCOL )
    check_panel_dfs_list(pnum, "after-busy", jcol, *nseg, segrep);
#endif

#ifdef PREDICT_OPT
    qsort(desc_eft, ndesc, sizeof(desc_eft_t), (int(*)())numcomp);
    pmod_eft = 0;
    for (j = 0; j < ndesc; ++j) {
	pmod_eft = SUPERLU_MAX( pmod_eft, desc_eft[j].eft ) + desc_eft[j].pmod;
    }

    if ( ndesc == 0 ) {
	/* No modifications from descendants */
	pmod_eft = 0;
	for (j = cp_firstkid[jcol]; j != EMPTY; j = cp_nextkid[j]) {
	    kid = (pxgstrf_shared->pan_status[j].size > 0) ? 
			j : (j + pxgstrf_shared->pan_status[j].size);
	    pmod_eft = SUPERLU_MAX( pmod_eft,
			   	cp_panel[kid].est + cp_panel[kid].pdiv );
	}
    }
    
    cp_panel[jcol].est = pmod_eft;
    
#endif

}

