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
#include "slu_mt_ddefs.h"

void dlsolve(int_t, int_t, double *, double *);
void dmatvec(int_t, int_t, int_t, double *, double *, double *);


void
pdgstrf_bmod2D(
	       const int_t pnum,   /* process number */
	       const int_t m,      /* number of columns in the matrix */
	       const int_t w,      /* current panel width */
	       const int_t jcol,   /* leading column of the current panel */
	       const int_t fsupc,  /* leading column of the updating supernode */
	       const int_t krep,   /* last column of the updating supernode */
	       const int_t nsupc,  /* number of columns in the updating s-node */
	       int_t nsupr,        /* number of rows in the updating s-node */
	       int_t nrow,         /* number of rows below the diagonal block of
				    the updating supernode */
	       int_t *repfnz,      /* in */
	       int_t *panel_lsub,  /* modified */
	       int_t *w_lsub_end,  /* modified */
	       int_t *spa_marker,  /* modified; size n-by-w */
	       double *dense,    /* modified */
	       double *tempv,    /* working array - zeros on entry/exit */
	       GlobalLU_t *Glu,  /* modified */
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
 *
 *    Performs numeric 2-D block updates (sup-panel) in topological order.
 *    Results are returned in SPA dense[*,w].
 *
 */
#if ( MACH==CRAY_PVP )
    _fcd ftcs1 = _cptofcd("L", strlen("L")),
         ftcs2 = _cptofcd("N", strlen("N")),
         ftcs3 = _cptofcd("U", strlen("U"));
#endif
#ifdef USE_VENDOR_BLAS
    int          incx = 1, incy = 1;
    double      alpha, beta;
#endif
    double      zero = 0.0;
    double      one = 1.0;

    double       ukj, ukj1, ukj2;
    int_t          luptr, luptr1, luptr2;
    int          segsze, nsupr32 = nsupr;
    int          block_nrow;  /* no of rows in a block row */
    register int_t lptr;   /* point_ts to the row subscripts of a supernode */
    int_t          kfnz, irow, no_zeros; 
    register int_t isub, isub1, i;
    register int_t jj;	      /* index through each column in the panel */
    int_t          krep_ind;
    int_t          *repfnz_col; /* repfnz[] for a column in the panel */
    int_t          *col_marker; /* each column of the spa_marker[*,w] */
    int_t          *col_lsub;   /* each column of the panel_lsub[*,w] */
    double       *dense_col;  /* dense[] for a column in the panel */
    double       *TriTmp, *MatvecTmp;
    register int_t ldaTmp;
    register int_t r_ind, r_hi;
    static   int_t first = 1, maxsuper, rowblk;
    int_t          *lsub, *xlsub_end;
    double       *lusup;
    int_t          *xlusup;
    register float flopcnt;
    
#ifdef TIMING    
    double *utime = Gstat->utime;
    double f_time;
#endif    
    
    if ( first ) {
	maxsuper = sp_ienv(3);
	rowblk   = sp_ienv(4);
	first = 0;
    }
    ldaTmp = maxsuper + rowblk;

    lsub      = Glu->lsub;
    xlsub_end = Glu->xlsub_end;
    lusup     = Glu->lusup;
    xlusup    = Glu->xlusup;
    lptr      = Glu->xlsub[fsupc];
    krep_ind  = lptr + nsupc - 1;
    repfnz_col= repfnz;
    dense_col = dense;
    TriTmp    = tempv;
    col_marker= spa_marker;
    col_lsub  = panel_lsub;
	
	
    /* ---------------------------------------------------------------
     * Sequence through each column in the panel -- triangular solves.
     * The results of the triangular solves of all columns in the
     * panel are temporaroly stored in TriTemp array.
     * For the unrolled small supernodes of size <= 3, we also perform
     * matrix-vector updates from below the diagonal block.
     * ---------------------------------------------------------------
     */
    for (jj = jcol; jj < jcol + w; ++jj, col_marker += m, col_lsub += m,
	 repfnz_col += m, dense_col += m, TriTmp += ldaTmp ) {

	kfnz = repfnz_col[krep];
	if ( kfnz == EMPTY ) continue;	/* Skip any zero segment */
	    
	segsze = krep - kfnz + 1;
	luptr = xlusup[fsupc];

        flopcnt = segsze * (segsze - 1) + 2 * nrow * segsze;
	Gstat->procstat[pnum].fcops += flopcnt;

/*	ops[TRSV] += segsze * (segsze - 1);
	ops[GEMV] += 2 * nrow * segsze;        */
	
#ifdef TIMING	    
	f_time = SuperLU_timer_();
#endif
	
	/* Case 1: Update U-segment of size 1 -- col-col update */
	if ( segsze == 1 ) {
	    ukj = dense_col[lsub[krep_ind]];
	    luptr += nsupr*(nsupc-1) + nsupc;
	    for (i = lptr + nsupc; i < xlsub_end[fsupc]; i++) {
		irow = lsub[i];
                dense_col[irow] -= ukj * lusup[luptr];
		++luptr;
#ifdef SCATTER_FOUND		
		if ( col_marker[irow] != jj ) {
		    col_marker[irow] = jj;
		    col_lsub[w_lsub_end[jj-jcol]++] = irow;
		}
#endif		
	    }
#ifdef TIMING
	    utime[FLOAT] += SuperLU_timer_() - f_time;
#endif	    
	} else if ( segsze <= 3 ) {
	    ukj = dense_col[lsub[krep_ind]];
	    ukj1 = dense_col[lsub[krep_ind - 1]];
	    luptr += nsupr*(nsupc-1) + nsupc-1;
	    luptr1 = luptr - nsupr;
	    if ( segsze == 2 ) {
                ukj -= ukj1 * lusup[luptr1];
		dense_col[lsub[krep_ind]] = ukj;
		for (i = lptr + nsupc; i < xlsub_end[fsupc]; ++i) {
		    irow = lsub[i];
		    luptr++; luptr1++;
                    dense_col[irow] -= (ukj * lusup[luptr]
                                                + ukj1 * lusup[luptr1]);
#ifdef SCATTER_FOUND		
		    if ( col_marker[irow] != jj ) {
			col_marker[irow] = jj;
			col_lsub[w_lsub_end[jj-jcol]++] = irow;
		    }
#endif		
		}
#ifdef TIMING
		utime[FLOAT] += SuperLU_timer_() - f_time;
#endif	    
	    } else {
		ukj2 = dense_col[lsub[krep_ind - 2]];
		luptr2 = luptr1 - nsupr;
                ukj1 -= ukj2 * lusup[luptr2-1];
                ukj = ukj - ukj1*lusup[luptr1] - ukj2*lusup[luptr2];
		dense_col[lsub[krep_ind]] = ukj;
		dense_col[lsub[krep_ind-1]] = ukj1;
		for (i = lptr + nsupc; i < xlsub_end[fsupc]; ++i) {
		    irow = lsub[i];
		    luptr++; luptr1++; luptr2++;
                    dense_col[irow] -= (ukj * lusup[luptr]
                             + ukj1*lusup[luptr1] + ukj2*lusup[luptr2]);
#ifdef SCATTER_FOUND		
		    if ( col_marker[irow] != jj ) {
			col_marker[irow] = jj;
			col_lsub[w_lsub_end[jj-jcol]++] = irow;
		    }
#endif		
		}
	    }
#ifdef TIMING
	    utime[FLOAT] += SuperLU_timer_() - f_time;
#endif
	} else  { /* segsze >= 4 */
	    /* Copy A[*,j] segment from dense[*] to TriTmp[*], which
	       holds the result of triangular solve.    */
	    no_zeros = kfnz - fsupc;
	    isub = lptr + no_zeros;
	    for (i = 0; i < segsze; ++i) {
		irow = lsub[isub];
		TriTmp[i] = dense_col[irow]; /* Gather */
		++isub;
	    }

	    /* start effective triangle */
	    luptr += nsupr * no_zeros + no_zeros;
	    
#ifdef TIMING	    
	    f_time = SuperLU_timer_();
#endif
	    
#ifdef USE_VENDOR_BLAS
#if ( MACH==CRAY_PVP )
	    STRSV( ftcs1, ftcs2, ftcs3, &segsze, &lusup[luptr], 
		   &nsupr32, TriTmp, &incx );
#else
	    dtrsv_( "L", "N", "U", &segsze, &lusup[luptr], 
		   &nsupr32, TriTmp, &incx );
#endif
#else		
	    dlsolve ( nsupr, segsze, &lusup[luptr], TriTmp );
#endif
		
#ifdef TIMING	    
	    utime[FLOAT] += SuperLU_timer_() - f_time;
#endif	    
	} /* else ... */
	    
    }  /* for jj ... end tri-solves */

    /* --------------------------------------------------------
     * Perform block row updates from below the diagonal block.
     * Push each block all the way into SPA dense[*].
     * --------------------------------------------------------
     */
    for ( r_ind = 0; r_ind < nrow; r_ind += rowblk ) {
	    
	r_hi = SUPERLU_MIN(nrow, r_ind + rowblk);
	block_nrow = SUPERLU_MIN(rowblk, r_hi - r_ind);
	luptr = xlusup[fsupc] + nsupc + r_ind;
	isub1 = lptr + nsupc + r_ind;
	    
	repfnz_col = repfnz;
	TriTmp = tempv;
	dense_col = dense;
	col_marker= spa_marker;
	col_lsub  = panel_lsub;
	
	/* Sequence through each column in the panel -- matrix-vector */
	for (jj = jcol; jj < jcol + w; ++jj, col_marker += m, col_lsub += m,
	     repfnz_col += m, dense_col += m, TriTmp += ldaTmp) {

	    kfnz = repfnz_col[krep];
	    if ( kfnz == EMPTY ) continue; /* skip any zero segment */
	    
	    segsze = krep - kfnz + 1;
	    if ( segsze <= 3 ) continue;   /* skip unrolled cases */
		
	    /* Perform a block update, and scatter the result of
	       matrix-vector into SPA dense[*].		 */
	    no_zeros = kfnz - fsupc;
	    luptr1 = luptr + nsupr * no_zeros;
	    MatvecTmp = &TriTmp[maxsuper];
	    
#ifdef TIMING
	    f_time = SuperLU_timer_();
#endif	    
	    
#ifdef USE_VENDOR_BLAS
            alpha = one;
            beta = zero;
#if ( MACH==CRAY_PVP )
	    SGEMV( ftcs2, &block_nrow, &segsze, &alpha, &lusup[luptr], 
		  &nsupr32, TriTmp, &incx, &beta, MatvecTmp, &incy );
#else
	    dgemv_( "N", &block_nrow, &segsze, &alpha, &lusup[luptr1], 
		   &nsupr32, TriTmp, &incx, &beta, MatvecTmp, &incy );
#endif /* _CRAY_PVP */
#else
	    dmatvec( nsupr, (int_t) block_nrow, (int_t) segsze, &lusup[luptr1],
		    TriTmp, MatvecTmp);
#endif
		
#ifdef TIMING
	    utime[FLOAT] += SuperLU_timer_() - f_time;
#endif	    

	    /* Scatter MatvecTmp[*] into SPA dense[*] temporarily,
	     * such that MatvecTmp[*] can be re-used for the
	     * the next block row update. dense[] will be copied into 
	     * global store after the whole panel has been finished.
	     */
	    isub = isub1;
	    for (i = 0; i < block_nrow; i++) {
		irow = lsub[isub];
                dense_col[irow] -= MatvecTmp[i]; /* Scatter-add */
#ifdef SCATTER_FOUND		
		if ( col_marker[irow] != jj ) {
		    col_marker[irow] = jj;
		    col_lsub[w_lsub_end[jj-jcol]++] = irow;
		}
#endif		
		MatvecTmp[i] = zero;
		++isub;
	    }
	    
	} /* for jj ... */

    } /* for each block row ... */

    
    /* ------------------------------------------------
       Scatter the triangular solves into SPA dense[*].
       ------------------------------------------------ */
    repfnz_col = repfnz;
    TriTmp = tempv;
    dense_col = dense;
    
    for (jj = 0; jj < w; ++jj, repfnz_col += m, dense_col += m, 
	 TriTmp += ldaTmp) {
	kfnz = repfnz_col[krep];
	if ( kfnz == EMPTY ) continue; /* skip any zero segment */
	
	segsze = krep - kfnz + 1;
	if ( segsze <= 3 ) continue; /* skip unrolled cases */
	
	no_zeros = kfnz - fsupc;		
	isub = lptr + no_zeros;
	for (i = 0; i < segsze; i++) {
	    irow = lsub[isub];
	    dense_col[irow] = TriTmp[i]; /* Scatter */
	    TriTmp[i] = zero;
	    ++isub;
	}
    } /* for jj ... */
	
}
