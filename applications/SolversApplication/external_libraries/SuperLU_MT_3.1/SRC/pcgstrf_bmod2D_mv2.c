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
#include "slu_mt_cdefs.h"

void clsolve(int_t, int_t, complex *, complex *);
void cmatvec(int_t, int_t, int_t, complex *, complex *, complex *);
void cmatvec2(int_t, int_t, int_t, complex*, complex*, complex*, complex*, complex*);


void
pcgstrf_bmod2D_mv2(
		   const int_t pnum, /* process number */
		   const int_t n,    /* number of rows in the matrix */
		   const int_t w,    /* current panel width */
		   const int_t jcol, /* leading column of the current panel */
		   const int_t fsupc,/* leading column of the updating supernode */
		   const int_t krep, /* last column of the updating s-node */
		   const int_t nsupc,/* number of columns in the updating s-node */
		   int_t nsupr, /* number of rows in the updating s-node */
		   int_t nrow,  /* number of rows below the diagonal block of
				 the updating supernode */
		   int_t *repfnz,     /* in */
		   int_t *panel_lsub, /* modified */
		   int_t *w_lsub_end, /* modified */
		   int_t *spa_marker, /* modified; size n-by-w */
		   complex *dense,   /* modified */
		   complex *tempv,   /* working array - zeros on entry/exit */
		   GlobalLU_t *Glu, /* modified */
		   Gstat_t *Gstat   /* modified */
		   )
{
/*
 * -- SuperLU MT routine (version 3.0) --
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

    complex      zero = {0.0, 0.0};
    complex      one = {1.0, 0.0};
    complex      comp_temp, comp_temp1;

#if ( MACH==CRAY_PVP )
    _fcd ftcs1 = _cptofcd("L", strlen("L")),
         ftcs2 = _cptofcd("N", strlen("N")),
         ftcs3 = _cptofcd("U", strlen("U"));
#endif
#ifdef USE_VENDOR_BLAS
    int          incx = 1, incy = 1;
    complex       alpha = one, beta = zero;
#endif

    complex       ukj, ukj1, ukj2;
    int_t          luptr, luptr1, luptr2;
    int          segsze, nsupr32 = nsupr;
    int          block_nrow;  /* no of rows in a block row */
    register int_t lptr; /* points to the row subscripts of a supernode */
    int_t          kfnz, irow, no_zeros; 
    register int_t isub, isub1, i, j;
    register int_t jj;	      /* index through each column in the panel */
    int_t          krep_ind;
    int_t          *repfnz_col; /* repfnz[] for a column in the panel */
    int_t          *col_marker; /* each column of the spa_marker[*,w] */
    int_t          *col_lsub;   /* each column of the panel_lsub[*,w] */
    complex       *dense_col;  /* dense[] for a column in the panel */
    complex       *TriTmp;
    register int_t ldaTmp;
    register int_t r_ind, r_hi;
    static   int_t first = 1, maxsuper, rowblk;
    register int_t twocols;
    int_t          kfnz2[2], jj2[2]; /* detect two identical columns */
    complex       *tri[2], *matvec[2];
    int_t          *lsub, *xlsub_end;
    complex       *lusup;
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

    
    /* ---------------------------------------------------------------
     * Sequence through each column in the panel -- triangular solves.
     * The results of the triangular solves of all columns in the
     * panel are temporaroly stored in TriTemp[*] array.
     * For the unrolled small supernodes of size <= 3, we also perform
     * matrix-vector updates from below the diagonal block.
     * ---------------------------------------------------------------
     */
    repfnz_col= repfnz;
    dense_col = dense;
    TriTmp    = tempv;
    col_marker= spa_marker;
    col_lsub  = panel_lsub;
    
    for (jj = jcol; jj < jcol + w; ++jj, col_marker += n, col_lsub += n,
	 repfnz_col += n, dense_col += n, TriTmp += ldaTmp ) {

	kfnz = repfnz_col[krep];
	if ( kfnz == EMPTY ) continue;	/* Skip any zero segment */
	    
	segsze = krep - kfnz + 1;
	luptr = xlusup[fsupc];

	Gstat->procstat[pnum].fcops += flopcnt;

#ifdef TIMING	    
	f_time = SuperLU_timer_();
#endif
	
	/* Case 1: Update U-segment of size 1 -- col-col update */
	if ( segsze == 1 ) {
	    ukj = dense_col[lsub[krep_ind]];
	    luptr += nsupr*(nsupc-1) + nsupc;
	    for (i = lptr + nsupc; i < xlsub_end[fsupc]; i++) {
		irow = lsub[i];
                cc_mult(&comp_temp, &ukj, &lusup[luptr]);
                c_sub(&dense_col[irow], &dense_col[irow], &comp_temp);
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
                cc_mult(&comp_temp, &ukj1, &lusup[luptr1]);
                c_sub(&ukj, &ukj, &comp_temp);
		dense_col[lsub[krep_ind]] = ukj;
		for (i = lptr + nsupc; i < xlsub_end[fsupc]; ++i) {
		    irow = lsub[i];
		    luptr++; luptr1++;
                    cc_mult(&comp_temp, &ukj, &lusup[luptr]);
                    cc_mult(&comp_temp1, &ukj1, &lusup[luptr1]);
                    c_add(&comp_temp, &comp_temp, &comp_temp1);
                    c_sub(&dense_col[irow], &dense_col[irow], &comp_temp);
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
                cc_mult(&comp_temp, &ukj2, &lusup[luptr2-1]);
                c_sub(&ukj1, &ukj1, &comp_temp);

                cc_mult(&comp_temp, &ukj1, &lusup[luptr1]);
                cc_mult(&comp_temp1, &ukj2, &lusup[luptr2]);
                c_add(&comp_temp, &comp_temp, &comp_temp1);
                c_sub(&ukj, &ukj, &comp_temp);
		dense_col[lsub[krep_ind]] = ukj;
		dense_col[lsub[krep_ind-1]] = ukj1;
		for (i = lptr + nsupc; i < xlsub_end[fsupc]; ++i) {
		    irow = lsub[i];
		    luptr++; luptr1++; luptr2++;
                    cc_mult(&comp_temp, &ukj, &lusup[luptr]);
                    cc_mult(&comp_temp1, &ukj1, &lusup[luptr1]);
                    c_add(&comp_temp, &comp_temp, &comp_temp1);
                    cc_mult(&comp_temp1, &ukj2, &lusup[luptr2]);
                    c_add(&comp_temp, &comp_temp, &comp_temp1);
                    c_sub(&dense_col[irow], &dense_col[irow], &comp_temp);
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
	    CTRSV( ftcs1, ftcs2, ftcs3, &segsze, &lusup[luptr], 
		   &nsupr32, TriTmp, &incx );
#else
	    ctrsv_( "L", "N", "U", &segsze, &lusup[luptr], &nsupr32, TriTmp, &incx );
#endif
#else		
	    clsolve ( (int_t) nsupr, (int_t) segsze, &lusup[luptr], TriTmp );
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
	luptr1 = xlusup[fsupc] + nsupc + r_ind;
	isub1 = lptr + nsupc + r_ind;
	repfnz_col = repfnz;
	twocols = 0;
	
	/* Sequence through each column in the panel -- matrix-vector */
	for (jj = jcol; jj < jcol + w; ++jj, repfnz_col += n) {

	    kfnz = repfnz_col[krep];
	    if ( kfnz == EMPTY ) continue; /* skip zero segment */
	    segsze = krep - kfnz + 1;
	    if ( segsze <= 3 ) continue;   /* skip unrolled cases */

	    /* Now segsze >= 4 ... */
	    
	    if ( twocols == 1 ) { /* got two columns */
		jj2[1] = jj;
		twocols = 0;
		for (j = 0; j < 2; ++j) {
		    i = n * (jj2[j] - jcol);
		    kfnz2[j] = repfnz[i + krep];
		    tri[j] = tempv + ldaTmp * (jj2[j] - jcol);
		    matvec[j] = tri[j] + maxsuper;
		}
		
		if ( kfnz2[0] < kfnz2[1] ) { /* First column is bigger */
		    no_zeros = kfnz2[0] - fsupc;
		    segsze = kfnz2[1] - kfnz2[0];
		    luptr = luptr1 + nsupr * no_zeros;
#ifdef USE_VENDOR_BLAS
#if ( MACH==CRAY_PVP )
		    CGEMV( ftcs2, &block_nrow, &segsze, &alpha, &lusup[luptr], 
			   &nsupr32, tri[0], &incx, &beta, matvec[0], &incy );
#else
		    cgemv_( "N", &block_nrow, &segsze, &alpha, &lusup[luptr], 
			   &nsupr32, tri[0], &incx, &beta, matvec[0], &incy );
#endif
#else
		    cmatvec ((int_t) nsupr, (int_t) block_nrow, (int_t) segsze, &lusup[luptr],
			     tri[0], matvec[0]);
#endif
		} else if ( kfnz2[0] > kfnz2[1] ) {
		    no_zeros = kfnz2[1] - fsupc;
		    segsze = kfnz2[0] - kfnz2[1];
		    luptr = luptr1 + nsupr * no_zeros;
#ifdef USE_VENDOR_BLAS
#if ( MACH==CRAY_PVP )
		    CGEMV( ftcs2, &block_nrow, &segsze, &alpha, &lusup[luptr], 
			   &nsupr, tri[1], &incx, &beta, matvec[1], &incy );
#else
		    cgemv_( "N", &block_nrow, &segsze, &alpha, &lusup[luptr], 
			   &nsupr32, tri[1], &incx, &beta, matvec[1], &incy );
#endif
#else
		    cmatvec (nsupr, block_nrow, segsze, &lusup[luptr],
			     tri[1], matvec[1]);
#endif
		}
		
		/* Do matrix-vector multiply with two destinations */
		kfnz = SUPERLU_MAX( kfnz2[0], kfnz2[1] );
		no_zeros = kfnz - fsupc;
		segsze = krep - kfnz + 1;
		luptr = luptr1 + nsupr * no_zeros;
#if ( MACH==DEC )
		cgemv2_ (&nsupr, &block_nrow, &segsze, &lusup[luptr],
			  &tri[0][kfnz-kfnz2[0]], &tri[1][kfnz-kfnz2[1]],
			  matvec[0], matvec[1]);
		/*#elif ( MACH==CRAY_PVP )
	        CGEMV2(&nsupr, &block_nrow, &segsze, &lusup[luptr],
		       &tri[0][kfnz-kfnz2[0]], &tri[1][kfnz-kfnz2[1]],
		       matvec[0], matvec[1]); */
#else
		cmatvec2 (nsupr, block_nrow, segsze, &lusup[luptr],
			  &tri[0][kfnz-kfnz2[0]], &tri[1][kfnz-kfnz2[1]],
			  matvec[0], matvec[1]);
#endif

#ifdef TIMING
		utime[FLOAT] += SuperLU_timer_() - f_time;
#endif	    
		/* end for two destination update */
	    } else { /* wait for a second column */
		jj2[0] = jj;
		twocols = 1;
	    }

	} /* for jj ... */

	if ( twocols == 1 ) { /* one more column left */
	    i = jj2[0] - jcol;
	    tri[0] = tempv + ldaTmp * i;
	    matvec[0] = tri[0] + maxsuper;
	    kfnz = repfnz[i*n + krep];
	    no_zeros = kfnz - fsupc;
	    segsze = krep - kfnz + 1;
	    luptr = luptr1 + nsupr * no_zeros;
	    
#ifdef USE_VENDOR_BLAS
#if ( MACH==CRAY_PVP )
	    CGEMV( ftcs2, &block_nrow, &segsze, &alpha, &lusup[luptr], 
		   &nsupr, tri[0], &incx, &beta, matvec[0], &incy );
#else
	    cgemv_( "N", &block_nrow, &segsze, &alpha, &lusup[luptr], 
		   &nsupr32, tri[0], &incx, &beta, matvec[0], &incy );
#endif
#else
	    cmatvec(nsupr, block_nrow, segsze, &lusup[luptr],
		    tri[0], matvec[0]);
#endif
	} /* if twocols == 1 */
    
	/* Scatter matvec[*] into SPA dense[*]. */
	repfnz_col = repfnz;
	dense_col = dense;
	col_marker = spa_marker;
	col_lsub = panel_lsub;
	matvec[0] = tempv + maxsuper;
	for (jj = jcol; jj < jcol + w; ++jj, repfnz_col += n, dense_col += n,
	     col_marker += n, col_lsub += n, matvec[0] += ldaTmp) {
	    kfnz = repfnz_col[krep];
	    if ( kfnz == EMPTY ) continue; /* skip zero segment */
	    segsze = krep - kfnz + 1;
	    if ( segsze <= 3 ) continue;   /* skip unrolled cases */

	    isub = isub1;
	    for (i = 0; i < block_nrow; ++i) {
		irow = lsub[isub];
                c_sub(&dense_col[irow], &dense_col[irow],
                              &matvec[0][i]); /* Scatter-add */
#ifdef SCATTER_FOUND		
		if ( col_marker[irow] != jj ) {
		    col_marker[irow] = jj;
		    col_lsub[w_lsub_end[jj-jcol]++] = irow;
		}
#endif		
		matvec[0][i] = zero;
		++isub;
	    }
	} /* for jj ... */
	
    } /* for each block row ... */

    
    /* ------------------------------------------------
       Scatter the triangular solves into SPA dense[*].
       ------------------------------------------------ */
    repfnz_col = repfnz;
    dense_col = dense;
    TriTmp = tempv;
    for (jj = 0; jj < w; ++jj, repfnz_col += n,
	 dense_col += n, TriTmp += ldaTmp) {
	kfnz = repfnz_col[krep];
	if ( kfnz == EMPTY ) continue; /* skip any zero segment */
	segsze = krep - kfnz + 1;
	if ( segsze <= 3 ) continue;   /* skip unrolled cases */
	
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
