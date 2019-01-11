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

void dlsolve(int_t, int_t, double*, double*);
void dmatvec(int_t, int_t, int_t, double*, double*, double*);


int_t
pdgstrf_column_bmod(
		    const int_t  pnum,   /* process number */
		    const int_t  jcol,   /* current column in the panel */
		    const int_t  fpanelc,/* first column in the panel */
		    const int_t  nseg,   /* number of s-nodes to update jcol */
		    int_t        *segrep,/* in */
		    int_t        *repfnz,/* in */
		    double     *dense, /* modified */
		    double     *tempv, /* working array */
		    pxgstrf_shared_t *pxgstrf_shared, /* modified */
		    Gstat_t *Gstat     /* modified */
		    )
{
/*
 * -- SuperLU MT routine (version 3.0) --
 * Lawrence Berkeley National Lab, Univ. of California Berkeley,
 * and Xerox Palo Alto Research Center.
 * September 10, 2007
 *
 * Purpose:
 * ========
 *    Performs numeric block updates (sup-col) in topological order.
 *    It features: col-col, 2cols-col, 3cols-col, and sup-col updates.
 *    Special processing on the supernodal portion of L\U[*,j].
 *
 * Return value:
 * =============
 *      0 - successful return
 *    > 0 - number of bytes allocated when run out of space
 *
 */
#if ( MACH==CRAY_PVP )
    _fcd ftcs1 = _cptofcd("L", strlen("L")),
         ftcs2 = _cptofcd("N", strlen("N")),
         ftcs3 = _cptofcd("U", strlen("U"));
#endif
    
#ifdef USE_VENDOR_BLAS    
    int         incx = 1, incy = 1;
    double      alpha, beta;
#endif
    GlobalLU_t *Glu = pxgstrf_shared->Glu;   /* modified */
    
    /* krep = representative of current k-th supernode
     * fsupc = first supernodal column
     * nsupc = no of columns in supernode
     * nsupr = no of rows in supernode (used as leading dimension)
     * luptr = location of supernodal LU-block in storage
     * kfnz = first nonz in the k-th supernodal segment
     * no_zeros = no of leading zeros in a supernodal U-segment
     */
    double	  ukj, ukj1, ukj2;
    register int_t lptr, kfnz, isub, irow, i, no_zeros;
    register int_t luptr, luptr1, luptr2;
    int_t          fsupc;
    int          nsupc, nsupr, segsze;
    int          nrow; /* No of rows in the matrix of matrix-vector */
    int_t          jsupno, k, ksub, krep, krep_ind, ksupno;
    int_t          ufirst, nextlu;
    int_t          fst_col; /* First column within small LU update */
    int_t          d_fsupc; /* Distance between the first column of the current
			     panel and the first column of the current snode.*/
    int_t          *xsup, *supno;
    int_t          *lsub, *xlsub, *xlsub_end;
    double       *lusup;
    int_t          *xlusup, *xlusup_end;
    double       *tempv1;
    int_t          mem_error;
    register float flopcnt;

    double      zero = 0.0;
    double      one = 1.0;
    double      none = -1.0;

    xsup       = Glu->xsup;
    supno      = Glu->supno;
    lsub       = Glu->lsub;
    xlsub      = Glu->xlsub;
    xlsub_end  = Glu->xlsub_end;
    lusup      = Glu->lusup;
    xlusup     = Glu->xlusup;
    xlusup_end = Glu->xlusup_end;
    jsupno     = supno[jcol];

    /* 
     * For each nonz supernode segment of U[*,j] in topological order 
     */
    k = nseg - 1;
    for (ksub = 0; ksub < nseg; ksub++) {

	krep = segrep[k];
	k--;
	ksupno = supno[krep];
#if ( DEBUGlvel>=2 )
if (jcol==BADCOL)
printf("(%d) pdgstrf_column_bmod[1]: %d, nseg %d, krep %d, jsupno %d, ksupno %d\n",
       pnum, jcol, nseg, krep, jsupno, ksupno);
#endif    
	if ( jsupno != ksupno ) { /* Outside the rectangular supernode */

	    fsupc = xsup[ksupno];
	    fst_col = SUPERLU_MAX ( fsupc, fpanelc );

  	    /* Distance from the current supernode to the current panel; 
	       d_fsupc=0 if fsupc >= fpanelc. */
  	    d_fsupc = fst_col - fsupc; 

	    luptr = xlusup[fst_col] + d_fsupc;
	    lptr = xlsub[fsupc] + d_fsupc;
	    kfnz = repfnz[krep];
	    kfnz = SUPERLU_MAX ( kfnz, fpanelc );
	    segsze = krep - kfnz + 1;
	    nsupc = krep - fst_col + 1;
	    nsupr = xlsub_end[fsupc] - xlsub[fsupc]; /* Leading dimension */
	    nrow = nsupr - d_fsupc - nsupc;
	    krep_ind = lptr + nsupc - 1;

            flopcnt = segsze * (segsze - 1) + 2 * nrow * segsze;
	    Gstat->procstat[pnum].fcops += flopcnt;

#if ( DEBUGlevel>=2 )
if (jcol==BADCOL)	    
printf("(%d) pdgstrf_column_bmod[2]: %d, krep %d, kfnz %d, segsze %d, d_fsupc %d,\
fsupc %d, nsupr %d, nsupc %d\n",
       pnum, jcol, krep, kfnz, segsze, d_fsupc, fsupc, nsupr, nsupc);

#endif


	    /* 
	     * Case 1: Update U-segment of size 1 -- col-col update 
	     */
	    if ( segsze == 1 ) {
	  	ukj = dense[lsub[krep_ind]];
		luptr += nsupr*(nsupc-1) + nsupc;

		for (i = lptr + nsupc; i < xlsub_end[fsupc]; ++i) {
		    irow = lsub[i];
		    dense[irow] -=  ukj*lusup[luptr];
		    luptr++;
		}
	    } else if ( segsze <= 3 ) {
		ukj = dense[lsub[krep_ind]];
		luptr += nsupr*(nsupc-1) + nsupc-1;
		ukj1 = dense[lsub[krep_ind - 1]];
		luptr1 = luptr - nsupr;
		if ( segsze == 2 ) { /* Case 2: 2cols-col update */
		    ukj -= ukj1 * lusup[luptr1];
		    dense[lsub[krep_ind]] = ukj;
		    for (i = lptr + nsupc; i < xlsub_end[fsupc]; ++i) {
		    	irow = lsub[i];
		    	luptr++;
		    	luptr1++;
		    	dense[irow] -= ( ukj*lusup[luptr]
					+ ukj1*lusup[luptr1] );
		    }
		} else { /* Case 3: 3cols-col update */
		    ukj2 = dense[lsub[krep_ind - 2]];
		    luptr2 = luptr1 - nsupr;
		    ukj1 -= ukj2 * lusup[luptr2-1];
		    ukj = ukj - ukj1*lusup[luptr1] - ukj2*lusup[luptr2];
		    dense[lsub[krep_ind]] = ukj;
		    dense[lsub[krep_ind-1]] = ukj1;
		    for (i = lptr + nsupc; i < xlsub_end[fsupc]; ++i) {
		    	irow = lsub[i];
		    	luptr++;
		    	luptr1++;
			luptr2++;
		    	dense[irow] -= ( ukj*lusup[luptr]
			     + ukj1*lusup[luptr1] + ukj2*lusup[luptr2] );
		    }
		}


	    } else {
	  	/*
		 * Case: sup-col update
		 * Perform a triangular solve and block update,
		 * then scatter the result of sup-col update to dense
		 */
		no_zeros = kfnz - fst_col;

	        /* Copy U[*,j] segment from dense[*] to tempv[*] */
	        isub = lptr + no_zeros;
	        for (i = 0; i < segsze; i++) {
	  	    irow = lsub[isub];
		    tempv[i] = dense[irow];
		    ++isub; 
	        }

	        /* Dense triangular solve -- start effective triangle */
		luptr += nsupr * no_zeros + no_zeros; 
#ifdef USE_VENDOR_BLAS
#if ( MACH==CRAY_PVP )
		STRSV( ftcs1, ftcs2, ftcs3, &segsze, &lusup[luptr], 
		       &nsupr, tempv, &incx );
#else
		dtrsv_( "L", "N", "U", &segsze, &lusup[luptr], 
		       &nsupr, tempv, &incx );
#endif
		
 		luptr += segsze;  /* Dense matrix-vector */
		tempv1 = &tempv[segsze];
		alpha = one;
		beta = zero;
#if ( MACH==CRAY_PVP )
		SGEMV( ftcs2, &nrow, &segsze, &alpha, &lusup[luptr], 
		       &nsupr, tempv, &incx, &beta, tempv1, &incy );
#else
		dgemv_( "N", &nrow, &segsze, &alpha, &lusup[luptr], 
		       &nsupr, tempv, &incx, &beta, tempv1, &incy );
#endif
#else
		dlsolve ( (int_t) nsupr, (int_t) segsze, &lusup[luptr], tempv );

 		luptr += segsze;  /* Dense matrix-vector */
		tempv1 = &tempv[segsze];
		dmatvec ((int_t)nsupr, (int_t)nrow, (int_t)segsze, &lusup[luptr], tempv, tempv1);
#endif
                /* Scatter tempv[] into SPA dense[*] */
                isub = lptr + no_zeros;
                for (i = 0; i < segsze; i++) {
                    irow = lsub[isub];
                    dense[irow] = tempv[i]; /* Scatter */
                    tempv[i] = zero;
                    isub++;
                }

		/* Scatter tempv1[] into SPA dense[*] */
		for (i = 0; i < nrow; i++) {
		    irow = lsub[isub];
                    dense[irow] -= tempv1[i];
		    tempv1[i] = zero;
		    ++isub;
		}
	    } /* else segsze >= 4 */
	    
	} /* if jsupno ... */

    } /* for each segment... */

    
    /* ------------------------------------------
       Process the supernodal portion of L\U[*,j]
       ------------------------------------------ */
    
    fsupc = SUPER_FSUPC (jsupno);
    nsupr = xlsub_end[fsupc] - xlsub[fsupc];
    if ( (mem_error = Glu_alloc(pnum, jcol, nsupr, LUSUP, &nextlu, 
			       pxgstrf_shared)) )
	return mem_error;
    xlusup[jcol] = nextlu;
    lusup = Glu->lusup;
    
    /* Gather the nonzeros from SPA dense[*,j] into L\U[*,j] */
    for (isub = xlsub[fsupc]; isub < xlsub_end[fsupc]; ++isub) {
  	irow = lsub[isub];
	lusup[nextlu] = dense[irow];
	dense[irow] = zero;
#ifdef DEBUG
if (jcol == -1)
    printf("(%d) pdgstrf_column_bmod[lusup] jcol %d, irow %d, lusup %.10e\n",
	   pnum, jcol, irow, lusup[nextlu]);
#endif	
	++nextlu;
    }
    xlusup_end[jcol] = nextlu; /* close L\U[*,jcol] */

#if ( DEBUGlevel>=2 )
if (jcol == -1) {
    nrow = xlusup_end[jcol] - xlusup[jcol];
    print_double_vec("before sup-col update", nrow, &lsub[xlsub[fsupc]],
		     &lusup[xlusup[jcol]]);
}
#endif    
    
    /*
     * For more updates within the panel (also within the current supernode), 
     * should start from the first column of the panel, or the first column 
     * of the supernode, whichever is bigger. There are 2 cases:
     *    (1) fsupc < fpanelc,  then fst_col := fpanelc
     *    (2) fsupc >= fpanelc, then fst_col := fsupc
     */
    fst_col = SUPERLU_MAX ( fsupc, fpanelc );

    if ( fst_col < jcol ) {

  	/* distance between the current supernode and the current panel;
	   d_fsupc=0 if fsupc >= fpanelc. */
  	d_fsupc = fst_col - fsupc;

	lptr = xlsub[fsupc] + d_fsupc;
	luptr = xlusup[fst_col] + d_fsupc;
	nsupr = xlsub_end[fsupc] - xlsub[fsupc]; /* Leading dimension */
	nsupc = jcol - fst_col;	/* Excluding jcol */
	nrow = nsupr - d_fsupc - nsupc;

	/* points to the beginning of jcol in supernode L\U[*,jsupno] */
	ufirst = xlusup[jcol] + d_fsupc;	

#if ( DEBUGlevel>=2 )
if (jcol==BADCOL)
printf("(%d) pdgstrf_column_bmod[3] jcol %d, fsupc %d, nsupr %d, nsupc %d, nrow %d\n",
       pnum, jcol, fsupc, nsupr, nsupc, nrow);
#endif    

        flopcnt = nsupc * (nsupc - 1) + 2 * nrow * nsupc;
	Gstat->procstat[pnum].fcops += flopcnt;

/*	ops[TRSV] += nsupc * (nsupc - 1);
	ops[GEMV] += 2 * nrow * nsupc;    */
	
#ifdef USE_VENDOR_BLAS
	alpha = none; beta = one; /* y := beta*y + alpha*A*x */
#if ( MACH==CRAY_PVP )
	STRSV( ftcs1, ftcs2, ftcs3, &nsupc, &lusup[luptr], 
	       &nsupr, &lusup[ufirst], &incx );
	SGEMV( ftcs2, &nrow, &nsupc, &alpha, &lusup[luptr+nsupc], &nsupr,
	       &lusup[ufirst], &incx, &beta, &lusup[ufirst+nsupc], &incy );
#else
	dtrsv_( "L", "N", "U", &nsupc, &lusup[luptr], 
	       &nsupr, &lusup[ufirst], &incx );
	dgemv_( "N", &nrow, &nsupc, &alpha, &lusup[luptr+nsupc], &nsupr,
	       &lusup[ufirst], &incx, &beta, &lusup[ufirst+nsupc], &incy );
#endif
#else
	dlsolve ( nsupr, nsupc, &lusup[luptr], &lusup[ufirst] );

	dmatvec ( nsupr, nrow, nsupc, &lusup[luptr+nsupc],
		 &lusup[ufirst], tempv );
	
        /* Copy updates from tempv[*] into lusup[*] */
	isub = ufirst + nsupc;
	for (i = 0; i < nrow; i++) {
            lusup[isub] -= tempv[i];
            tempv[i] = 0.0;
	    ++isub;
	}
#endif
    } /* if fst_col < jcol ... */ 

    return 0;
}
