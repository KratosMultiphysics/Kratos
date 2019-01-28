/*! \file
Copyright (c) 2003, The Regents of the University of California, through
Lawrence Berkeley National Laboratory (subject to receipt of any required 
approvals from U.S. Dept. of Energy) 

All rights reserved. 

The source code is distributed under BSD license, see the file License.txt
at the top-level directory.
*/

#include "slu_mt_cdefs.h"

int_t
pcgstrf_snode_bmod(
		   const int_t  pnum,   /* process number */
		   const int_t  jcol,   /* in - current column in the s-node */
		   const int_t  jsupno, /* in */
		   const int_t  fsupc,  /* in - first column in the s-node */
		   complex     *dense, /* in */
		   complex     *tempv, /* working array */
		   GlobalLU_t *Glu,   /* modified */
		   Gstat_t *Gstat     /* modified */
		   )
{
/*
 * -- SuperLU MT routine (version 2.0) --
 * Lawrence Berkeley National Lab, Univ. of California Berkeley,
 * and Xerox Palo Alto Research Center.
 * September 10, 2007
 *
 * Performs numeric block updates within the relaxed supernode. 
 */

    complex      zero = {0.0, 0.0};
    complex      one = {1.0, 0.0};
    complex      none = {-1.0, 0.0};

#if ( MACH==CRAY_PVP )
    _fcd ftcs1, ftcs2, ftcs3;
#endif
#ifdef USE_VENDOR_BLAS    
    int            incx = 1, incy = 1;
    complex         alpha = none, beta = one;
#endif
    
    int_t            luptr;
    int              nsupc, nsupr, nrow;
    int_t            isub, irow, i, iptr; 
    register int_t   ufirst, nextlu;
    complex         *lusup;
    int_t            *lsub, *xlsub, *xlsub_end, *xlusup, *xlusup_end;
    register float flopcnt;

    lsub       = Glu->lsub;
    xlsub      = Glu->xlsub;
    xlsub_end  = Glu->xlsub_end;
    lusup      = Glu->lusup;
    xlusup     = Glu->xlusup;
    xlusup_end = Glu->xlusup_end;

    nextlu = xlusup[jcol];
    
    /*
     *	Process the supernodal portion of L\U[*,j]
     */
    for (isub = xlsub[fsupc]; isub < xlsub_end[fsupc]; isub++) {
  	irow = lsub[isub];
	lusup[nextlu] = dense[irow];
	dense[irow] = zero;
	++nextlu;
    }

    xlusup_end[jcol] = nextlu;
    
    if ( fsupc < jcol ) {

	luptr = xlusup[fsupc];
	nsupr = xlsub_end[fsupc] - xlsub[fsupc];
	nsupc = jcol - fsupc;	/* Excluding jcol */
	ufirst = xlusup[jcol];	/* Points to the beginning of column
				   jcol in supernode L\U(jsupno). */
	nrow = nsupr - nsupc;
	
	Gstat->procstat[pnum].fcops += flopcnt;

/*	ops[TRSV] += nsupc * (nsupc - 1);
	ops[GEMV] += 2 * nrow * nsupc;    */

#ifdef USE_VENDOR_BLAS
#if ( MACH==CRAY_PVP )
	ftcs1 = _cptofcd("L", strlen("L"));
	ftcs2 = _cptofcd("N", strlen("N"));
	ftcs3 = _cptofcd("U", strlen("U"));
	CTRSV( ftcs1, ftcs2, ftcs3, &nsupc, &lusup[luptr], &nsupr, 
	      &lusup[ufirst], &incx );
	CGEMV( ftcs2, &nrow, &nsupc, &alpha, &lusup[luptr+nsupc], &nsupr, 
	      &lusup[ufirst], &incx, &beta, &lusup[ufirst+nsupc], &incy );
#else
	ctrsv_( "L", "N", "U", &nsupc, &lusup[luptr], &nsupr, 
	      &lusup[ufirst], &incx );
	cgemv_( "N", &nrow, &nsupc, &alpha, &lusup[luptr+nsupc], &nsupr, 
		&lusup[ufirst], &incx, &beta, &lusup[ufirst+nsupc], &incy );
#endif
#else
	clsolve ( nsupr, nsupc, &lusup[luptr], &lusup[ufirst] );
	cmatvec ( nsupr, nrow, nsupc, &lusup[luptr+nsupc], 
		 &lusup[ufirst], &tempv[0] );

        /* Scatter tempv[*] into lusup[*] */
	iptr = ufirst + nsupc;
	for (i = 0; i < nrow; i++) {
            c_sub(&lusup[iptr], &lusup[iptr], &tempv[i]);
            ++iptr;
            tempv[i] = zero;
	}
#endif

    }

    return 0;
}
