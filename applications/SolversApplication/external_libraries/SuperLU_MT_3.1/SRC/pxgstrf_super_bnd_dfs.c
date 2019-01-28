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
pxgstrf_super_bnd_dfs(
		      const int_t  pnum, /* process number */
		      const int_t  m,    /* number of rows in the matrix */
		      const int_t  n,    /* number of columns in the matrix */
		      const int_t  jcol, /* first column of the H-supernode */
		      const int_t  w,    /* size of the H-supernode */
		      SuperMatrix *A,  /* original matrix */
		      int_t        *perm_r,   /* in */
		      int_t        *iperm_r,  /* in; inverse of perm_r */
		      int_t        *xprune,   /* in */
		      int_t        *ispruned, /* in */
		      int_t        *marker,   /* modified */
		      int_t        *parent,   /* working array */
		      int_t        *xplore,   /* working array */
		      pxgstrf_shared_t *pxgstrf_shared /* modified */	      
		      )
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
 * Performs a symbolic structure prediction on a supernode in the Householder
 * matrix H, with jcol being the leading column.
 *
 */
    GlobalLU_t *Glu = pxgstrf_shared->Glu; /* modified */	      
    register int_t krep, chperm, chrep, kchild;
    register int_t invp_rep; /* "krep" numbered in the original A */
    register int_t krow, kperm, xdfs, maxdfs, kpar;
    register int_t fsupc, k, jj, found;
    register int_t nrow;  /* union of the nonzero rows in a supernode */
    NCPformat    *Astore;
    int_t          *asub, *xa_begin, *xa_end;
    int_t          *xsup, *xsup_end, *supno, *lsub, *xlsub, *xlsub_end;

    /* Initialize pointers */
    xsup       = Glu->xsup;
    xsup_end   = Glu->xsup_end;
    supno      = Glu->supno;
    lsub       = Glu->lsub;
    xlsub      = Glu->xlsub;
    xlsub_end  = Glu->xlsub_end;
    Astore   = A->Store;
    asub     = Astore->rowind;
    xa_begin = Astore->colbeg;
    xa_end   = Astore->colend;
    nrow = 0;
    found = n + jcol;

    /* For each column in the H-supernode */
    for (jj = jcol; jj < jcol + w; ++jj) {

      /* For each nonz in A[*,jj] do dfs */
      for (k = xa_begin[jj]; k < xa_end[jj]; ++k) {
	  krow = asub[k];
	  
	  /* krow was visited before, go to the next nonzero. */
	  if ( marker[krow] == found ) continue;
	  
  	  /* For each unmarked nbr krow of jj ...   */
	  kperm = perm_r[krow];
	  
	  if ( kperm == EMPTY ) { /* krow is in L */
	      marker[krow] = found;
	      ++nrow;
	  } else {
	      /* krow is in U: if its supernode-rep krep has been explored,
		               skip the search.       */
	      krep = SUPER_REP( supno[kperm] );
	      invp_rep = iperm_r[krep];
		  
	      /* Perform dfs starting at krep */
              if ( marker[invp_rep] != found ) {
		  marker[invp_rep] = found;
		  parent[krep] = EMPTY;
		  if ( ispruned[krep] ) {
		      if ( SINGLETON( supno[krep] ) )
			  xdfs = xlsub_end[krep];
		      else xdfs = xlsub[krep];
		      maxdfs = xprune[krep];
		  } else {
		      fsupc = SUPER_FSUPC( supno[krep] );
		      xdfs = xlsub[fsupc] + krep-fsupc+1;
		      maxdfs = xlsub_end[fsupc];
		  }

		  do {
		      /* For each unmarked kchild of krep ... */
		      while ( xdfs < maxdfs ) {
			  kchild = lsub[xdfs];
			  xdfs++;
			  if (marker[kchild] != found) { /* Not reached yet */
			      chperm = perm_r[kchild];
			      
			      if ( chperm == EMPTY ) { /* kchild is in L */
				  marker[kchild] = found;
				  ++nrow;
			      } else {
				  /* kchild is in U: 
				   *   chrep = its supernode-rep. If its rep
				   *   has been explored, skip the search.
				   */
				  chrep = SUPER_REP( supno[chperm] );
				  invp_rep = iperm_r[chrep];
				  
				  /* Continue dfs at snode-rep of kchild */
				  if ( marker[invp_rep] != found ) {
				      marker[invp_rep] = found;
				      xplore[krep] = xdfs;
				      xplore[m + krep] = maxdfs;
				      parent[chrep] = krep;
				      krep = chrep;/* Go deeper down G(L^t) */
				      xdfs = xlsub[krep];     
				      maxdfs = xprune[krep];
				      if ( ispruned[krep] ) {
					  if ( SINGLETON( supno[krep] ) )
					      xdfs = xlsub_end[krep];
					  else xdfs = xlsub[krep];
					  maxdfs = xprune[krep];
				      } else {
					  fsupc = SUPER_FSUPC(supno[krep]);
					  xdfs = xlsub[fsupc] + krep-fsupc+1;
					  maxdfs = xlsub_end[fsupc];
				      }
				  } /* if */
			      } /* else */
			  } /* if... */
		      } /* while xdfs < maxdfs */
		      
		      /* krow has no more unexplored nbrs:
		       *    Place snode-rep krep in postorder dfs, if this 
		       *    segment is seen for the first time. Note that
		       *    the "repfnz[krep]" may change later.
		       *    Backtrack dfs to its parent.
		       */
		      kpar = parent[krep]; /* Pop stack, mimic recursion */
		      if ( kpar == EMPTY ) break; /* dfs done */
		      krep = kpar;
		      xdfs = xplore[krep];
		      maxdfs = xplore[m+krep];
		  } while ( kpar != EMPTY ); /* do-while - until empty stack */
	      } /* if */
	  } /* else */
      } /* for each nonz in A[*,jj] */
    } /* for jj ... */

    DynamicSetMap(pnum, jcol, nrow * w, pxgstrf_shared);
    
/*    for (i = 1; i < w; ++i) Glu->map_in_sup[jcol + i] = -i;*/
    
#if ( DEBUGlevel>=1 )
    printf("(%d) pxgstrf_super_bnd_dfs(): jcol= %d, w= %d, nrow= %d\n",
	   pnum, jcol, w, nrow);
#endif
}
