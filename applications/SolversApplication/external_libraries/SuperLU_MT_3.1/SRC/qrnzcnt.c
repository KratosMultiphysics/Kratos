/*! \file
Copyright (c) 2003, The Regents of the University of California, through
Lawrence Berkeley National Laboratory (subject to receipt of any required 
approvals from U.S. Dept. of Energy) 

All rights reserved. 

The source code is distributed under BSD license, see the file License.txt
at the top-level directory.
*/
#include <stdlib.h>
#include "slu_mt_ddefs.h"

#define EMPTY    (-1)
#define ROOT     (neqns)   /* dummy root of the e-tree */


int_t
qrnzcnt(int_t neqns, int_t adjlen, int_t *xadj, int_t *adjncy, int_t *zfdperm,
	int_t *perm, int_t *invp, int_t *etpar, int_t *colcnt_h,
	int_t *nlnz, int_t *part_super_ata, int_t *part_super_h)
{
/*
     o 5/20/95 Xiaoye S. Li:
         Translated from fcnthn.f using f2c;
         Modified to use 0-based indexing in C;
         Initialize xsup = 0 as suggested by B. Peyton to handle singletons.

     o 5/24/95 Xiaoye S. Li:
         Modified to compute row/column counts of R in QR factorization
           1. Compute row counts of A, and f(i) in a separate pass
	                        def
	   2. Re-define hadj[k] ===   U    { j | j in Struct(A_i*), j>k}
	                          i:f(i)==k
	 Record supernode partition in part_super_ata[*] of size neqns:
	   part_super_ata[k] = size of the supernode beginning at column k;
 	                     = 0, elsewhere.

     o 1/16/96 Xiaoye S. Li:
         Modified to incorporate row/column counts of the Householder
	 Matrix H in the QR factorization A --> H , R.
	 
	 Record supernode partition in part_super_h[*] of size neqns:
	   part_super_h[k] = size of the supernode beginning at column k;
 	                   = 0, elsewhere.
	 
   ***********************************************************************   
     Version:        0.3   
     Last modified:  January 12, 1995   
     Authors:        Esmond G. Ng and Barry W. Peyton   

     Mathematical Sciences Section, Oak Ridge National Laboratoy   

   ***********************************************************************   
   **************     FCNTHN  ..... FIND NONZERO COUNTS    ***************   
   ***********************************************************************   

     PURPOSE:   
         THIS SUBROUTINE DETERMINES THE ROW COUNTS AND COLUMN COUNTS IN   
         THE CHOLESKY FACTOR.  IT USES A DISJOINT SET UNION ALGORITHM.   

         TECHNIQUES:   
         1) SUPERNODE DETECTION.   
         2) PATH HALVING.   
         3) NO UNION BY RANK.   

     NOTES:   
         1) ASSUMES A POSTORDERING OF THE ELIMINATION TREE.   

     INPUT PARAMETERS:   
         (I) NEQNS       -   NUMBER OF EQUATIONS.   
         (I) ADJLEN      -   LENGTH OF ADJACENCY STRUCTURE.   
         (I) XADJ(*)     -   ARRAY OF LENGTH NEQNS+1, CONTAINING POINTERS   
                             TO THE ADJACENCY STRUCTURE.
         (I) ADJNCY(*)   -   ARRAY OF LENGTH ADJLEN, CONTAINING   
                             THE ADJACENCY STRUCTURE.
         (I) ZFDPERM(*)  -   THE ROW PERMUTATION VECTOR THAT PERMUTES THE
	                     MATRIX TO HAVE ZERO-FREE DIAGONAL.
			     ZFDPERM(I) = J MEANS ROW I OF THE ORIGINAL
			     MATRIX IS IN ROW J OF THE PERMUTED MATRIX.
         (I) PERM(*)     -   ARRAY OF LENGTH NEQNS, CONTAINING THE   
                             POSTORDERING.   
         (I) INVP(*)     -   ARRAY OF LENGTH NEQNS, CONTAINING THE   
                             INVERSE OF THE POSTORDERING.   
         (I) ETPAR(*)    -   ARRAY OF LENGTH NEQNS, CONTAINING THE   
                             ELIMINATION TREE OF THE POSTORDERED MATRIX.   

     OUTPUT PARAMETERS:   
         (I) ROWCNT(*)   -   ARRAY OF LENGTH NEQNS, CONTAINING THE NUMBER   
                             OF NONZEROS IN EACH ROW OF THE FACTOR,   
                             INCLUDING THE DIAGONAL ENTRY.   
         (I) COLCNT(*)   -   ARRAY OF LENGTH NEQNS, CONTAINING THE NUMBER   
                             OF NONZEROS IN EACH COLUMN OF THE FACTOR,   
                             INCLUDING THE DIAGONAL ENTRY.   
         (I) NLNZ        -   NUMBER OF NONZEROS IN THE FACTOR, INCLUDING   
                             THE DIAGONAL ENTRIES.
         (I) PART_SUPER_ATA  SUPERNODE PARTITION IN THE CHOLESKY FACTOR
	                     OF A'A.
	 (I) PART_SUPER_H    SUPERNODE PARTITION IN THE HOUSEHOLDER
	                     MATRIX H.

     WORK PARAMETERS:   
         (I) SET(*)      -   ARRAY OF LENGTH NEQNS USED TO MAINTAIN THE   
                             DISJOINT SETS (I.E., SUBTREES).   
         (I) PRVLF(*)    -   ARRAY OF LENGTH NEQNS USED TO RECORD THE   
                             PREVIOUS LEAF OF EACH ROW SUBTREE.   
         (I) LEVEL(*)    -   ARRAY OF LENGTH NEQNS+1 CONTAINING THE LEVEL   
                             (DISTANCE FROM THE ROOT).   
         (I) WEIGHT(*)   -   ARRAY OF LENGTH NEQNS+1 CONTAINING WEIGHTS   
                             USED TO COMPUTE COLUMN COUNTS.   
         (I) FDESC(*)    -   ARRAY OF LENGTH NEQNS+1 CONTAINING THE   
                             FIRST (I.E., LOWEST-NUMBERED) DESCENDANT.   
         (I) NCHILD(*)   -   ARRAY OF LENGTH NEQNS+1 CONTAINING THE   
                             NUMBER OF CHILDREN.   
         (I) PRVNBR(*)   -   ARRAY OF LENGTH NEQNS USED TO RECORD THE   
                             PREVIOUS ``LOWER NEIGHBOR'' OF EACH NODE.   

     FIRST CREATED ON    APRIL 12, 1990.   
     LAST UPDATED ON     JANUARY 12, 1995.   

   ***********************************************************************   
*/

    /* Local variables */
    int_t  temp, last1, last2, i, j, k, lflag, pleaf, hinbr, jstop,
	 jstrt, ifdesc, oldnbr, parent, lownbr, lca;
    int_t  xsup;        /* the ongoing supernode */
    int_t  *set, *prvlf, *level, *weight, *fdesc, *nchild, *prvnbr;
    int_t  *fnz;        /* first nonzero column subscript in each row */
    int_t  *marker;     /* used to remove duplicate indices */
    int_t  *fnz_hadj;   /* higher-numbered neighbors of the first nonzero
			 (higher adjacency set of A'A) */
    int_t  *hadj_begin; /* pointers to the fnz_hadj[] structure */
    int_t  *hadj_end;   /* pointers to the fnz_hadj[] structure */

    /* Locally malloc'd room for QR purpose */

    /* ----------------------------------------------------------
       FIRST set is defined as first[j] := { i : f[i] = j } ,
       which is a collection of disjoint sets of integers between
       0 and n-1.
       ---------------------------------------------------------- */    
    int_t  *first;    /* header pointing to FIRST set */
    int_t  *firstset; /* linked list to describe FIRST set */
    int_t  *weight_h; /* weights for H */
    int_t  *rowcnt;   /* row colunts for Lc */ 
    int_t  *colcnt;   /* column colunts for Lc */ 
    int_t  *rowcnt_h; /* row colunts for H */ 
    int_t  nsuper;    /* total number of fundamental supernodes in Lc */
    int_t  nhnz;
    
    set    = intMalloc(neqns);
    prvlf  = intMalloc(neqns);
    level  = intMalloc(neqns + 1);    /* length n+1 */
    weight = intMalloc(neqns + 1);    /* length n+1 */
    fdesc  = intMalloc(neqns + 1);    /* length n+1 */
    nchild = intMalloc(neqns + 1);    /* length n+1 */
    prvnbr = intMalloc(neqns);
    fnz_hadj   = intMalloc(adjlen + 2*neqns + 1);
    hadj_begin = fnz_hadj + adjlen;        /* neqns+1 */
    hadj_end   = hadj_begin + neqns + 1;   /* neqns */
    fnz        = set;    /* aliasing for the time being */
    marker     = prvlf;  /*     "    "    "             */
    
    first    = intMalloc(neqns);
    firstset = intMalloc(neqns);
    weight_h = intCalloc(neqns + 1);  /* length n+1 */
    rowcnt_h = intMalloc(neqns);
    rowcnt   = intMalloc(neqns);
    colcnt   = intMalloc(neqns);
    
    /* -------------------------------------------------------
     * Compute fnz[*], first[*], nchild[*] and row counts of A.
     * Also find supernodes in H.
     *
     * Note that the structure of each row of H forms a simple path in
     * the etree between fnz[i] and i (George, Liu & Ng (1988)).
     * The "first vertices" of the supernodes in H are characterized
     * by the following conditions:
     *     1) first nonzero in each row of A, i.e., fnz(i);
     *  or 2) nchild >= 2;
     * ------------------------------------------------------- */
    for (k = 0; k < neqns; ++k) {
	fnz[k] = first[k] = marker[k] = EMPTY;
	rowcnt[k] = part_super_ata[k] = 0;
	part_super_h[k] = 0;
	nchild[k] = 0;
    }
    nchild[ROOT] = 0;
    xsup = 0;
    for (k = 0; k < neqns; ++k) {
	parent = etpar[k];
	++nchild[parent];
	if ( k != 0 && nchild[k] >= 2 ) {
	    part_super_h[xsup] = k - xsup;
	    xsup = k;
	}
	oldnbr = perm[k];
	for (j = xadj[oldnbr]; j < xadj[oldnbr+1]; ++j) {
	    /*
	     * Renumber vertices of G(A) by postorder
	     */
/*	    i = invp[zfdperm[adjncy[j]]];*/
	    i = zfdperm[adjncy[j]];
	    ++rowcnt[i];
	    if (fnz[i] == EMPTY) {
		/*
		 * Build linked list to describe FIRST sets
		 */
		fnz[i] = k;
		firstset[i] = first[k];
		first[k] = i;
		if ( k != 0 && xsup != k ) {
		    part_super_h[xsup] = k - xsup;
		    xsup = k;
		}
	    }
	}
    }
    part_super_h[xsup] = neqns - xsup;

#ifdef CHK_NZCNT
    printf("%8s%8s%8s\n", "k", "fnz", "first");
    for (k = 0; k < neqns; ++k)
	printf("%8d%8d%8d\n", k, fnz[k], first[k]);
#endif
    
    /* Set up fnz_hadj[*] structure. */
    hadj_begin[0] = 0;
    for (k = 0; k < neqns; ++k) {
	temp = 0;
	oldnbr = perm[k];
	hadj_end[k] = hadj_begin[k];
	for (j = xadj[oldnbr]; j < xadj[oldnbr+1]; ++j) {
/*	    hinbr = invp[zfdperm[adjncy[j]]];*/
	    hinbr = zfdperm[adjncy[j]];
	    jstrt = fnz[hinbr];    /* first nonzero must be <= k */
	    if ( jstrt != k && marker[jstrt] < k ) {
		/* ----------------------------------
		   filtering k itself and duplicates
		   ---------------------------------- */
		fnz_hadj[hadj_end[jstrt]] = k;
		++hadj_end[jstrt];
		marker[jstrt] = k;
	    }
	    if ( jstrt == k ) temp += rowcnt[hinbr];
	}
	hadj_begin[k+1] = hadj_begin[k] + temp;
    }

#ifdef CHK_NZCNT
    printf("%8s%8s\n", "k", "hadj");
    for (k = 0; k < neqns; ++k) {
	printf("%8d", k);
	for (j = hadj_begin[k]; j < hadj_end[k]; ++j)
	    printf("%8d", fnz_hadj[j]);
	printf("\n");
    }
#endif
	
    /*   --------------------------------------------------   
         COMPUTE LEVEL(*), FDESC(*), NCHILD(*).   
         INITIALIZE ROWCNT(*), COLCNT(*),   
                    SET(*), PRVLF(*), WEIGHT(*), PRVNBR(*).   
         --------------------------------------------------   */
    level[ROOT] = 0;
    for (k = neqns-1; k >= 0; --k) {
	rowcnt[k] = 1;
	colcnt[k] = 0;
	set[k] = k;
	prvlf[k] = EMPTY;
	level[k] = level[etpar[k]] + 1;
	weight[k] = 1;
	fdesc[k] = k;
	prvnbr[k] = EMPTY;
    }
    fdesc[ROOT] = EMPTY;
    for (k = 0; k < neqns; ++k) {
	parent = etpar[k];
	weight[parent] = 0;
	colcnt_h[k] = 0;
	ifdesc = fdesc[k];
	if (ifdesc < fdesc[parent]) {
	    fdesc[parent] = ifdesc;
	}
    }

    xsup    = 0;      /* BUG FIX */
    nsuper = 0;
    
    /*   ------------------------------------   
         FOR EACH ``LOW NEIGHBOR'' LOWNBR ...   
         ------------------------------------ */
    for (lownbr = 0; lownbr < neqns; ++lownbr) {
	for (i = first[lownbr]; i != EMPTY; i = firstset[i]) {
	    rowcnt_h[i] = 1 + ( level[lownbr] - level[i] );
	    ++weight_h[lownbr];
	    parent = etpar[i];
	    --weight_h[parent];
	}
	
	lflag  = 0;
	ifdesc = fdesc[lownbr];
	jstrt  = hadj_begin[lownbr];
	jstop  = hadj_end[lownbr];
	/*   -----------------------------------------------   
             FOR EACH ``HIGH NEIGHBOR'', HINBR OF LOWNBR ...   
             ----------------------------------------------- */
	for (j = jstrt; j < jstop; ++j) {	    
	    hinbr = fnz_hadj[j];
	    if (hinbr > lownbr) {
                if (ifdesc > prvnbr[hinbr]) {
		    /*  -------------------------   
			INCREMENT WEIGHT(LOWNBR).   
			------------------------- */
		    ++weight[lownbr];
		    pleaf = prvlf[hinbr];
		    /*  -----------------------------------------   
			IF HINBR HAS NO PREVIOUS ``LOW NEIGHBOR'' THEN ...   
			----------------------------------------- */
		    if (pleaf == EMPTY) {
			/* -----------------------------------------   
			   ... ACCUMULATE LOWNBR-->HINBR PATH LENGTH   
			       IN ROWCNT(HINBR).   
			   ----------------------------------------- */
			rowcnt[hinbr] = rowcnt[hinbr] +
			                level[lownbr] - level[hinbr];
		    } else {
			/* -----------------------------------------   
			   ... OTHERWISE, LCA <-- FIND(PLEAF), WHICH   
                               IS THE LEAST COMMON ANCESTOR OF PLEAF   
                               AND LOWNBR. (PATH HALVING.)   
			   ----------------------------------------- */
			last1 = pleaf;
			last2 = set[last1];
			lca = set[last2];
			while ( lca != last2 ) {
			    set[last1] = lca;
			    last1 = lca;
			    last2 = set[last1];
			    lca = set[last2];
			}
			/* -------------------------------------   
			   ACCUMULATE PLEAF-->LCA PATH LENGTH IN   
			   ROWCNT(HINBR). DECREMENT WEIGHT(LCA).   
			   ------------------------------------- */
			rowcnt[hinbr] = rowcnt[hinbr] + 
					level[lownbr] - level[lca];
			--weight[lca];
		    }
		    /* ----------------------------------------------   
		       LOWNBR NOW BECOMES ``PREVIOUS LEAF'' OF HINBR.   
		       ---------------------------------------------- */
		    prvlf[hinbr] = lownbr;
		    lflag = 1;
        	}
		/* --------------------------------------------------   
		   LOWNBR NOW BECOMES ``PREVIOUS NEIGHBOR'' OF HINBR.   
		   -------------------------------------------------- */
		prvnbr[hinbr] = lownbr;
	    }
	} /* for j ... */
	
	/* ----------------------------------------------------   
	   DECREMENT WEIGHT ( PARENT(LOWNBR) ).   
	   SET ( P(LOWNBR) ) <-- SET ( P(LOWNBR) ) + SET(XSUP).   
	   ---------------------------------------------------- */
	parent = etpar[lownbr];
	--weight[parent];
	if (lflag == 1 || nchild[lownbr] >= 2) {
	    /* lownbr is detected as the beginning of the new supernode */
	    if ( lownbr != 0 ) part_super_ata[xsup] = lownbr - xsup;
	    ++nsuper;
	    xsup = lownbr;
	} else {
	    if ( parent == ROOT && ifdesc == lownbr ) {
		/* lownbr is a singleton, and begins a new supernode
		   but is not detected as doing so -- BUG FIX */
		part_super_ata[lownbr] = 1;
		++nsuper;
		xsup = lownbr;
	    }
	}
	set[xsup] = parent;
    } /* for lownbr ... */
    
    /* ---------------------------------------------------------   
       USE WEIGHTS TO COMPUTE COLUMN (AND TOTAL) NONZERO COUNTS.   
       --------------------------------------------------------- */
    *nlnz = nhnz = 0;
    for (k = 0; k < neqns; ++k) {
	/* for R */
	temp = colcnt[k] + weight[k];
	colcnt[k] = temp;
	*nlnz += temp;
	parent = etpar[k];
	if (parent != ROOT) {
	    colcnt[parent] += temp;
	}

	/* for H */
	temp = colcnt_h[k] + weight_h[k];
	colcnt_h[k] = temp;
	nhnz += temp;
	if (parent != ROOT) {
	    colcnt_h[parent] += temp;	    
	}
    }
    part_super_ata[xsup] = neqns - xsup;

    /* Fix the supernode partition in H. */
    
    free (set);
    free (prvlf);
    free (level);
    free (weight);
    free (fdesc);
    free (nchild);
    free (prvnbr);
    free (fnz_hadj);

    free (first);
    free (firstset);
    free (weight_h);
    free (rowcnt_h);
    free (rowcnt);
    free (colcnt);
    
#if ( PRNTlevel==1 )
    printf(".. qrnzcnt() nlnz %d, nhnz %d, nlnz/nhnz %.2f\n", 
		*nlnz, nhnz, (float) *nlnz/nhnz);
#endif

#if ( DEBUGlevel>=2 )
    print_int_vec("part_super_h", neqns, part_super_h);
#endif    

    return 0;
    
} /* qrnzcnt_ */

