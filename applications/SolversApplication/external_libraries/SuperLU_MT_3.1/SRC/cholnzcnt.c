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

#define ROOT     (neqns)

/*
     5/20/95  X.S. Li:
         Translated from fcnthn.f using f2c
         Modified to use 0-based indexing in C
         Initialize xsup = 0 as suggested by B. Peyton

     10/17/07 X.S. Li:
	 Modified to record supernode partition of L in part_super_L[*]:
	     part_super_L[k] = size of the supernode beginning at column k;
 	                     = 0, elsewhere.
         
   
   ***********************************************************************   
   ***********************************************************************   

     Version:        0.3   
     Last modified:  January 12, 1995   
     Authors:        Esmond G. Ng and Barry W. Peyton   

     Mathematical Sciences Section, Oak Ridge National Laboratoy   

   ***********************************************************************   
   ***********************************************************************   
   **************     FCNTHN  ..... FIND NONZERO COUNTS    ***************   
   ***********************************************************************   
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
         (I) XADJ(*)     -   ARRAY OF LENGTH NEQNS+1, CONTAINING POINTERS   
                             TO THE ADJACENCY STRUCTURE.   
         (I) ADJNCY(*)   -   ARRAY OF LENGTH XADJ(NEQNS+1)-1, CONTAINING   
                             THE ADJACENCY STRUCTURE.   
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
	 (I) PART_SUPER_L   SUPERNODE PARTITION IN THE FACTORED MATRIX L.

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

int_t cholnzcnt(int_t neqns, int_t *xadj, int_t *adjncy,
	      int_t *perm, int_t *invp, int_t *etpar, 
	      int_t *colcnt, int_t *nlnz, int_t *part_super_L)
{
    /* Local variables */
    int_t temp, xsup, last1, last2, j, k, lflag, pleaf, hinbr, jstop,
	  jstrt, ifdesc, oldnbr, parent, lownbr, lca;

    /* Work arrays */
    int_t *rowcnt = intMalloc(neqns);
    int_t *set    = intMalloc(neqns);
    int_t *prvlf  = intMalloc(neqns);
    int_t *level  = intMalloc(neqns + 1);    /* length n+1 */
    int_t *weight = intMalloc(neqns + 1);    /* length n+1 */
    int_t *fdesc  = intMalloc(neqns + 1);    /* length n+1 */
    int_t *nchild = intMalloc(neqns + 1);    /* length n+1 */
    int_t *prvnbr = intMalloc(neqns);
    

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
	nchild[k] = 0;
	prvnbr[k] = EMPTY;
	part_super_L[k] = 0;
    }
    nchild[ROOT] = 0;
    fdesc[ROOT] = EMPTY;
    for (k = 0; k < neqns; ++k) {
	parent = etpar[k];
	weight[parent] = 0;
	++nchild[parent];
	ifdesc = fdesc[k];
	if (ifdesc < fdesc[parent]) {
	    fdesc[parent] = ifdesc;
	}
    }

    xsup = 0;  /* BUG FIX */
    
    /*   ------------------------------------   
         FOR EACH ``LOW NEIGHBOR'' LOWNBR ...   
         ------------------------------------ */
    for (lownbr = 0; lownbr < neqns; ++lownbr) {
	lflag = 0;
	ifdesc = fdesc[lownbr];
	oldnbr = perm[lownbr];
	jstrt = xadj[oldnbr];
	jstop = xadj[oldnbr + 1];
	/*   -----------------------------------------------   
             FOR EACH ``HIGH NEIGHBOR'', HINBR OF LOWNBR ...   
             ----------------------------------------------- */
	for (j = jstrt; j < jstop; ++j) {
	    hinbr = invp[adjncy[j]];
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
	}
	/* ----------------------------------------------------   
	   DECREMENT WEIGHT ( PARENT(LOWNBR) ).   
	   SET ( P(LOWNBR) ) <-- SET ( P(LOWNBR) ) + SET(XSUP).   
	   ---------------------------------------------------- */
	parent = etpar[lownbr];
	--weight[parent];
	if (lflag == 1 || nchild[lownbr] >= 2) {
	    /* lownbr is the first vertex of a supernode */
	    part_super_L[xsup] = lownbr - xsup;
	    xsup = lownbr;
	}
	set[xsup] = parent;
    }

    part_super_L[xsup] = neqns - xsup;
    
    /* ---------------------------------------------------------   
       USE WEIGHTS TO COMPUTE COLUMN (AND TOTAL) NONZERO COUNTS.   
       --------------------------------------------------------- */
    *nlnz = 0;
    for (k = 0; k < neqns; ++k) {
	temp = colcnt[k] + weight[k];
	colcnt[k] = temp;
	*nlnz += temp;
	parent = etpar[k];
	if (parent != 0) {
	    colcnt[parent] += temp;
	}
    }

    free (rowcnt);
    free (set);
    free (prvlf);
    free (level);
    free (weight);
    free (fdesc);
    free (nchild);
    free (prvnbr);

#if ( PRNTlevel>=1 )
    printf(".. cholnzcnt() nlnz %d\n", *nlnz);
#endif

    return 0;
} /* cholnzcnt */

