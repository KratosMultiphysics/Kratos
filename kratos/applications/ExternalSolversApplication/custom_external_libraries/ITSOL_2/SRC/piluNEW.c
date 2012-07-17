#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "globheads.h"
#include "protos.h"

int pilu(p4ptr amat, csptr B, csptr C, double *droptol, 
	 int *lfil, csptr schur) {
/*---------------------------------------------------------------------- 
| PARTIAL ILUT -
| Converted to C so that dynamic memory allocation may be implememted
| in order to have no dropping in block LU factors.
|----------------------------------------------------------------------
| Partial block ILU factorization with dual truncation. 
|                                                                      
| |  B   F  |        |    L      0  |   |  U   L^{-1} F |
| |         |   =    |              | * |               |
| |  E   C  |        | E U^{-1}  I  |   |  0       S    |                   
|                                                                      
| where B is a sub-matrix of dimension B->n.
| 
|----------------------------------------------------------------------
|
| on entry:
|========== 
| ( amat ) = Permuted matrix stored in a PerMat4 struct on entry -- 
|            Individual matrices stored in SpaFmt structs.
|            On entry matrices have C (0) indexing.
|            on return contains also L and U factors.
|            Individual matrices stored in SpaFmt structs.
|            On return matrices have C (0) indexing.
|
| lfil[0]  =  number nonzeros in L-part
| lfil[1]  =  number nonzeros in U-part
| lfil[2]  =  number nonzeros in L^{-1} F
| lfil[3]  =  not used
| lfil[4]  =  number nonzeros in Schur complement
|
| droptol[0] = threshold for dropping small terms in L during
|              factorization.
| droptol[1] = threshold for dropping small terms in U.
| droptol[2] = threshold for dropping small terms in L^{-1} F during
|              factorization.
| droptol[3] = threshold for dropping small terms in E U^{-1} during
|              factorization.
| droptol[4] = threshold for dropping small terms in Schur complement
|              after factorization is completed.
|
| On return:
|===========
|
| (schur)  = contains the Schur complement matrix (S in above diagram)
|            stored in SpaFmt struct with C (0) indexing.
|
|
|       integer value returned:
|
|             0   --> successful return.
|             1   --> Error.  Input matrix may be wrong.  (The 
|                         elimination process has generated a
|                         row in L or U whose length is > n.)
|             2   --> Memory allocation error.
|             5   --> Illegal value for lfil or last.
|             6   --> zero row in B block encountered.
|             7   --> zero row in [E C] encountered.
|             8   --> zero row in new Schur complement
|----------------------------------------------------------------------- 
| work arrays:
|=============
| jw, jwrev = integer work arrays of length B->n.
| w         = real work array of length B->n. 
| jw2, jwrev2 = integer work arrays of length C->n.
| w2          = real work array of length C->n. 
|----------------------------------------------------------------------- 
|     All processing is done using C indexing.
|--------------------------------------------------------------------*/
   int i, ii, j, jj, jcol, jpos, jrow, k, *jw, *jwrev;
   int **lfja, *lflen, len, len2, lenu, lenl, rmax;
   int *jw2, *jwrev2, lsize, rsize;
   int fil0=lfil[0],fil1=lfil[1],fil2=lfil[2],fil4=lfil[4];
   double tnorm, tabs, tmax, t, s, fact, *w, *w2, **lfma;
   double drop0=droptol[0], drop1=droptol[1], drop2=droptol[2];
   double drop3=droptol[3], drop4=droptol[4];
   int lrowz, *lrowj, rrowz, *rrowj;
   double *lrowm, *rrowm;
/*-----------------------------------------------------------------------*/
   lsize = amat->nB;
   rsize = C->n;
   rmax = lsize > rsize ? lsize : rsize;
   jw = (int *) Malloc(rmax*sizeof(int), "pilu:1" );
   w = (double *) Malloc(rmax*sizeof(double), "pilu:2" );
   jwrev = (int *) Malloc(rmax*sizeof(int), "pilu:3" );
   jw2 = (int *) Malloc(rmax*sizeof(int), "pilu:4" );
   w2 = (double *) Malloc(rmax*sizeof(double), "pilu:5" );
   jwrev2 = (int *) Malloc(rmax*sizeof(int), "pilu:6" );
   if (fil0 < 0 || fil1<0 || amat->L->n<=0) goto label9995;
   lfma = (double **) Malloc(lsize*sizeof(double *), "pilu:7" ); 
   lfja = (int **) Malloc(lsize*sizeof(int *), "pilu:8" ); 
   lflen = (int *) Malloc(lsize*sizeof(int), "pilu:9" ); 
/*---------------------------------------------------------------------
|    beginning of first main loop - L, U, L^{-1}F calculations
|--------------------------------------------------------------------*/
   for (j=0; j<rmax; j++)
      jwrev[j] = -1;
   for (j=0; j<rmax; j++)
      jwrev2[j] = -1;
   for (ii=0; ii<lsize; ii++) {
     lrowj = B->ja[ii];
     lrowm = B->ma[ii];
     lrowz = B->nzcount[ii];
     rrowj = amat->F->ja[ii];
     rrowm = amat->F->ma[ii];
     rrowz = amat->F->nzcount[ii];
/*---------------------------------------------------------------------
|   check for zero row in B block
|--------------------------------------------------------------------*/
     for (k=0; k<lrowz; k++)
       if (lrowm[k] != 0.0) goto label41;
     goto label9996;
/*---------------------------------------------------------------------
|     unpack B-block in arrays w, jw, jwrev
|     WE ASSUME THERE IS A DIAGONAL ELEMENT
|--------------------------------------------------------------------*/
   label41:
     lenu = 1;
     lenl = 0;
     w[ii] = 0.0;
     jw[ii] = ii;
     jwrev[ii] = ii;
     for (j=0; j<lrowz; j++) {
       jcol = lrowj[j];
       t = lrowm[j];
       if (jcol < ii) {
	 jw[lenl] = jcol;
	 w[lenl] = t;
	 jwrev[jcol] = lenl;
	 lenl++;
       }
       else if (jcol == ii)
	 w[ii] = t;
       else {
	 jpos = ii+lenu;
	 jw[jpos] = jcol;
	 w[jpos] = t;
	 jwrev[jcol] = jpos;
	 lenu++;
       }
     }
/*---------------------------------------------------------------------
|     unpack F-block in arrays w2, jw2, jwrev2 
|     (all entries are in U portion)
|--------------------------------------------------------------------*/
     len2 = 0;
     for (j=0; j<rrowz; j++) {
       jcol = rrowj[j];
       jw2[len2] = jcol;
       w2[len2] = rrowm[j];
       jwrev2[jcol] = len2;
       len2++;
     }
/*---------------------------------------------------------------------
|     Eliminate previous rows -  
|--------------------------------------------------------------------*/
     len = 0;
     for (jj=0; jj<lenl; jj++) {
/*---------------------------------------------------------------------
|    in order to do the elimination in the correct order we must select
|    the smallest column index among jw(k), k=jj+1, ..., lenl.
|--------------------------------------------------------------------*/
       jrow = jw[jj];
       k = jj;
/*---------------------------------------------------------------------
|     determine smallest column index
|--------------------------------------------------------------------*/
       for (j=jj+1; j<lenl; j++) {
	 if (jw[j] < jrow) {
	   jrow = jw[j];
	   k = j;
	 }
       }
       if (k != jj) {    
	 /*   exchange in jw   */
	 j = jw[jj];
	 jw[jj] = jw[k];
	 jw[k] = j;
	 /*   exchange in jwrev   */
	 jwrev[jrow] = jj;
	 jwrev[j] = k;
	 /*   exchange in w   */
	 s = w[jj];
	 w[jj] = w[k];
	 w[k] = s;
       }
/*---------------------------------------------------------------------
|     zero out element in row.
|--------------------------------------------------------------------*/
       jwrev[jrow] = -1;
/*---------------------------------------------------------------------
|     get the multiplier for row to be eliminated (jrow).
|--------------------------------------------------------------------*/
       lrowm = amat->U->ma[jrow];
       fact = w[jj] * lrowm[0];
       if (fabs(fact) > drop0 ) {   /*   DROPPING IN L   */
	 lrowj = amat->U->ja[jrow];
	 lrowz = amat->U->nzcount[jrow];
	 rrowj = lfja[jrow];
	 rrowm = lfma[jrow];
	 rrowz = lflen[jrow];
/*---------------------------------------------------------------------
|     combine current row and row jrow
|--------------------------------------------------------------------*/
	 for (k=1; k<lrowz; k++) {
	   s = fact * lrowm[k];
	   j = lrowj[k];
	   jpos = jwrev[j];
/*---------------------------------------------------------------------
|     dealing with U
|--------------------------------------------------------------------*/
	   if (j >= ii) {
/*---------------------------------------------------------------------
|     this is a fill-in element
|--------------------------------------------------------------------*/
	     if (jpos == -1) {
	       if (lenu > lsize) {printf("U  row = %d\n",ii);
	       goto label9991;}
	       i = ii + lenu;
	       jw[i] = j;
	       jwrev[j] = i;
	       w[i] = - s;
	       lenu++;
	     }
/*---------------------------------------------------------------------
|     this is not a fill-in element 
|--------------------------------------------------------------------*/
	     else
	       w[jpos] -= s;
	   }
/*---------------------------------------------------------------------
|     dealing  with L
|--------------------------------------------------------------------*/
	   else {
/*---------------------------------------------------------------------
|     this is a fill-in element
|--------------------------------------------------------------------*/
	     if (jpos == -1) {
	       if (lenl > lsize) {printf("L  row = %d\n",ii);
	       goto label9991;}
	       jw[lenl] = j;
	       jwrev[j] = lenl;
	       w[lenl] = - s;
	       lenl++;
	     }
/*---------------------------------------------------------------------
|     this is not a fill-in element 
|--------------------------------------------------------------------*/
	     else
	       w[jpos] -= s;
	   }
	 }
/*---------------------------------------------------------------------
|     dealing  with  L^{-1} F
|--------------------------------------------------------------------*/
	 for (k=0; k<rrowz; k++) {
	   s = fact * rrowm[k];
	   j = rrowj[k];
	   jpos = jwrev2[j];
/*---------------------------------------------------------------------
|     this is a fill-in element
|--------------------------------------------------------------------*/
	   if (jpos == -1) {
	     jw2[len2] = j;
	     jwrev2[j] = len2;
	     w2[len2] = - s;
	     len2++;
	   }
/*---------------------------------------------------------------------
|     this is not a fill-in element
|--------------------------------------------------------------------*/
	   else
	     w2[jpos] -= s;
	 }
/*---------------------------------------------------------------------
|     store this pivot element
|--------------------------------------------------------------------*/
	 w[len] = fact;
	 jw[len]  = jrow;
	 len++;
       }
     }
/*---------------------------------------------------------------------
|     reset nonzero indicators
|--------------------------------------------------------------------*/
     for (j=0; j<len2; j++)    /*  L^{-1} F block  */
       jwrev2[jw2[j]] = -1;
     for (j=0; j<lenl; j++)    /*  L block  */
       jwrev[jw[j]] = -1;
     for (j=0; j<lenu; j++)    /*  U block  */
       jwrev[jw[ii+j]] = -1;
/*---------------------------------------------------------------------
|     done reducing this row, now store L
|--------------------------------------------------------------------*/
     lenl = len > fil0 ? fil0 : len;
     amat->L->nzcount[ii] = lenl;
     if (lenl < len) 
       qsplitC(w, jw, len, lenl);
     if (len > 0) {
       amat->L->ja[ii] = (int *) Malloc(lenl*sizeof(int), "pilu:10" ); 
       amat->L->ma[ii] = (double *) Malloc(lenl*sizeof(double), "pilu:11" ); 
       memcpy(amat->L->ja[ii], jw, lenl*sizeof(int));
       memcpy(amat->L->ma[ii], w, lenl*sizeof(double));
      }
/*---------------------------------------------------------------------
|     store the diagonal element of U
|     dropping in U if size is less than drop1 * diagonal entry
|--------------------------------------------------------------------*/
     t = w[ii];
     tnorm = fabs(t);
     len = 0;
     for (j=1; j<lenu; j++) {
       if ( fabs(w[ii+j]) > drop1*tnorm ) {
	 w[len] = w[ii+j];
	 jw[len] = jw[ii+j];
	 len++;
       }
     }
     lenu = len+1 > fil1 ? fil1 : len+1;
     amat->U->nzcount[ii] = lenu;
     jpos = lenu-1;
     if (jpos < len) 
       qsplitC(w, jw, len, jpos);
     amat->U->ma[ii] = (double *) Malloc(lenu*sizeof(double), "pilu:12" );
     amat->U->ja[ii] = (int *) Malloc(lenu*sizeof(int), "pilu:13" );
     if (t == 0.0) t=(0.0001+drop1);
     amat->U->ma[ii][0] = 1.0 / t;
     amat->U->ja[ii][0] = ii;
/*---------------------------------------------------------------------
|     copy the rest of U
|--------------------------------------------------------------------*/
     memcpy(&amat->U->ja[ii][1], jw, jpos*sizeof(int));
     memcpy(&amat->U->ma[ii][1], w, jpos*sizeof(double));
/*---------------------------------------------------------------------
|     copy  L^{-1} F
|--------------------------------------------------------------------*/
     len = 0;
     for (j=0; j<len2; j++) {
       if ( fabs(w2[j]) > drop2*tnorm ) {
	 w[len] = w2[j];
	 jw[len] = jw2[j];
	 len++;
       }
     }
     lenu = len > fil2 ? fil2 : len;
     if (lenu < len)
       qsplitC(w, jw, len, lenu);
     lflen[ii] = lenu;
     
     if (lenu > 0) {
       lfja[ii]  = (int *) Malloc(lenu*sizeof(int), "pilu:14" ); 
       lfma[ii]  = (double *) Malloc(lenu*sizeof(double), "pilu:15" ); 
       memcpy(lfma[ii], w, lenu*sizeof(double));
       memcpy(lfja[ii], jw, lenu*sizeof(int)); 
     }
   }
/*---------------------------------------------------------------------
|    beginning of second main loop   E U^{-1} and Schur complement
|--------------------------------------------------------------------*/
   for (ii=0; ii<rsize; ii++) {
     lrowj = amat->E->ja[ii];
     lrowm = amat->E->ma[ii];
     lrowz = amat->E->nzcount[ii];
     rrowj = C->ja[ii];
     rrowm = C->ma[ii];
     rrowz = C->nzcount[ii];
/*---------------------------------------------------------------------
|    determine if there is a zero row in [ E C ]
|--------------------------------------------------------------------
     for (k=0; k<lrowz; k++)
       if (lrowm[k] != 0.0) goto label42;
     for (k=0; k<rrowz; k++)
       if (rrowm[k] != 0.0) goto label42;
     goto label9997;
     label42:
*/
/*---------------------------------------------------------------------
|     unpack E in arrays w, jw, jwrev
|--------------------------------------------------------------------*/
     lenl = 0;
     for (j=0; j<lrowz; j++) {
       jcol = lrowj[j];
       jw[lenl] = jcol;
       w[lenl] = lrowm[j];
       jwrev[jcol] = lenl;
       lenl++;
     }
/*---------------------------------------------------------------------
|     unpack C in arrays w2, jw2, jwrev2    
|--------------------------------------------------------------------*/
     lenu = 0;
     for (j=0; j<rrowz; j++) {
       jcol = rrowj[j];
       jw2[lenu] = jcol;
       w2[lenu] = rrowm[j];
       jwrev2[jcol] = lenu;
       lenu++;
     }
/*---------------------------------------------------------------------
|     eliminate previous rows
|--------------------------------------------------------------------*/
     len = 0;
     for (jj=0; jj<lenl; jj++) {
/*---------------------------------------------------------------------
|    in order to do the elimination in the correct order we must select
|    the smallest column index among jw(k), k=jj+1, ..., lenl.
|--------------------------------------------------------------------*/
       jrow = jw[jj];
       k = jj;
/*---------------------------------------------------------------------
|     determine smallest column index
|--------------------------------------------------------------------*/
       for (j=jj+1; j<lenl; j++) {
	 if (jw[j] < jrow) {
	   jrow = jw[j];
	   k = j;
	 }
       }
       if (k != jj) {    
	 /*   exchange in jw   */
	 j = jw[jj];
	 jw[jj] = jw[k];
	 jw[k] = j;
	 /*   exchange in jwrev   */
	 jwrev[jrow] = jj;
	 jwrev[j] = k;
	 /*   exchange in w   */
	 s = w[jj];
	 w[jj] = w[k];
	 w[k] = s;
       }
/*---------------------------------------------------------------------
|     zero out element in row.
|--------------------------------------------------------------------*/
       jwrev[jrow] = -1;
/*---------------------------------------------------------------------
|     get the multiplier for row to be eliminated (jrow).
|--------------------------------------------------------------------*/
       lrowm = amat->U->ma[jrow];
       fact = w[jj] * lrowm[0];
       if ( fabs(fact) > drop3 ) {      /*  DROPPING IN E U^{-1}   */
	 lrowj = amat->U->ja[jrow];
	 lrowz = amat->U->nzcount[jrow];
	 rrowj = lfja[jrow];
	 rrowm = lfma[jrow];
	 rrowz = lflen[jrow];
/*---------------------------------------------------------------------
|     combine current row and row jrow   -   first  E U^{-1}
|--------------------------------------------------------------------*/
	 for (k=1; k<lrowz; k++) {
	   s = fact * lrowm[k];
	   j = lrowj[k];
	   jpos = jwrev[j];
/*---------------------------------------------------------------------
|     fill-in element
|--------------------------------------------------------------------*/
	   if (jpos == -1) {
	     if (lenl > lsize) {printf(" E U^{-1}  row = %d\n",ii);
	     goto label9991;}
	     jw[lenl] = j;
	     jwrev[j] = lenl;
	     w[lenl] = - s;
	     lenl++;
	   }
/*---------------------------------------------------------------------
|     this is not a fill-in element 
|--------------------------------------------------------------------*/
	   else
	     w[jpos] -= s;
	 }
/*---------------------------------------------------------------------
|     incorporate into Schur complement   C - (E U^{-1}) (L^{-1} F)
|--------------------------------------------------------------------*/
	 for (k=0; k<rrowz; k++) {
	   s = fact * rrowm[k];
	   j = rrowj[k];
	   jpos = jwrev2[j];
/*---------------------------------------------------------------------
|     this is not a fill-in element 
|--------------------------------------------------------------------*/
	   if (jpos == -1) {
	     jw2[lenu] = j;
	     jwrev2[j] = lenu;
	     w2[lenu] = - s;
	     lenu++;
	   }
/*---------------------------------------------------------------------
|     this is not a fill-in element
|--------------------------------------------------------------------*/
	   else
	     w2[jpos] -= s;
	 }
/*---------------------------------------------------------------------
|     store this pivot element
|--------------------------------------------------------------------*/
	 w[len] = fact;
	 jw[len] = jrow;
	 len++;
       }
     }
/*---------------------------------------------------------------------
|     reset nonzero indicators
|--------------------------------------------------------------------*/
     for (j=0; j<lenu; j++)    /*  Schur complement  */
       jwrev2[jw2[j]] = -1;
     for (j=0; j<lenl; j++)    /*  E U^{-1} block  */
       jwrev[jw[j]] = -1;
/*---------------------------------------------------------------------
|     done reducing this row, now throw away row of E U^{-1}
|     and apply a dropping strategy to the Schur complement.
|
|     store the diagonal element of Schur Complement
|
|     drop in Schur complement if size less than drop4*tnorm
|     where tnorm is the size of the maximum entry in the row
|--------------------------------------------------------------------*/
     tnorm = 0.0; 
     tmax  = 0.0;
     for (j=0; j<lenu; j++) {
       tabs = fabs(w2[j]) ;
       if (tmax < tabs) tmax = tabs;
       tnorm += tabs;
     }
       /* if (fabs(w2[j]) > tnorm) tnorm =  fabs(w2[j]); */
     if (tnorm == 0.0) {
       len = 1;
       w[0] = 1.0; 
       jw[0] = ii;
     } 
     else {
       len = 0;
       /*     tabs = drop4*tmax*(tmax/tnorm); */
       tabs = drop4*tmax*tmax/( tnorm * (double) lenu);
       for (j=0; j<lenu; j++) {
	 if (fabs(w2[j]) > tabs) {
	   w[len] = w2[j];
	   jw[len] = jw2[j];
	   len++;
	 }
       }
     }
     lenu = len > fil4 ? fil4 : len;
     schur->nzcount[ii] = lenu;
     jpos = lenu;
     if (jpos < len)
       qsplitC(w, jw, len, jpos);
     schur->ma[ii] = (double *) Malloc(lenu*sizeof(double), "pilu:16" );
     schur->ja[ii] = (int *) Malloc(lenu*sizeof(int), "pilu:17" );
/*---------------------------------------------------------------------
|     copy ---
|--------------------------------------------------------------------*/
     memcpy(&schur->ja[ii][0], jw, jpos*sizeof(int));
     memcpy(&schur->ma[ii][0], w, jpos*sizeof(double));
   }
/*---------------------------------------------------------------------
|     end main loop - now do cleanup
|--------------------------------------------------------------------*/
   free(jw);
   free(w);
   free(jwrev);
   free(jw2);
   free(w2);
   free(jwrev2);
   for (i=0; i<lsize; i++) {
     if (lflen[i] > 0) {
       free(lfma[i]);
       free(lfja[i]);
     }
   }
   free(lfma);
   free(lfja);
   free(lflen);
/*---------------------------------------------------------------------
|     done  --  correct return
|--------------------------------------------------------------------*/
   return 0;
label9991:
/*  Incomprehensible error. Matrix must be wrong.  */
   return 1;
/* label9992:
   Memory allocation error.  
   return 2; */
label9995:
/*  illegal value for lfil or last entered  */
   return 5;
label9996:
/*  zero row encountered  */
/* printf("row = %d\n",ii+1); */
   return 6;
}
/*---------------------------------------------------------------------
|     end of pilut
|--------------------------------------------------------------------*/

