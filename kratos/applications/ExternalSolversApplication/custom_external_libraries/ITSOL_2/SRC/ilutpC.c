#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "globheads.h"
#include "protos.h"

#ifdef ILUTIME
void msg_timer_clear(int *);
void msg_timer_start(int *);
void msg_timer_stop(int *); 
double msg_timer(int *);
#endif

int ilutpC(csptr amat, double *droptol, int *lfil, double permtol,
		int mband, ilutptr ilusch)
{
/*---------------------------------------------------------------------- 
| ILUTP -- ILUT with column pivoting -- adapted from ILUTP [Sparskit] 
| Converted to C so that dynamic memory allocation may be implememted.
| All indexing is in C format.
|----------------------------------------------------------------------
| ILUT factorization with dual truncation. 
|----------------------------------------------------------------------
|
| on entry:
|========== 
| ( amat ) =  Matrix stored in SpaFmt struct.
|
| lfil[5]  =  number nonzeros in L-part
| lfil[6]  =  number nonzeros in U-part     ( lfil >= 0 )
|
| droptol[5] = threshold for dropping small terms in L during
|              factorization.
| droptol[6] = threshold for dropping small terms in U.
|
| permtol  =  tolerance ratio used to  determine whether or not to permute
|             two columns.  At step i columns i and j are permuted when
|
|                     abs(a(i,j))*permtol > abs(a(i,i))
|
|           [0 --> never permute; good values 0.1 to 0.01]
|
| mband    =  permuting is done within a band extending to mband
|	      diagonals only. 
|             mband = 0 --> no pivoting. 
|	      mband = n --> pivot is searched in whole column 
|
|
| On return:
|===========
|
| (ilusch) =  Contains L and U factors in an LUfact struct.
|             Individual matrices stored in SpaFmt structs.
|             On return matrices have C (0) indexing.
|
| iperm    =  reverse permutation array.
|
|       integer value returned:
|
|             0   --> successful return.
|             1   --> Error.  Input matrix may be wrong.  (The 
|                         elimination process has generated a
|                         row in L or U whose length is > n.)
|             2   --> Memory allocation error.
|             5   --> Illegal value for lfil.
|             6   --> zero row encountered.
|----------------------------------------------------------------------- 
| work arrays:
|=============
| jw, jwrev = integer work arrays of length n
| w         = real work array of length n. 
|----------------------------------------------------------------------- 
|     All processing is done using C indexing.
|--------------------------------------------------------------------*/
   int i, ii, j, jj, jcol, jpos, jrow, k, *jw=NULL, *jwrev=NULL;
   int len, lenu, lenl, rmax, *iprev=NULL, fil5=lfil[5], fil6=lfil[6];
   double tnorm, t, s, fact, *w=NULL, drop5=droptol[5], drop6=droptol[6];
   int rowz, *rowj, imax, icut, *iperm=ilusch->perm2;
   double *rowm, xmax, xmax0, tmp;
   int dec, decnt=0;

#ifdef ILUTIME
   int tt, t22=22, t23=23, t24=24, t25=25, t26=26, t27=27, t28=28,
       t29=29, t30=30, t31=31, t32=32, t33=33, t34=34, t35=35;
   for (tt=22; tt<36; tt++)
      msg_timer_clear(&tt);
#endif
   ilusch->n = rmax = amat->n;
   if (rmax == 0) return(0); 
   if (rmax > 0) {
      jw = (int *) Malloc(rmax*sizeof(int), "ilutpC:1" );
      w = (double *) Malloc(rmax*sizeof(double), "ilutpC:2" );
      jwrev = (int *) Malloc(rmax*sizeof(int), "ilutpC:3" );
      iprev = (int *) Malloc(rmax*sizeof(int), "ilutpC:4" );
   }
   if (fil5<0 || fil6<0 || rmax<=0) goto label998;
/*---------------------------------------------------------------------
|    beginning of main loop - L, U calculations
|--------------------------------------------------------------------*/
   for (j=0; j<rmax; j++) {
      jwrev[j] = -1;
      iprev[j] = j;
   }
   for (ii=0; ii<rmax; ii++) {
#ifdef ILUTIME
msg_timer_start(&t22);
#endif
      rowj = amat->ja[ii];
      rowm = amat->ma[ii];
      rowz = amat->nzcount[ii];
      tnorm = 0.0;
      for (k=0; k<rowz; k++)
	 tnorm += fabs(rowm[k]);
      if (tnorm == 0.0) goto label9991;
      tnorm = tnorm / ((double) rowz);
#ifdef ILUTIME
msg_timer_stop(&t22);
msg_timer_start(&t23);
#endif
/*---------------------------------------------------------------------
|     unpack amat in arrays w, jw, jwrev
|     WE ASSUME THERE IS A DIAGONAL ELEMENT
|--------------------------------------------------------------------*/
      lenu = 1;
      lenl = 0;
      w[ii] = 0.0;
      jw[ii] = ii;
      jwrev[ii] = ii;
      for (j=0; j<rowz; j++) {
	 jcol = iprev[rowj[j]];
	 t = rowm[j];
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
#ifdef ILUTIME
msg_timer_stop(&t23);
#endif
/*---------------------------------------------------------------------
|     eliminate previous rows
|--------------------------------------------------------------------*/
      len = 0;
      for (jj=0; jj<lenl; jj++) {
/*---------------------------------------------------------------------
|    in order to do the elimination in the correct order we must select
|    the smallest column index among jw(k), k=jj+1, ..., lenl.
|--------------------------------------------------------------------*/
#ifdef ILUTIME
msg_timer_start(&t24);
#endif
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
/*
#ifdef ILUTIME
msg_timer_stop(&t24);
#endif
*/
/*---------------------------------------------------------------------
|     get the multiplier for row to be eliminated (jrow).
|--------------------------------------------------------------------*/
	 rowm = ilusch->U->ma[jrow];
	 fact = w[jj] * rowm[0];
	 if ( fabs(fact) > drop5 ) {    /*  DROPPING IN L  */
#ifdef ILUTIME
msg_timer_start(&t25);
#endif
	    rowj = ilusch->U->ja[jrow];
	    rowz = ilusch->U->nzcount[jrow];
/*---------------------------------------------------------------------
|     combine current row and row jrow
|--------------------------------------------------------------------*/
	    for (k=1; k<rowz; k++) {
	       s = fact * rowm[k];
	       j = iprev[rowj[k]];   /*  new column number  */
	       jpos = jwrev[j];
/*---------------------------------------------------------------------
|     dealing with U
|--------------------------------------------------------------------*/
	       if (j >= ii) {
/*---------------------------------------------------------------------
|     this is a fill-in element
|--------------------------------------------------------------------*/
		  if (jpos == -1) {
		     if (lenu > rmax) goto label994;
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
		     if (lenl > rmax) goto label994;
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
|     store this pivot element
|--------------------------------------------------------------------*/
#ifdef ILUTIME
msg_timer_stop(&t25);
#endif
	    w[len] = fact;
	    jw[len]  = jrow;
	    len++;
	 }
      }
/*---------------------------------------------------------------------
|     reset nonzero indicators
|--------------------------------------------------------------------*/
#ifdef ILUTIME
msg_timer_start(&t26);
#endif
      for (j=0; j<lenu; j++)    /*  U block  */
	 jwrev[jw[ii+j]] = -1;
#ifdef ILUTIME
msg_timer_stop(&t26);
#endif
/*---------------------------------------------------------------------
|     done reducing this row, now store L
|--------------------------------------------------------------------*/
#ifdef ILUTIME
msg_timer_start(&t27);
#endif
      lenl = len;
      len = lenl > fil5 ? fil5 : lenl;
      ilusch->L->nzcount[ii] = len;
      if (lenl > len)
         qsplitC(w, jw, lenl, len);
/*
  printf("  row %d   length of L = %d",ii,len);
*/
      if (len > 0) {
	 ilusch->L->ja[ii] = (int *) Malloc(len*sizeof(int), "ilutpC:5" ); 
	 ilusch->L->ma[ii] = (double *) Malloc(len*sizeof(double), "ilutpC:6"); 
	 memcpy(ilusch->L->ma[ii], w, len*sizeof(double));
	 for (j=0; j<len; j++)
	    ilusch->L->ja[ii][j] = iperm[jw[j]];
      }
/*---------------------------------------------------------------------
|     apply dropping strategy to U  (entries after diagonal)
|--------------------------------------------------------------------*/
#ifdef ILUTIME
msg_timer_stop(&t27);
msg_timer_start(&t28);
#endif
      len = 0;
      for (j=1; j<lenu; j++) {
	 if ( fabs(w[ii+j]) > drop6*tnorm ) {
	    len++;
	    w[ii+len] = w[ii+j];
	    jw[ii+len] = jw[ii+j];
	 }
      }
      lenu = len+1;
      len = lenu > fil6 ? fil6 : lenu;
      ilusch->U->nzcount[ii] = len;
      if (lenu > len+1)
         qsplitC(&w[ii+1], &jw[ii+1], lenu-1, len);
      ilusch->U->ma[ii] = (double *) Malloc(len*sizeof(double), "ilutpC:7" ); 
      ilusch->U->ja[ii] = (int *) Malloc(len*sizeof(int), "ilutpC:8" ); 
/*---------------------------------------------------------------------
|     determine next pivot
|--------------------------------------------------------------------*/
/* HERE  - all lines with dec included for counting pivots */
      imax = ii;
      xmax = fabs(w[imax]);
      xmax0 = xmax;
      /* icut = ii - 1 + mband - ii % mband;  */ 
      icut = ii-1 + mband; 
      dec = 0;
      for (k=ii+1; k<ii+len; k++) {
	 t = fabs(w[k]);
	 if ( (t > xmax) && (t*permtol > xmax0) && (jw[k] <= icut) ) {
	    imax = k;
	    xmax = t;
	    dec = 1;
	 }
      }
      if (dec == 1) decnt++;
/*---------------------------------------------------------------------
|     exchange w's
|--------------------------------------------------------------------*/
      tmp = w[ii];
      w[ii] = w[imax];
      w[imax] = tmp;
/*---------------------------------------------------------------------
|     update iperm and reverse iperm
|--------------------------------------------------------------------*/
      j = jw[imax];
      i = iperm[ii];
      iperm[ii] = iperm[j];
      iperm[j] = i;
      iprev[iperm[ii]] = ii;
      iprev[iperm[j]] = j;
/*---------------------------------------------------------------------
|     now store U in original coordinates
|--------------------------------------------------------------------*/
/*
  printf("   length of U = %d  ",len);
 */
      if (w[ii] == 0.0)  w[ii] = (0.0001+drop6)*tnorm;
      ilusch->U->ma[ii][0] = 1.0 / w[ii];
      ilusch->U->ja[ii][0] = ii;
      memcpy(&ilusch->U->ma[ii][1], &w[ii+1], (len-1)*sizeof(double));
      lenu = 0;
      for (k=ii+1; k<ii+len; k++)
	 ilusch->U->ja[ii][++lenu] = iperm[jw[k]];
#ifdef ILUTIME
msg_timer_stop(&t28);
#endif
   }
   cpermC(ilusch->U, iprev);
   cpermC(ilusch->L, iprev);
/* DO NOT PERMUTE ORIGINAL MATRIX
   cpermC(amat, iprev);
*/
#ifdef ILUTIME
   tnorm = 0;
   for (tt=22; tt<29; tt++) {
      printf(" loop %d      %e \n",tt,msg_timer(&tt));
      tnorm += msg_timer(&tt);
   }
   printf(" total        %e \n",tnorm);
#endif
/*---------------------------------------------------------------------
|     end main loop - now do clean up
|--------------------------------------------------------------------*/
   if (rmax > 0) {
      free(jw);
      free(w);
      free(jwrev);
      free(iprev);
   }
   printf("There were %d pivots\n",decnt);
/*---------------------------------------------------------------------
|     done  --  correct return
|--------------------------------------------------------------------*/
   return(0);
label994:
/*  Incomprehensible error. Matrix must be wrong.  */
   return(1);
/* label997:
   Memory allocation error.  
   return(2);  */
label998:
/*  illegal value for lfil or last entered  */
   return(5);
label9991:
/*  zero row encountered  */
   return(6);
}
/*---------------------------------------------------------------------
|     end of ilutpC
|--------------------------------------------------------------------*/

int ilutD(csptr amat, double *droptol, int *lfil, ilutptr ilusch)
{
/*---------------------------------------------------------------------- 
| ILUT -
| Converted to C so that dynamic memory allocation may be implememted.
| All indexing is in C format.
|----------------------------------------------------------------------
| ILUT factorization with dual truncation. 
|
| March 1, 2000 - dropping in U: keep entries > tau * diagonal entry
|----------------------------------------------------------------------
|
| on entry:
|========== 
| ( amat ) =  Matrix stored in SpaFmt struct.
| (ilusch) =  Pointer to ILUTfac struct
|
| lfil[5]  =  number nonzeros in L-part 
| lfil[6]  =  number nonzeros in U-part     ( lfil >= 0 )
|
| droptol[5] = threshold for dropping small terms in L during 
|              factorization. 
| droptol[6] = threshold for dropping small terms in U.
|
| On return:
|===========
|
| (ilusch) =  Contains L and U factors in an LUfact struct.
|             Individual matrices stored in SpaFmt structs.
|             On return matrices have C (0) indexing.
|
|       integer value returned:
|
|             0   --> successful return.
|             1   --> Error.  Input matrix may be wrong.  (The 
|                         elimination process has generated a
|                         row in L or U whose length is > n.)
|             2   --> Memory allocation error.
|             5   --> Illegal value for lfil or last.
|             6   --> zero row encountered.
|----------------------------------------------------------------------- 
| work arrays:
|=============
| jw, jwrev = integer work arrays of length n
| w         = real work array of length n. 
|----------------------------------------------------------------------- 
|     All processing is done using C indexing.
|--------------------------------------------------------------------*/
   int i, ii, j,jj, jcol, jpos, jrow, k, *jw=NULL, *jwrev=NULL;
   int len, lenu, lenl, rmax, fil5=lfil[5], fil6=lfil[6];
   double tnorm, t, s, fact, *w=NULL, drop5=droptol[5], drop6=droptol[6];
   int rowz, *rowj;
   double *rowm;

   ilusch->n = rmax = amat->n;
   if (rmax == 0) return(0); 
   if (rmax > 0) {
      jw = (int *) Malloc(rmax*sizeof(int), "ilutD:1" );
      w = (double *) Malloc(rmax*sizeof(double), "ilutD:2" );
      jwrev = (int *) Malloc(rmax*sizeof(int), "ilutD:3" );
   }
   if (fil5<0 || fil6<0 || rmax<=0) goto label9995;
/*---------------------------------------------------------------------
|    beginning of first main loop - L, U, L^{-1}F calculations
|--------------------------------------------------------------------*/
   for (j=0; j<rmax; j++)
      jwrev[j] = -1;
   for (ii=0; ii<rmax; ii++) {
      rowj = amat->ja[ii];
      rowm = amat->ma[ii];
      rowz = amat->nzcount[ii];
      for (k=0; k<rowz; k++)
	 if (rowm[k] != 0.0) goto label40;
      goto label9996;
/*---------------------------------------------------------------------
|     unpack amat in arrays w, jw, jwrev
|     WE ASSUME THERE IS A DIAGONAL ELEMENT
|--------------------------------------------------------------------*/
label40:
      lenu = 1;
      lenl = 0;
      w[ii] = 0.0;
      jw[ii] = ii;
      jwrev[ii] = ii;
      for (j=0; j<rowz; j++) {
	 jcol = rowj[j];
	 t = rowm[j];
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
	 rowm = ilusch->U->ma[jrow];
	 fact = w[jj] * rowm[0];
	 if ( fabs(fact) > drop5 ) {    /*  DROPPING IN L  */
	    rowj = ilusch->U->ja[jrow];
	    rowz = ilusch->U->nzcount[jrow];
/*---------------------------------------------------------------------
|     combine current row and row jrow
|--------------------------------------------------------------------*/
	    for (k=1; k<rowz; k++) {
	       s = fact * rowm[k];
	       j = rowj[k];
	       jpos = jwrev[j];
/*---------------------------------------------------------------------
|     dealing with U
|--------------------------------------------------------------------*/
	       if (j >= ii) {
/*---------------------------------------------------------------------
|     this is a fill-in element
|--------------------------------------------------------------------*/
		  if (jpos == -1) {
		     if (lenu > rmax) goto label9991;
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
		     if (lenl > rmax) goto label9991;
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
      for (j=0; j<lenl; j++)    /*  L block  */
	 jwrev[jw[j]] = -1;
      for (j=0; j<lenu; j++)    /*  U block  */
	 jwrev[jw[ii+j]] = -1;
/*---------------------------------------------------------------------
|     done reducing this row, now store L
|--------------------------------------------------------------------*/
      lenl = len > fil5 ? fil5 : len;
      ilusch->L->nzcount[ii] = lenl;
      if (len > lenl)
         qsplitC(w, jw, len, lenl);
      if (len > 0) {
	 ilusch->L->ja[ii] = (int *) Malloc(lenl*sizeof(int), "ilutD:4" ); 
	 ilusch->L->ma[ii] = (double *) Malloc(lenl*sizeof(double), "ilutD:5"); 
	 memcpy(ilusch->L->ja[ii], jw, lenl*sizeof(int));
	 memcpy(ilusch->L->ma[ii], w, lenl*sizeof(double));
      }
/*---------------------------------------------------------------------
|     store the diagonal element of U
|
|     dropping in U if size is less than drop1 * diagonal entry
|--------------------------------------------------------------------*/
      t = w[ii];
      tnorm = fabs(t);
      len = 0;
      for (j=1; j<lenu; j++) {
	 if ( fabs(w[ii+j]) > drop6*tnorm ) {
	    w[len] = w[ii+j];
	    jw[len] = jw[ii+j];
	    len++;
	 }
      }
      lenu = len+1 > fil6 ? fil6 : len+1;
      ilusch->U->nzcount[ii] = lenu;
      jpos = lenu-1;
      if (len > jpos)
         qsplitC(w, jw, len, jpos);
      ilusch->U->ma[ii] = (double *) Malloc(lenu*sizeof(double), "ilutD:6" ); 
      ilusch->U->ja[ii] = (int *) Malloc(lenu*sizeof(int), "ilutD:7" ); 
      if (t == 0.0) t=(0.0001+drop6)*tnorm;
      ilusch->U->ma[ii][0] = 1.0 / t;
      ilusch->U->ja[ii][0] = ii;
/*---------------------------------------------------------------------
|     copy the rest of U
|--------------------------------------------------------------------*/
      memcpy(&ilusch->U->ja[ii][1], jw, jpos*sizeof(int));
      memcpy(&ilusch->U->ma[ii][1], w, jpos*sizeof(double));
   }
/*---------------------------------------------------------------------
|     end main loop - now do clean up
|--------------------------------------------------------------------*/
   if (rmax > 0) {
      free(jw);
      free(w);
      free(jwrev);
   }
/*---------------------------------------------------------------------
|     done  --  correct return
|--------------------------------------------------------------------*/
   return(0);
label9991:
/*  Incomprehensible error. Matrix must be wrong.  */
   return(1);
/* label9992:
   Memory allocation error.  
   return(2);  */
label9995:
/*  illegal value for lfil or last entered  */
   return(5);
label9996:
/*  zero row encountered  */
   return(6);
}
/*---------------------------------------------------------------------
|     end of ilutNEW
|--------------------------------------------------------------------*/
