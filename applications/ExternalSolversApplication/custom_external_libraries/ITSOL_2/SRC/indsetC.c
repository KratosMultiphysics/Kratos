#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "globheads.h"
#include "protos.h"

void *Malloc( int, char * );

int add2is(int *last, int nod, int *iord, int *riord)
{
/*----------------------------------------------------------------------
|   adds element nod to independent set
|---------------------------------------------------------------------*/
   (*last)++;
   iord[nod] = *last;
   riord[*last] = nod;
   return 0;
}
/*---------------------------------------------------------------------
|---- end of add2is ---------------------------------------------------
|--------------------------------------------------------------------*/
int add2com(int *nback, int nod, int *iord, int *riord) 
{
/*----------------------------------------------------------------------
|   adds element nod to independent set
|---------------------------------------------------------------------*/
  iord[nod] = *nback;
  riord[*nback] = nod;
  (*nback)--;
  return 0;
}
/*---------------------------------------------------------------------
|---- end of add2com --------------------------------------------------
|--------------------------------------------------------------------*/
int indsetC(csptr mat, int bsize, int *iord, int *nnod, double tol) 
{
/*--------------------------------------------------------------------- 
| greedy algorithm for independent set ordering -- 
|----------------------------------------------------------------------
|     Input parameters:
|     -----------------
|     (mat)  =  matrix in SpaFmt format
|     
|     bsize  =  integer (input) the target size of each block.
|               each block is of size >= bsize. 
|
|     w      =  weight factors for the selection of the elements in the
|               independent set. If w(i) is small i will be left for the
|               vertex cover set. 
|
|     tol    =  a tolerance for excluding a row from independent set.
|
|     Output parameters:
|     ------------------ 
|     iord   = permutation array corresponding to the independent set 
|     ordering.  Row number i will become row number iord[i] in 
|     permuted matrix.
|     
|     nnod   = (output) number of elements in the independent set. 
|     
|----------------------------------------------------------------------- 
|     the algorithm searches nodes in lexicographic order and groups
|     the (BSIZE-1) nearest nodes of the current to form a block of
|     size BSIZE. The current algorithm does not use values of the matrix.
|---------------------------------------------------------------------*/ 
/*   local variables   */
   int nod, jcount, lastlev, begin, last0, last, nback, mid,
     j1, j2, jcol, inod, jnod, j, k, jcount0, begin0, *rowj;
   int prog, n=mat->n, *riord;
   double *w;
   csptr matT,gmat;  
   int weightsC(csptr, double *); 
/*-----------------------------------------------------------------------*/
   riord = (int *) Malloc(n*sizeof(int), "indsetC:1" );
   w     = (double *) Malloc(n*sizeof(double), "indsetC:2" );
   matT  = (csptr) Malloc(sizeof(SparMat), "indsetC:3" );
/*  	 call weights to compute the weights for  input matrix.. */
   setupCS(matT, mat->n, 1);
   SparTran(mat, matT, 1, 0);
   SparTran(matT, mat, 1, 1); 
   weightsC(mat, w); 
/*---------------------------------------------------------------------- 
| scan all nodes first to eliminate those not satisfying DD criterion 
+----------------------------------------------------------------------*/
   nback = n-1; 
   nod = 0;
   for(j=0; j<n; j++)
     iord[j] = -1; 
   for(j=0; j<n; j++) {
     if (w[j] < tol) {
       add2com(&nback, j, iord, riord);
       nod++;
     }
   }
   last = -1;
   for (nod=0; nod<n; nod++) {	
     while (iord[nod] != -1)   {
       if (++nod >= mat->n) goto label50;
     }
/*-------------------- initialize level-set - contains nod (only)*/
     add2is(&last, nod, iord, riord);
     begin   = last;
     begin0  = begin; 
     lastlev = begin;
     jcount  = 1;
/*----------------------------------------------------------------------
|     put all the nearest neighbor nodes of the current node into
|     the block until the number is BSIZE.
|---------------------------------------------------------------------*/
     prog = 1;
     while (jcount < bsize && prog) {
/*--------------------   traverse all the current level-set   */
       last0 = last;
       jcount0 = jcount;
       for (inod=begin; inod<=last0; inod++) {
	 jnod = riord[inod]; 
/*--------------------   This assumes A is not symmetric.   */
	 gmat = mat; 
	 for (k=0; k<2; k++) {
	   rowj = gmat->ja[jnod];
	   for (j=0; j<gmat->nzcount[jnod]; j++) {
	     jcol = rowj[j];
	     if (iord[jcol] == -1 ) {	
	       add2is(&last, jcol, iord, riord);
	       jcount++;
	     }
	   }
	   gmat = matT; 	
	 }
       }
       prog = jcount > jcount0 ? 1 : 0;
       lastlev = begin;
       begin = last0+1;
     }
/*-----------------------------------------------------------------------
| the neighbors of elements of last level go to the complement   
| gmat loop over original matrix and its transpose 
+-----------------------------------------------------------------------*/ 
     gmat = mat; 
     for (k=0; k<2; k++) {
       for (inod=lastlev; inod<=last; inod++)  {	
	 jnod = riord[inod]; 
	 rowj = gmat->ja[jnod];
	 for (j=0; j<gmat->nzcount[jnod]; j++){	
	   jcol = rowj[j];
	   if (iord[jcol] == -1) 
	     add2com(&nback, jcol, iord, riord);
	 }
       }
       gmat = matT; 	
     }
/*   reverse ordering for this level   */
     mid = (begin0+last) / 2;
     for (inod=begin0; inod<=mid; inod++) {
       j = last - inod + begin0;
       jnod = riord[inod];
       riord[inod] = riord[j];
       riord[j] = jnod;
     }
   }
/*--------------------------------------------------
|  end-main-loop
|-------------------------------------------------*/
/*-------------------- relabel nodes of vertex cover   */
label50:
   *nnod = last;
   j1 = *nnod;
   for (j2=*nnod+1; j2<n; j2++) { 
     if (iord[riord[j2]] > -1) {
       if (++j1 != j2) {
	 j = riord[j2];
	 riord[j2] = riord[j1];
	 riord[j1] = j;
       }
     }
   }
/*-------------------- obtain reverse permutation array   */
   for (j=0; j<n; j++)
     iord[riord[j]] = j;
   (*nnod)++;
   cleanCS(matT); 
   free(riord);
   free(w);
   return 0;
}
/*---------------------------------------------------------------------
|-----end-of-indsetC---------------------------------------------------
|--------------------------------------------------------------------*/
int weightsC(csptr mat, double *w)       
{
/*---------------------------------------------------------------------
|     defines weights based on diagonal dominance ratios
|--------------------------------------------------------------------*/
   int irow, k, n=mat->n, *kj, kz;
   double tdia, wmax=0.0, tnorm, *kr;
   for (irow=0; irow<n; irow++) {
      kz = mat->nzcount[irow];
      kr = mat->ma[irow];
      kj = mat->ja[irow];
      tnorm = 0.0;
      tdia = 0.0;
      for (k=0; k<kz; k++) {
	 if (kj[k] == irow) tdia = fabs(kr[k]); 
	 tnorm += fabs(kr[k]);
      }
      if (tnorm > 0.0)
	 tnorm =  tdia / tnorm;
      w[irow] = tnorm;
      if (tnorm > wmax) wmax = tnorm;
   }
   for (irow=0; irow<n; irow++)
      w[irow] = w[irow]/wmax;
   return 0;
}
/*---------------------------------------------------------------------
|---- end of weightsC -------------------------------------------------
|--------------------------------------------------------------------*/
