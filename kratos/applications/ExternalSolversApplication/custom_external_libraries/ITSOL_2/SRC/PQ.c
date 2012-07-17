#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "globheads.h"
#include "protos.h"

#define  ALPHA  0.00001

/*--------------------end protos */
int PQperm(csptr mat, int bsize, int *Pord, int *Qord, int *nnod, 
	     double tol) {
/*--------------------------------------------------------------------- 
| algorithm for nonsymmetric  block selection - 
|----------------------------------------------------------------------
|     Input parameters:
|     -----------------
|     (mat)  =  matrix in SpaFmt format
|     
|     tol    =  a tolerance for excluding a row from B block 
|
|     bsize not used here - it is used in arms2.. 
|
|     Output parameters:
|     ------------------ 
|     Pord   = row permutation array.  Row number i will become row number 
|		Pord[i] in permuted matrix. (Old to new labels) 
|     Qord   = column permutation array.  Column number j will become 
|		row number Qord[j] in permuted matrix.
|             [destination lists] - i.e., old to new arrays= 
|     
|     nnod   = number of elements in the B-block 
|     
|---------------------------------------------------------------------*/ 
/*--------------------   local variables   */
  int *icor, *jcor, *row; 
  int i, j, ii, k, col, jj, rnz, nzi, n=mat->n, count, numnode;
  double  aij, rn, *mrow;  
/*-----------------------------------------------------------------------*/  
  for (j=0; j<n; j++) {
    Pord[j] = -1;
    Qord[j] = -1; 
  }
  icor = (int *) malloc(n*sizeof(int));
  jcor = (int *) malloc(n*sizeof(int));
  if ( (icor==NULL) || (jcor==NULL) ) return 1;
  numnode = 0;
  count = 0;
/*-------------------- wDiag selects candidate entries in a sorted oder */
  i = 1; 
  preSel(mat, icor, jcor, i, tol, &count) ;
/*-------------------- add entries one by one to diagnl */
  /* needs recoding so as to scan rows only once instead of 2 */
  for (i = 0; i<count; i++){ 
    ii = icor[i];
    jj = jcor[i]; 
    if (Qord[jj] != -1) continue;
/*-------------------- */
    row = mat->ja[ii]; 
    mrow = mat->ma[ii];
    nzi = mat->nzcount[ii] ;
/*-------------------- rnz = already assigned cols (either in L or F) */
    rn = fabs(mrow[0]); 
    rnz = (nzi-1) ; 
    for (k=0; k < nzi; k++) {
      aij = fabs(mrow[k]);
      col = row[k];
      if (Qord[col] >=0 ) {
	rn -= aij; 
	rnz-- ; 
      }
      else if (Qord[col] == -2) {
	rnz--;
      }
    } 
    if (rn < 0.0) continue;   
    Pord[ii] = numnode;
    Qord[jj] = numnode;
    numnode++; 
/*-------------------- acceptance test among others */    
    for (k=0; k < nzi; k++) {
      col = row[k];
      if (Qord[col] != -1) continue;
      aij = fabs(mrow[k]);
      if (rnz*aij > rn) 
	Qord[col] = -2;
      else 
	rn -= aij;
      rnz--;
    }
  }
  /*-------------------- number of B nodes */
  *nnod = numnode; 
  /* printf(" nnod found = %d \n",*nnod);  */ 
/*--------------------------------------------------
|    end-main-loop - complete permutation arrays
|-------------------------------------------------*/
  for (i=0; i<n; i++)
    if (Pord[i] < 0) 
      Pord[i] = numnode++;
  
  if (numnode != n) {
    printf("  ** counting error - type 1 \n"); return 1; }   
  numnode = *nnod;
  for (j=0; j<n; j++)
    if (Qord[j] < 0)
      Qord[j] = numnode++;
/*--------------------              */ 
  if (numnode != n) {
    printf(" **  counting error type 2 \n"); return 1; }

/*-------------------- debugging - check permutations */
  /* 
     printf(" checking P  and Q  :    -->  \n") ;
     check_perm(n, Pord) ;
     check_perm(n, Qord) ;
  */
/*--------------------------------------------------
|  clean up before returning
|-------------------------------------------------*/
  free(icor);
  free(jcor);
  return 0;
}
/*---------------------------------------------------------------------
|-----end-of-indsetPQ--------------------------------------------------
|--------------------------------------------------------------------*/

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

/*-----------------------------------------------------------------------*/
   riord = (int *) Malloc(n*sizeof(int), "indsetC:1" );
   w     = (double *) Malloc(n*sizeof(double), "indsetC:2" );
   matT  = (csptr) Malloc(sizeof(SparMat), "indsetC:3" );
/*  	 call weights to compute the weights for  input matrix.. */
   setupCS(matT, mat->n,1);
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

int preSel(csptr mat, int *icor, int *jcor, int job, double tol, int *count) 
{
/*---------------------------------------------------------------------
| does a preselection of possible diagonal entries. will return a list
| of "count" bi-indices representing "good" entries to be selected as 
| diagonal elements -- the order is important (from best to
| to worst). The list is in the form (icor(ii), jcor(ii)) 
|
|      ON ENTRY: 
|       mat   = matrix in csptr format 
|       tol   = tolerance used for selecting best rows -|
|       job   = indicates whether or not to permute the max entry in 
|               each row to first position 
|        NOTE: CAN RECODE BY HAVING JCOR CARRY THE INDEX IN ROW[I] WHERE
|              MAX IS LOCATED.. 
|       
|      ON RETURN: 
|       icor  = list of row indices of entries selected 
|       jcor  = list of column indices of entries selected 
|       count = number of entries selected (size of B block) 
|--------------------------------------------------------------------*/
  int i, k, kmax, n=mat->n, col, jmax, countL;
  int *nz, *jcol; 
  double *mrow, rownorm, *weight, t, tmax, wmax;

/*--------------------begin */
/*-----------------------------------------------------------------------*/
  nz =mat->nzcount;
  weight = (double *) malloc(n*sizeof(double));
  if ( weight==NULL) return 1;  
  /*-------------------- compute max entry for each row */
  wmax = 0.0; 
  for (i=0; i<n; i++) {
    jcol = mat->ja[i];
    mrow = mat->ma[i];
    tmax = 0.0; 
    kmax = 0; 
    rownorm = 0.0; 
    for (k = 0; k<nz[i]; k++) {
      col = jcol[k] ; 
      t = fabs(mrow[k]); 
      if (t != 0.0) {
	rownorm += t; 
	if (tmax < t) {
	  tmax = t;
	  kmax = k;
	}
      }
    }
    jmax = jcol[kmax];
    jcor[i] = jmax; 
    if (job && kmax != 0) {
      t = mrow[kmax];
      mrow[kmax] = mrow[0];
      mrow[0] = t;
      jcol[kmax] = jcol[0];
      jcol[0] = jmax;
    }
/*-------------------- save max diag. dominance ratio  */
    t = tmax / rownorm;  
    if (wmax < t)  wmax = t; 
    weight[i] = t;
    /* remove!! ALREADY ASSIGNED  */
    jcor[i] = jmax;
  }
/*-------------------- now select according to tol */
  countL = 0; 
  for (i=0; i<n; i++) {
    t = weight[i] ;
    col = jcor[i];
    if (t < wmax*tol) continue ;
    weight[countL] =  t /((double) nz[i]) ; 
    icor[countL] = i; 
    jcor[countL] = col;
    countL++;
  }
/*-------------------- sort them  */
  qsortR2I(weight, icor, jcor, 0, countL-1);
  *count = countL;
  free(weight);
  return 0;
}
/*---------------------------------------------------------------------
|---- end of preSel ---------------------------------------------------
|--------------------------------------------------------------------*/

