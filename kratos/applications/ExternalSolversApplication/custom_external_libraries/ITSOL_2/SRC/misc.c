#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "globheads.h"
#include "protos.h"

int dumpCooMat(csptr A, int nglob, int,  FILE *ft);

int qsplitC(double *a, int *ind, int n, int Ncut)
{
/*----------------------------------------------------------------------
|     does a quick-sort split of a real array.
|     on input a[0 : (n-1)] is a real array
|     on output is permuted such that its elements satisfy:
|
|     abs(a[i]) >= abs(a[Ncut-1]) for i < Ncut-1 and
|     abs(a[i]) <= abs(a[Ncut-1]) for i > Ncut-1
|
|     ind[0 : (n-1)] is an integer array permuted in the same way as a.
|---------------------------------------------------------------------*/
   double tmp, abskey;
   int j, itmp, first, mid, last, ncut;
   ncut = Ncut - 1;

   first = 0;
   last = n-1;
   if (ncut<first || ncut>last) return 0;
/* outer loop -- while mid != ncut */
do{
   mid = first;
   abskey = fabs(a[mid]);
  for (j=first+1; j<=last; j++) {
     if (fabs(a[j]) > abskey) {
        mid = mid+1;
	 tmp = a[mid];
	 itmp = ind[mid];
	 a[mid] = a[j];
	 ind[mid] = ind[j];
	 a[j]  = tmp;
	 ind[j] = itmp;
      }
   }
/*-------------------- interchange */
   tmp = a[mid];
   a[mid] = a[first];
   a[first]  = tmp;
   itmp = ind[mid];
   ind[mid] = ind[first];
   ind[first] = itmp;
/*-------------------- test for while loop */
   if (mid == ncut) break;
   if (mid > ncut) 
      last = mid-1;
   else
      first = mid+1;
   }while(mid != ncut);
   
   return 0;
}
/*--------------- end of qsplitC ----------------------------------------
|---------------------------------------------------------------------*/
int SparTran(csptr amat, csptr bmat, int job, int flag)
{
/*----------------------------------------------------------------------
| Finds the transpose of a matrix stored in SpaFmt format.
|
|-----------------------------------------------------------------------
| on entry:
|----------
| (amat) = a matrix stored in SpaFmt format.
|
| job    = integer to indicate whether to fill the values (job.eq.1)
|          of the matrix (bmat) or only the pattern.
|
| flag   = integer to indicate whether the matrix has been filled
|          0 - no filled
|          1 - filled
|
| on return:
| ----------
| (bmat) = the transpose of (mata) stored in SpaFmt format.
|
| integer value returned:
|             0   --> successful return.
|             1   --> memory allocation error.
|---------------------------------------------------------------------*/
  int i, j, *ind, pos, size=amat->n, *aja;
  double *ama=NULL;
  ind = (int *) Malloc(size*sizeof(int), "SparTran:1" );
  for (i=0; i<size; i++)
    ind[i] = 0;
  if(!flag) {
/*--------------------  compute lengths  */
    for (i=0; i<size; i++) {
      aja = amat->ja[i];
      for (j=0; j<amat->nzcount[i]; j++)
	ind[aja[j]]++;
    }
/*--------------------  allocate space  */
    for (i=0; i<size; i++) {
      bmat->ja[i] = (int *) Malloc(ind[i]*sizeof(int), "SparTran:2" );
      bmat->nzcount[i] = ind[i];
      if (job == 1) {
	bmat->ma[i] = (double *) Malloc(ind[i]*sizeof(double), "SparTran:3" );
      }
      ind[i] = 0;
    }
  }
/*--------------------  now do the actual copying  */
  for (i=0; i<size; i++) {
    aja = amat->ja[i];
    if (job == 1)
      ama = amat->ma[i];
    for (j=0; j<amat->nzcount[i]; j++) {
      pos = aja[j];
      bmat->ja[pos][ind[pos]] = i;
      if (job == 1)
	bmat->ma[pos][ind[pos]] = ama[j];
      ind[pos]++;
    }
  }
  free(ind);
  return 0;
}
/*-------------- end of SparTran ---------------------------------------
|---------------------------------------------------------------------*/

void swapj(int v[], int i, int j){
  int temp;
  temp = v[i];
  v[i] = v[j];
  v[j] = temp;
}
void swapm(double v[], int i, int j) {
  double temp;
  temp = v[i];
  v[i] = v[j];
  v[j] = temp;
}


int roscalC(csptr mata, double *diag, int nrm)
{
/*---------------------------------------------------------------------
|
| This routine scales each row of mata so that the norm is 1.
|
|----------------------------------------------------------------------
| on entry:
| mata  = the matrix (in SpaFmt form)
| nrm   = type of norm
|          0 (\infty),  1 or 2
|
| on return
| diag  = diag[j] = 1/norm(row[j])
|
|     0 --> normal return
|     j --> row j is a zero row
|--------------------------------------------------------------------*/
/*   local variables    */
   int i, k;
   double *kr, scal;
   for (i=0; i<mata->n; i++) {
      scal = 0.0;
      kr = mata->ma[i];
      if (nrm == 0) {
	 for (k=0; k<mata->nzcount[i]; k++)
	    if (fabs(kr[k]) > scal) scal = fabs(kr[k]);
      }
      else if (nrm == 1) {
         for (k=0; k<mata->nzcount[i]; k++)
            scal += fabs(kr[k]);
      }
      else {  /* nrm = 2 */
         for (k=0; k<mata->nzcount[i]; k++)
            scal += kr[k]*kr[k];
      }
      if (nrm == 2) scal = sqrt(scal);
      if (scal == 0.0) {
	scal = 1.0; 
	/* YS. return i+1; */
      }
      else 
	scal = 1.0 / scal;
      diag[i] = scal;
      for (k=0; k<mata->nzcount[i]; k++)
	 kr[k] = kr[k] * scal;
   }
   return 0;
}
/*---------------end of roscalC-----------------------------------------
----------------------------------------------------------------------*/
int coscalC(csptr mata, double *diag, int nrm)
{
/*---------------------------------------------------------------------
|
| This routine scales each column of mata so that the norm is 1.
|
|----------------------------------------------------------------------
| on entry:
| mata  = the matrix (in SpaFmt form)
| nrm   = type of norm
|          0 (\infty),  1 or 2
|
| on return
| diag  = diag[j] = 1/norm(row[j])
|
|     0 --> normal return
|     j --> column j is a zero column
|--------------------------------------------------------------------*/
/*   local variables    */
   int i, j, k;
   double *kr;
   int *ki;
   for (i=0; i<mata->n; i++)
      diag[i] = 0.0;
/*---------------------------------------
|   compute the norm of each column
|--------------------------------------*/
   for (i=0; i<mata->n; i++) {
      kr = mata->ma[i];
      ki = mata->ja[i];
      if (nrm == 0) {
	 for (k=0; k<mata->nzcount[i]; k++) {
	    j = ki[k];
	    if (fabs(kr[k]) > diag[j]) diag[j] = fabs(kr[k]);
	 }
      }
      else if (nrm == 1) {
         for (k=0; k<mata->nzcount[i]; k++)
            diag[ki[k]] += fabs(kr[k]);
      }
      else {  /*  nrm = 2 */
         for (k=0; k<mata->nzcount[i]; k++)
            diag[ki[k]] += kr[k]*kr[k];
      }
   }
   if (nrm == 2) {
      for (i=0; i<mata->n; i++)
	 diag[i] = sqrt(diag[i]);
   }
/*---------------------------------------
|   invert
|--------------------------------------*/
   for (i=0; i<mata->n; i++) {
      if (diag[i] == 0.0)
	/* return i+1;*/
	diag[i] = 1.0; 
      else 
	 diag[i] = 1.0 / diag[i];
   }
/*---------------------------------------
|   C = A * D
|--------------------------------------*/
   for (i=0; i<mata->n; i++) {
      kr = mata->ma[i];
      ki = mata->ja[i];
      for (k=0; k<mata->nzcount[i]; k++)
	 kr[k] = kr[k] * diag[ki[k]];
   }
   return 0;
}
/*---------------end of coscalC-----------------------------------------
----------------------------------------------------------------------*/
void dscale(int n, double *dd, double *x, double * y)
{ 
/* Computes  y == DD * x                               */
/* scales the vector x by the diagonal dd - output in y */
  int k;
  for (k=0; k<n; k++) 
    y[k] = dd[k]*x[k];
}

  
void qsortC(int *ja, double *ma, int left, int right, int abval)
{
/*----------------------------------------------------------------------
|
| qqsort: sort ma[left]...ma[right] into decreasing order
| from Kernighan & Ritchie
|
| ja holds the column indices
| abval = 1: consider absolute values
|         0: values
|
|---------------------------------------------------------------------*/
  int i, last;

  if (left >= right)  return;
  if (abval) {
    swapj(ja, left, (left+right)/2);
    swapm(ma, left, (left+right)/2);
    last = left;
    for (i=left+1; i<=right; i++) {
      if (fabs(ma[i]) > fabs(ma[left])) {
	swapj(ja, ++last, i);
	swapm(ma, last, i);
      }
    }
    swapj(ja, left, last);
    swapm(ma, left, last);
    qsortC(ja, ma, left, last-1, abval);
    qsortC(ja, ma, last+1, right, abval);
  }
  else {
    swapj(ja, left, (left+right)/2);
    swapm(ma, left, (left+right)/2);
    last = left;
    for (i=left+1; i<=right; i++) {
      if (ma[i] > ma[left]) {
	swapj(ja, ++last, i);
	swapm(ma, last, i);
      }
    }
    swapj(ja, left, last);
    swapm(ma, left, last);
    qsortC(ja, ma, left, last-1, abval);
    qsortC(ja, ma, last+1, right, abval);
  }
}

void printmat(FILE *ft, csptr A, int i0, int i1){
/*-------------------------------------------------------------+
| to dump rows i0 to i1 of matrix for debugging purposes       |
|--------------------------------------------------------------*/
  int i, k, nzi;
  int *row;
  double *rowm;
  for (i=i0; i<i1; i++)    {
    nzi = A->nzcount[i];
    row = A->ja[i];
    rowm = A->ma[i];
    for (k=0; k< nzi; k++){
      fprintf(ft," row %d  a  %e ja %d \n", i+1, rowm[k], row[k]+1);
    }
  }
}

void qsortR2I(double *wa, int *cor1, int *cor2, int left, int right){
/*----------------------------------------------------------------------
|
| qqsort: sort wa[left]...wa[right] into decreasing order
| from Kernighan & Ritchie
|
|---------------------------------------------------------------------*/
  int i, last;
  
  if (left >= right)  return;
  
  swapm(wa, left, (left+right)/2);
  swapj(cor1, left, (left+right)/2);
  swapj(cor2, left, (left+right)/2);
  last = left;
  for (i=left+1; i<=right; i++) {
    if (wa[i] > wa[left]) {
      swapm(wa, ++last, i);
      swapj(cor1, last, i);
      swapj(cor2, last, i);
    }
  }
  swapm(wa, left, last);
  swapj(cor1, left, last);
  swapj(cor2, left, last);
  qsortR2I(wa, cor1, cor2, left, last-1);
  qsortR2I(wa, cor1, cor2, last+1, right);
}


void qsort2C(int *ja, double *ma, int left, int right, int abval){
/*----------------------------------------------------------------------
|
| qqsort: sort ma[left]...ma[right] into increasing order
| from Kernighan & Ritchie
|
| ja holds the column indices
| abval = 1: consider absolute values
|         0: values
|
|---------------------------------------------------------------------*/
  int i, last;
  if (left >= right)  return;
  if (abval) {
    swapj(ja, left, (left+right)/2);
    swapm(ma, left, (left+right)/2);
    last = left;
    for (i=left+1; i<=right; i++) {
      if (fabs(ma[i]) < fabs(ma[left])) {
	swapj(ja, ++last, i);
	swapm(ma, last, i);
      }
    }
    swapj(ja, left, last);
    swapm(ma, left, last);
    qsort2C(ja, ma, left, last-1, abval);
    qsort2C(ja, ma, last+1, right, abval);
  }
  
  else {
    swapj(ja, left, (left+right)/2);
    swapm(ma, left, (left+right)/2);
    last = left;
    for (i=left+1; i<=right; i++) {
      if (ma[i] < ma[left]) {
	swapj(ja, ++last, i);
	swapm(ma, last, i);
      }
    }
    swapj(ja, left, last);
    swapm(ma, left, last);
    qsort2C(ja, ma, left, last-1, abval);
    qsort2C(ja, ma, last+1, right, abval);
  }
}

void qqsort(int *ja, double *ma, int left, int right){
/*----------------------------------------------------------------------
|
| qqsort: sort ja[left]...ja[right] into increasing order
| from Kernighan & Ritchie
|
| ma holds the real values
|
|---------------------------------------------------------------------*/
  int i, last;
  if (left >= right)  return;
  swapj(ja, left, (left+right)/2);
  swapm(ma, left, (left+right)/2);
  last = left;
  for (i=left+1; i<=right; i++) {
    if (ja[i] < ja[left]) {
      swapj(ja, ++last, i);
      swapm(ma, last, i);
    }
  }
  swapj(ja, left, last);
  swapm(ma, left, last);
  qqsort(ja, ma, left, last-1);
  qqsort(ja, ma, last+1, right);
}

void hilosort(csptr mat, int abval, int hilo){
/*----------------------------------------------------------------------
|
| This routine sorts the entries in each row of a matrix from hi to low.
|
|-----------------------------------------------------------------------
| on entry:
|----------
| (mat) = a matrix stored in SpaFmt format.
|
| abval =   1: use absolute values of entries
|           0: use values
|
| hilo  =   1: sort in decreasing order
|           0: sort in increasing order
|
|
| on return:
| ----------
| (mat) = (mat) where each row is sorted.
|
|---------------------------------------------------------------------*/
  int j, n=mat->n, *nnz=mat->nzcount;
  
  if (hilo)
    for (j=0; j<n; j++)
      qsortC(mat->ja[j], mat->ma[j], 0, nnz[j]-1, abval);
  
  else
    for (j=0; j<n; j++)
      qsort2C(mat->ja[j], mat->ma[j], 0, nnz[j]-1, abval);
  
  return;
}
/*------- end of hilosort ----------------------------------------------
|---------------------------------------------------------------------*/

void qsort3i(int *wa, int *cor1, int *cor2, int left, int right)
/*----------------------------------------------------------------------
|
| qqsort: sort wa[left]...wa[right] into increasing order
| from Kernighan & Ritchie
|
|---------------------------------------------------------------------*/
{
   int i, last;

   if (left >= right)  return;

   swapj(wa, left, (left+right)/2);
   swapj(cor1, left, (left+right)/2);
   swapj(cor2, left, (left+right)/2);
   last = left;
   for (i=left+1; i<=right; i++) {
      if (wa[i] < wa[left]) {
	 swapj(wa, ++last, i);
	 swapj(cor1, last, i);
	 swapj(cor2, last, i);
      }
   }
   swapj(wa, left, last);
   swapj(cor1, left, last);
   swapj(cor2, left, last);
   qsort3i(wa, cor1, cor2, left, last-1);
   qsort3i(wa, cor1, cor2, last+1, right);
}

int dumpArmsMat(arms PreSt, FILE *ft){
  int lev, nnz, nglob=0, old=0;
  p4ptr levmat = PreSt->levmat; 
  ilutptr ilus = PreSt->ilus; 
  int n = levmat->n;
  int nlev=PreSt->nlev; 
  FILE *dummy=NULL;

  nnz = nnz_arms(PreSt,dummy)-nnz_cs(ilus->C);

  fprintf(ft," %d %d %d \n", n, n,nnz);

  old = 0;
  for (lev = 0; lev<nlev; lev++) {
/*-------------------- Dump L U of B */      
    nglob+= levmat->nB;
    dumpCooMat(levmat->L, old, old, ft);
    dumpCooMat(levmat->U, old, old, ft);
    dumpCooMat(levmat->E, nglob, old, ft);
    dumpCooMat(levmat->F, old, nglob, ft);
    levmat = levmat->next;
    if (levmat == NULL) break;
    old = nglob;
  }
  dumpCooMat(ilus->L, nglob, nglob, ft);
  dumpCooMat(ilus->U, nglob, nglob, ft);
  return(0);
}

int dumpCooMat(csptr A, int shiftR, int shiftC, FILE *ft){
  /*-------------------- dump matrix in coo format
    -------------------- matlab style */
  int n, i, k, nzi;
  int *row;
  n = A->n;
  double *rowm;
  for (i=0; i<n; i++) {
    nzi = A->nzcount[i];
    row = A->ja[i];
    rowm = A->ma[i];
    for (k=0; k< nzi; k++){
/*   colmnidx  rowidx  val */
      fprintf(ft," %d  %d  %e \n", row[k]+shiftC, i+shiftR, rowm[k]);
    }
  }
  return(0);
}

int outputLU( iluptr lu, char *filename ){
/*----------------------------------------------------------------------
| Output the pattern of L\U, which can be loaded by matlab
----------------------------------------------------------------------*/
  FILE *fmatlab = fopen( filename, "w" );
  int n = lu->n, i, j, nzcount;
  csptr L = lu->L, U = lu->U;
  
  if( !fmatlab ) return -1;
  fprintf( fmatlab, "%d %d 0\n", n, n );
  for( i = 0; i < n; i++ ) {
    nzcount = L->nzcount[i];
    for( j = 0; j < nzcount; j++ )
      fprintf( fmatlab, "%d %d 1\n", i+1, L->ja[i][j]+1 );
  }
  for( i = 0; i < n; i++ ) {
    nzcount = U->nzcount[i];
    for( j = 0; j < nzcount; j++ )
      fprintf( fmatlab, "%d %d 1\n", i+1, U->ja[i][j]+1 );
  }
  for( i = 0; i < n; i++ )
    fprintf( fmatlab, "%d %d 1\n", i+1, i+1 );
  fclose( fmatlab );
  return 0;
}

int checkperm(int *p, int n) {
  /*-------------------- checks the validity of a permutation [for 
    debugging purposes.] 
    Return codes:
    0  -- permutation is valid
    1  -- a value perm[?] is outside the range 0--(n-1)
    2  -- perm[i] hit a value between 0--(n-1) more than once.
  */
  int *work;
  int k,i ;
  work = Malloc(n*sizeof(int)," check perm:work ");
  for (k=0; k<n ;k++)
    work[k] = -1;
  for (k=0; k<n ;k++){
    i = p[k];
    if ( (i <0) | (i>=n) ) 
      return(1);
    if  (work[i]>=0) 
      return (2);
    work[i] = k;
  }
  free(work);
  return(0);
}

void qsortR1I(double *wa, int *cor1, int left, int right)
/*----------------------------------------------------------------------
|
| qqsort: sort wa[left]...wa[right] into decreasing order
| from Kernighan & Ritchie
|
|---------------------------------------------------------------------*/
{
   int i, last;
   void swapj(int *, int, int);
   void swapm(double *, int, int);

   if (left >= right)  return;

   swapm(wa, left, (left+right)/2);
   swapj(cor1, left, (left+right)/2);
   last = left;
   for (i=left+1; i<=right; i++) {
      if (wa[i] > wa[left]) {
	 swapm(wa, ++last, i);
	 swapj(cor1, last, i);
      }
   }
   swapm(wa, left, last);
   swapj(cor1, left, last);
   qsortR1I(wa, cor1, left, last-1);
   qsortR1I(wa, cor1, last+1, right);
}



  
  
