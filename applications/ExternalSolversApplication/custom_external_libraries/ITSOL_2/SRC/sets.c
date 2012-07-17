#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include "globheads.h"
#include "protos.h"

void errexit( char *f_str, ... ){
  va_list argp;
  char out1[256], out2[256];

  va_start(argp, f_str);
  vsprintf(out1, f_str, argp);
  va_end(argp);

  sprintf(out2, "Error! %s\n", out1);

  fprintf(stdout, out2);
  fflush(stdout);

  exit( -1 );
}

void *Malloc( int nbytes, char *msg )
{
  void *ptr;

  if (nbytes == 0)
    return NULL;

  ptr = (void *)malloc(nbytes);
  if (ptr == NULL)
    errexit( "Not enough mem for %s. Requested size: %d bytes", msg, nbytes );

  return ptr;
}

int setupCS(csptr amat, int len, int job)
{
/*----------------------------------------------------------------------
| Initialize SpaFmt structs.
|----------------------------------------------------------------------
| on entry:
|==========
| ( amat )  =  Pointer to a SpaFmt struct.
|     len   =  size of matrix
|     job   =  0: pattern only
|              1: data and pattern
|
| On return:
|===========
|
|  amat->n
|      ->*nzcount
|      ->**ja
|      ->**ma
|
| integer value returned:
|             0   --> successful return.
|             1   --> memory allocation error.
|--------------------------------------------------------------------*/
   amat->n = len;
   amat->nzcount = (int *)Malloc( len*sizeof(int), "setupCS" );
   amat->ja = (int **) Malloc( len*sizeof(int *), "setupCS" );
   if( job == 1 ) 
       amat->ma = (double **) Malloc( len*sizeof(double *), "setupCS" );
   else
       amat->ma = NULL;
   return 0;
}
/*---------------------------------------------------------------------
|     end of setupCS
|--------------------------------------------------------------------*/

int cleanCS(csptr amat)
{
/*----------------------------------------------------------------------
| Free up memory allocated for SpaFmt structs.
|----------------------------------------------------------------------
| on entry:
|==========
| ( amat )  =  Pointer to a SpaFmt struct.
|--------------------------------------------------------------------*/
   /*   */
  int i;
  if (amat == NULL) return 0;
  if (amat->n < 1) return 0;
  for (i=0; i<amat->n; i++) {
    if (amat->nzcount[i] > 0) {
      if( amat->ma ) free(amat->ma[i]);
      free(amat->ja[i]);
    }
  }    
  if (amat->ma) free(amat->ma);
  free(amat->ja);
  free(amat->nzcount);
  free(amat);
  return 0;
}
/*---------------------------------------------------------------------
|     end of cleanCS
|--------------------------------------------------------------------*/

int cscpy(csptr amat, csptr bmat){
/*----------------------------------------------------------------------
| Convert CSR matrix to SpaFmt struct
|----------------------------------------------------------------------
| on entry:
|==========
| ( amat )   = Matrix stored in SpaFmt format
|
|
| On return:
|===========
|
| ( bmat )  =  Matrix stored as SpaFmt struct containing a copy
|              of amat 
|
|       integer value returned:
|             0   --> successful return.
|             1   --> memory allocation error.
|--------------------------------------------------------------------*/
  int j, len, size=amat->n;
  double *bma;
  int *bja;
/*------------------------------------------------------------*/
  for (j=0; j<size; j++) {
    len = bmat->nzcount[j] = amat->nzcount[j];
    if (len > 0) {
      bja = (int *) Malloc(len*sizeof(int), "cscpy:1" );
      bma = (double *) Malloc(len*sizeof(double), "cscpy:2" );
      memcpy(bja,amat->ja[j],len*sizeof(int));
      memcpy(bma,amat->ma[j],len*sizeof(double));
      bmat->ja[j] = bja;
      bmat->ma[j] = bma;
    }	
  }
  return 0;
}
/*-----------------------------------------------------------------------*/

int setupILU( iluptr lu, int n )
{
/*----------------------------------------------------------------------
| Initialize ILUSpar structs.
|----------------------------------------------------------------------
| on entry:
|==========
|   ( lu )  =  Pointer to a ILUSpar struct.
|       n   =  size of matrix
|
| On return:
|===========
|
|    lu->n
|      ->L     L matrix, SpaFmt format
|      ->D     Diagonals
|      ->U     U matrix, SpaFmt format
|      ->work  working buffer of length n
|      ->bf    buffer
|
| integer value returned:
|             0   --> successful return.
|            -1   --> memory allocation error.
|--------------------------------------------------------------------*/
    lu->n  = n;
    lu->D = (double *)Malloc( sizeof(double) * n, "setupILU" );
    lu->L = (csptr)Malloc( sizeof(SparMat), "setupILU" );
    setupCS( lu->L, n, 1 );
    lu->U = (csptr)Malloc( sizeof(SparMat), "setupILU" );
    setupCS( lu->U, n, 1 );
    lu->work = (int *)Malloc( sizeof(int) * n, "setupILU" );
    return 0;
}
/*---------------------------------------------------------------------
|     end of setupILU
|--------------------------------------------------------------------*/


int cleanILU( iluptr lu )
{
/*----------------------------------------------------------------------
| Free up memory allocated for ILUSpar structs.
|----------------------------------------------------------------------
| on entry:
|==========
|   ( lu )  =  Pointer to a ILUSpar struct.
|--------------------------------------------------------------------*/
  if( NULL == lu ) return 0;
  if( lu->D ) {
    free( lu->D );
  }
  cleanCS( lu->L );
  cleanCS( lu->U );
  if( lu->work ) free( lu->work );
  free( lu );
  return 0;
}
/*---------------------------------------------------------------------
|     end of cleanILU
|--------------------------------------------------------------------*/

int setupVBMat( vbsptr vbmat, int n, int *nB )
{
/*----------------------------------------------------------------------
| Initialize VBSpaFmt structs.
|----------------------------------------------------------------------
| on entry:
|==========
| ( vbmat ) =  Pointer to a VBSpaFmt struct.
|       n   =  size of block matrix
|      nB   =  size of diagonal block, so the real size of the matrix
|              is nB[0] + nB[1] + ... + nB[n-1]
|              do nothing if nB is NULL
|
| On return:
|===========
|
| vbmat->n
|      ->*bsz
|      ->*nzcount
|      ->**ja
|      ->**ba
|      ->*D
|
| integer value returned:
|             0   --> successful return.
|            -1   --> memory allocation error.
|--------------------------------------------------------------------*/
    int i;
    vbmat->n = n;
    if( nB ) {
        vbmat->bsz = (int *)Malloc( sizeof(int)*(n+1), "setupVBMat" );
        vbmat->bsz[0] = 0;
        for( i = 1; i <= n; i++ ) {
            vbmat->bsz[i] = vbmat->bsz[i-1] + nB[i-1];
        }
    } else
        vbmat->bsz = NULL;
    vbmat->nzcount = (int *)Malloc( sizeof(int)*n, "setupVBMat" );
    vbmat->ja = (int **)Malloc( sizeof(int *)*n, "setupVBMat" );
    vbmat->ba = (BData **)Malloc( sizeof(BData *) * n, "setupVBMat" );
    vbmat->D = NULL;
    return 0;
}
/*---------------------------------------------------------------------
|     end of setupVBMat
|--------------------------------------------------------------------*/

int cleanVBMat( vbsptr vbmat )
{
/*----------------------------------------------------------------------
| Free up memory allocated for VBSpaFmt structs.
|----------------------------------------------------------------------
| on entry:
|==========
| ( vbmat )  =  Pointer to a VBSpaFmt struct.
|--------------------------------------------------------------------*/
    int i, j; 
    if( vbmat == NULL ) return 0;
    if( vbmat->n < 1 ) return 0;
    
    for( i = 0; i < vbmat->n; i++ ) {
        if( vbmat->nzcount[i] > 0 ) {
            free( vbmat->ja[i] );
            if( vbmat->ba && vbmat->ba[i] ) {
                for( j = 0; j < vbmat->nzcount[i]; j++ ) {
                    free( vbmat->ba[i][j] );
                }
                free( vbmat->ba[i] );
            }
        }
        if( vbmat->D && vbmat->D[i] ) free( vbmat->D[i] );
    }
    if( vbmat->D ) free( vbmat->D );
    free( vbmat->ja );
    if( vbmat->ba ) free( vbmat->ba );
    free( vbmat->nzcount );
    if( vbmat->bsz ) free( vbmat->bsz );
    free( vbmat );
    return 0;
}
/*---------------------------------------------------------------------
|     end of cleanVBMat
|--------------------------------------------------------------------*/

int nnzVBMat( vbsptr vbmat )
{
    int nnz = 0, i, n = vbmat->n;
    for( i = 0; i < n; i++ ) {
        nnz += vbmat->nzcount[i];
    }
    return nnz;
}

int memVBMat( vbsptr vbmat )
{
    int mem = 0, nnz, i, j, n = vbmat->n, *bsz = vbmat->bsz, dm;
    for( i = 0; i < n; i++ ) {
        nnz = vbmat->nzcount[i];
        dm = 0;
        for( j = 0; j < nnz; j++ ) {
            dm += B_DIM(bsz,vbmat->ja[i][j]);
        }
        mem += dm * B_DIM(bsz,i);
    }
    return mem;
}

int setupVBILU( vbiluptr lu, int n, int *bsz )
{
/*----------------------------------------------------------------------
| Initialize VBILUSpar structs.
|----------------------------------------------------------------------
| on entry:
|==========
|   ( lu )  =  Pointer to a VBILUSpar struct.
|       n   =  size of block matrix
|     bsz   =  the row/col of the first element of each diagonal block
|
| On return:
|===========
|
|    lu->n
|      ->bsz
|      ->L     L matrix, VBSpaFmt format
|      ->D     Diagonals
|      ->U     U matrix, VBSpaFmt format
|      ->work  working buffer of length n
|      ->bf    buffer
|
| integer value returned:
|             0   --> successful return.
|            -1   --> memory allocation error.
|--------------------------------------------------------------------*/
    int i;
    int max_block_size = sizeof(double)*MAX_BLOCK_SIZE*MAX_BLOCK_SIZE;
    lu->n  = n;
    lu->bsz = (int *)Malloc( sizeof(int) * (n+1), "setupVBILU" );
    for( i = 0; i <= n; i++ ) lu->bsz[i] = bsz[i];
    lu->D = (BData *)Malloc( sizeof(BData) * n, "setupVBILU" );
    lu->L = (vbsptr)Malloc( sizeof(VBSparMat), "setupVBILU" );
    setupVBMat( lu->L, n, NULL );
    lu->U = (vbsptr)Malloc( sizeof(VBSparMat), "setupVBILU" );
    setupVBMat( lu->U, n, NULL );
    lu->work = (int *)Malloc( sizeof(int) * n, "setupVBILU" );
    lu->bf = (BData)Malloc( max_block_size, "setupVBILU" );
    return 0;
}
/*---------------------------------------------------------------------
|     end of setupVBILU
|--------------------------------------------------------------------*/

int cleanVBILU( vbiluptr lu ){
/*----------------------------------------------------------------------
| Free up memory allocated for VBILUSpar structs.
|----------------------------------------------------------------------
| on entry:
|==========
|   ( lu )  =  Pointer to a VBILUSpar struct.
|--------------------------------------------------------------------*/
    int n = lu->n, i;
    if( NULL == lu ) return 0;
    if( lu->D ) {
        for( i = 0; i < n; i++ ) {
            if( lu->D[i] ) free( lu->D[i] );
        }
        free( lu->D );
    }
    if( lu->bsz ) free( lu->bsz );
    cleanVBMat( lu->L );
    cleanVBMat( lu->U );
    if( lu->work ) free( lu->work );
    if( lu->bf ) free( lu->bf );
    free( lu );
    return 0;
}
/*---------------------------------------------------------------------
|     end of cleanVBILU
|--------------------------------------------------------------------*/


int mallocRow( iluptr lu, int nrow )
{
/*----------------------------------------------------------------------
| Prepare space of a row according to the result of level structure
|----------------------------------------------------------------------
| on entry:
|==========
|   ( lu )  =  Pointer to a ILUSpar struct.
|     nrow  =  the current row to deal with
|
| On return:
|===========
|
|    lu->L->ma[nrow][...]
|      ->U->ma[nrow][...]
|
| integer value returned:
|             0   --> successful return.
|            -1   --> memory allocation error.
|--------------------------------------------------------------------*/
    int nzcount = lu->L->nzcount[nrow];
    lu->L->ma[nrow] = (double *)Malloc( sizeof(double)*nzcount, "mallocRow" );
    nzcount = lu->U->nzcount[nrow];
    lu->U->ma[nrow] = (double *)Malloc( sizeof(double)*nzcount, "mallocRow" );
    return 0;
}
/*---------------------------------------------------------------------
|     end of mallocVBRow
|--------------------------------------------------------------------*/

int mallocVBRow( vbiluptr lu, int nrow )
{
/*----------------------------------------------------------------------
| Prepare space of a row according to the result of level structure
|----------------------------------------------------------------------
| on entry:
|==========
|   ( lu )  =  Pointer to a VBILUSpar struct.
|     nrow  =  the current row to deal with
|
| On return:
|===========
|
|    lu->L->ba[nrow][...]
|      ->D[nrow]
|      ->U->ba[nrow][...]
|
| integer value returned:
|             0   --> successful return.
|            -1   --> memory allocation error.
|--------------------------------------------------------------------*/
    int j, nzcount, ncol, szOfBlock;
    int *bsz = lu->bsz;
    
    nzcount = lu->L->nzcount[nrow];
    lu->L->ba[nrow] = (BData *)Malloc( sizeof(BData)*nzcount, "mallocVBRow" );
    for( j = 0; j < nzcount; j++ ) {
        ncol = lu->L->ja[nrow][j];
        szOfBlock = B_DIM(bsz,nrow)*B_DIM(bsz,ncol)*sizeof(double);
        lu->L->ba[nrow][j] = (BData)Malloc( szOfBlock, "mallocVBRow" );
    }

    szOfBlock = sizeof(double) * B_DIM(bsz,nrow) * B_DIM(bsz,nrow);
    lu->D[nrow] = (BData)Malloc( szOfBlock, "mallocVBRow" );

    nzcount = lu->U->nzcount[nrow];
    lu->U->ba[nrow] = (BData *)Malloc( sizeof(BData)*nzcount, "mallocVBRow" );
    for( j = 0; j < nzcount; j++ ) {
        ncol = lu->U->ja[nrow][j];
        szOfBlock = B_DIM(bsz,nrow)*B_DIM(bsz,ncol)*sizeof(double);
        lu->U->ba[nrow][j] = (BData)Malloc( szOfBlock, "mallocVBRow" );
    }
    return 0;
}
/*---------------------------------------------------------------------
|     end of mallocVBRow
|--------------------------------------------------------------------*/

void zrmC( int m, int n, BData data ){
  int mn = m * n, i;
  for( i = 0; i < mn; i++ ) data[i] = 0;
}

void copyBData( int m, int n, BData dst, BData src, int isig ){
  int mm = m * n, i;
  if( isig == 0 )
    for( i = 0; i < mm; i++ ) dst[i] = src[i];
  else
    for( i = 0; i < mm; i++ ) dst[i] = -src[i];
}

int setupP4 (p4ptr amat, int Bn, int Cn,  csptr F,  csptr E) 
{
/*----------------------------------------------------------------------
| initialize PerMat4 struct given the F, E, blocks.  
|----------------------------------------------------------------------
| on entry:
|==========
| ( amat )  =  Pointer to a PerMat4 struct.
|     Bn    =  size of B block
|     Cn    =  size of C block
|     F, E  = the two blocks to be assigned to srtruct - without the
|
| On return:
|===========
|
|  amat->L                for each block: amat->M->n
|      ->U                                       ->nzcount
|      ->E                                       ->ja
|      ->F                                       ->ma
|      ->perm
|      ->rperm       (if meth[1] > 0)
|      ->D1          (if meth[2] > 0)
|      ->D2          (if meth[3] > 0)
|
|  Scaling arrays are initialized to 1.0.
|
|       integer value returned:
|             0   --> successful return.
|             1   --> memory allocation error.
|--------------------------------------------------------------------*/
   int n;
   /* size n */
   n = amat->n = Bn + Cn;
   amat->nB = Bn; 
/* amat->perm = (int *) Malloc(n*sizeof(int), "setupP4:1" ); */
/*   assign space for wk -- note that this is only done at 1st level
     at other levels, copy pointer of wk from previous level */
   if (amat->prev == NULL)  /* wk has 2 * n entries now */
     amat->wk   = (double *) Malloc(2*n*sizeof(double), "setupP4:2" );
   else 
     amat->wk = (amat->prev)->wk; 

/*-------------------- L and U */ 
   amat->L = (csptr) Malloc(sizeof(SparMat), "setupP4:3" );
   if (setupCS(amat->L, Bn,1)) return 1;
   /*    fprintf(stdout,"  -- BN %d   Cn   %d \n", Bn,Cn);  */
   amat->U = (csptr) Malloc(sizeof(SparMat), "setupP4:4" );
   if (setupCS(amat->U, Bn,1)) return 1;

   amat->F = F; 
   amat->E = E; 
   return 0;
}
/*---------------------------------------------------------------------
|     end of setupP4 
|--------------------------------------------------------------------*/

int cleanP4(p4ptr amat)
{
/*----------------------------------------------------------------------
| Free up memory allocated for Per4Mat structs.
|----------------------------------------------------------------------
| on entry:
|==========
| ( amat )  =  Pointer to a Per4Mat struct.
|--------------------------------------------------------------------*/

/*  -------------------------- */
  if (amat == NULL) return 0;
  if (amat->n < 1) return 0;
   

  if (amat->perm) {
    if (amat->perm) free(amat->perm); 
    amat->perm = NULL;
  }
  
  if (!amat->symperm) { 
    if (amat->rperm) free(amat->rperm); 
    amat->rperm = NULL;
  } 
  
  if (amat->F) {
    cleanCS(amat->F); 
    amat->F = NULL;
  }
  if (amat->E) {
    cleanCS(amat->E); 
    amat->E = NULL;
  }
  if (amat->L) {
    cleanCS(amat->L);
    amat->L = NULL;
   }
  if (amat->U) {
    cleanCS(amat->U);
    amat->U = NULL;
  }
  
  if (amat->prev == NULL) 
    if (amat->wk) free(amat->wk);  
  
  if (amat->D1) free(amat->D1);
  if (amat->D2) free(amat->D2);
  return 0;
}
/*---------------------------------------------------------------------
|     end of cleanP4
|--------------------------------------------------------------------*/


int setupILUT(ilutptr amat, int len)
{
/*----------------------------------------------------------------------
| Allocate pointers for ILUTfac structs.
|----------------------------------------------------------------------
| on entry:
|==========
| ( amat )  =  Pointer to a ILUTfac struct.
|     len   =  size of L U  blocks
|
| On return:
|===========
|
|  amat->L                for each block: amat->M->n
|      ->U                                       ->nzcount
|                                                ->ja
|                                                ->ma
|      ->rperm       (if meth[0] > 0)
|      ->perm2       (if meth[1] > 0)
|      ->D1          (if meth[2] > 0)
|      ->D2          (if meth[3] > 0)
|
|  Permutation arrays are initialized to the identity.
|  Scaling arrays are initialized to 1.0.
|
|       integer value returned:
|             0   --> successful return.
|             1   --> memory allocation error.
|--------------------------------------------------------------------*/
  amat->n = len;
 amat->wk = (double *) Malloc(2*len*sizeof(double), "setupILUT:5" );
  amat->L = (csptr) Malloc(sizeof(SparMat), "setupILUT:6" );
  if (setupCS(amat->L, len,1)) return 1;
  amat->U = (csptr) Malloc(sizeof(SparMat), "setupILUT:7" );
  if (setupCS(amat->U, len,1)) return 1;
  return 0;
}    
/*---------------------------------------------------------------------
|     end of setupILUT
|--------------------------------------------------------------------*/
int cleanILUT(ilutptr amat, int indic)
{
/*----------------------------------------------------------------------
| Free up memory allocated for IluSpar structs.
|----------------------------------------------------------------------
| on entry:
|==========
| ( amat )  =  Pointer to a IluSpar struct.
|  indic    = indicator for number of levels.  indic=0 -> zero level.
|--------------------------------------------------------------------*/

  /*----------------*/
   
  if (amat->wk) {
    free(amat->wk); 
    amat->wk = NULL;
  }
  cleanCS(amat->L);
  cleanCS(amat->U);

  if (indic) cleanCS(amat->C);  
/*-------------------- nonsymmetric permutation */
  if (amat->rperm) {
    free(amat->rperm);
    amat->rperm = NULL;
  }
  if (amat->perm) {
    free(amat->perm); 
    amat->perm = NULL;
  }
  
/*-------------------- ilutp permutation */
  if (amat->perm2) free(amat->perm2);
/*-------------------- diagonal scalings */
   if (amat->D1) free(amat->D1);
   if (amat->D2) free(amat->D2);
   return 0;
}
/*---------------------------------------------------------------------
|     end of cleanILUT
|--------------------------------------------------------------------*/


void setup_arms (arms Levmat) {
  Levmat->ilus = (ilutptr) Malloc(sizeof(IluSpar), "setup_arms:ilus" );
  Levmat->levmat = (p4ptr) Malloc(sizeof(Per4Mat), "setup_arms:levmat" );
}
int cleanARMS(arms ArmsPre)
{
  p4ptr amat = ArmsPre->levmat;
  ilutptr cmat = ArmsPre->ilus;
/*----------------------------------------------------------------------
| Free up memory allocated for entire ARMS preconditioner.
|----------------------------------------------------------------------
| on entry:
|==========
| ( amat )  =  Pointer to a Per4Mat struct.
| ( cmat )  =  Pointer to a IluSpar struct.
|--------------------------------------------------------------------*/
/* case when nlev == 0 */  
  int indic=(amat->nB != 0) ;
    /*  && amat->next !=NULL) ; */
  
  p4ptr levc, levn;

  levc = amat; 

  if (indic) { 
    while (levc) {
      if (cleanP4(levc)) return(1) ; 
      levn = levc->next;
      free(levc);
      levc = levn;
    }		
  }	
   else 	
     if (amat) {
       free(amat) ; 
       amat = NULL;
     }
  
  cleanILUT(cmat,indic); 
  
  
  if (cmat) {
    free(cmat);	
    cmat = NULL;
  }
  
  return 0;
}
/*---------------------------------------------------------------------
|     end of cleanARMS 
|--------------------------------------------------------------------*/


int csSplit4(csptr amat, int bsize, int csize, csptr B, csptr F,
	     csptr E, csptr C)
{
/*---------------------------------------------------------------------
| Convert permuted csrmat struct to PerMat4 struct 
|                - matrix already permuted
|----------------------------------------------------------------------
| on entry:
|========== 
| ( amat )  =  Matrix stored in SpaFmt format.
|              Internal pointers (and associated memory) destroyed before
|              return.
|
| On return:
|===========
|
| B, E, F, C = 4 blocks in 
| 
|          | B   F |      
|   Amat = |       | 
|          | E   C | 
| 
|
|       integer value returned:
|             0   --> successful return.
|             1   --> memory allocation error.
|--------------------------------------------------------------------*/
   int j, j1, numr, numl, ind, newj, rowz, *rowj, *new1j, *new2j;
   double *rowm, *new1m, *new2m;
/*---------------------------------------------------------------------
|     Sort the matrix and separate into   |  B  F  |
|                                         |        |
|                                         |  E  C  |
|--------------------------------------------------------------------*/
   if (setupCS(B,bsize,1)) goto label111; 
   if (setupCS(F,bsize,1)) goto label111;
   if (setupCS(E,csize,1)) goto label111;
   if (setupCS(C,csize,1)) goto label111;
   new1j = (int *) Malloc(bsize*sizeof(int), "csSplit4:1" );
   new2j = (int *) Malloc(csize*sizeof(int), "csSplit4:2" );
   new1m = (double *) Malloc(bsize*sizeof(double), "csSplit4:3" );
   new2m = (double *) Malloc(csize*sizeof(double), "csSplit4:4" );
/*    B and F blocks */ 
   for (j=0; j<bsize; j++) {
      numl = numr = 0;
      rowz = amat->nzcount[j];
      rowj = amat->ja[j];
      rowm = amat->ma[j];
      for (j1=0; j1<rowz; j1++) {
	 if (rowj[j1]<bsize) numl++;
	 else numr++;
      }
      B->nzcount[j] = numl;
      F->nzcount[j] = numr;
      if (numl>0) {
	 B->ja[j] = (int *) Malloc(numl*sizeof(int), "csSplit4:5" );
	 B->ma[j] = (double *) Malloc(numl*sizeof(double), "csSplit4:6" );
      }
      if (numr>0) {
	 F->ja[j] = (int *) Malloc(numr*sizeof(int), "csSplit4:7" );
	 F->ma[j] = (double *) Malloc(numr*sizeof(double), "csSplit4:8" );
      }
      numl = numr = 0;
      for (j1=0; j1<rowz; j1++) {
	 newj = rowj[j1];
	 if (newj<bsize) {
	    new1j[numl] = newj;
	    new1m[numl] = rowm[j1];
	    numl++;
	 }
	 else {
	    new2j[numr] = newj - bsize;
	    new2m[numr] = rowm[j1];
	    numr++;
	 }
      }
      memcpy(B->ja[j], new1j, numl*sizeof(int));
      memcpy(B->ma[j], new1m, numl*sizeof(double));
      memcpy(F->ja[j], new2j, numr*sizeof(int));
      memcpy(F->ma[j], new2m, numr*sizeof(double));
   }
/*    E and C blocks */
   for (j=0; j<csize; j++) {
      numl = numr = 0;
      ind = bsize + j;
      rowz = amat->nzcount[ind];
      rowj = amat->ja[ind];
      rowm = amat->ma[ind];
      for (j1=0; j1<rowz; j1++) {
	 if (rowj[j1]<bsize) numl++;
	 else numr++;
      }
      E->nzcount[j] = numl;
      C->nzcount[j] = numr;
      if (numl>0) {
	E->ja[j] = (int *) Malloc(numl*sizeof(int), "csSplit4:9" );
	E->ma[j] = (double *) Malloc(numl*sizeof(double), "csSplit4:10" );
      }	
      if (numr>0) {
	C->ja[j] = (int *) Malloc(numr*sizeof(int), "csSplit4:11" );
	C->ma[j] = (double *) Malloc(numr*sizeof(double), "csSplit4:12" );
      }		
      numl = numr = 0;
      for (j1=0; j1<rowz; j1++) {
	newj = rowj[j1];
	if (newj<bsize) {
	  new1j[numl] = newj;
	  new1m[numl] = rowm[j1];
	  numl++;
	}
	else {
	  new2j[numr] = newj - bsize;
	  new2m[numr] = rowm[j1];
	  numr++;
	}
      }
      memcpy(E->ja[j], new1j, numl*sizeof(int));
      memcpy(E->ma[j], new1m, numl*sizeof(double));
      memcpy(C->ja[j], new2j, numr*sizeof(int));
      memcpy(C->ma[j], new2m, numr*sizeof(double));
   }

   if (new1j) free(new1j);
   if (new2j) free(new2j);
   if (new1m) free(new1m);
   if (new2m) free(new2m);
   return 0;
label111:
   return 1;
}
/*---------------------------------------------------------------------
|     end of csSplit4
|--------------------------------------------------------------------*/


int CSRcs( int n, double *a, int *ja, int *ia, csptr mat, int rsa )
{
/*----------------------------------------------------------------------
| Convert CSR matrix to SpaFmt struct
|----------------------------------------------------------------------
| on entry:
|==========
| a, ja, ia  = Matrix stored in CSR format (with FORTRAN indexing).
| rsa        = source file is symmetric HB matrix 
|
| On return:
|===========
|
| ( mat )  =  Matrix stored as SpaFmt struct. (C indexing)
|
|       integer value returned:
|             0   --> successful return.
|             1   --> memory allocation error.
|--------------------------------------------------------------------*/
  int i, j, j1, len, col, nnz;
  double *bra;
  int *bja;
  /*    setup data structure for mat (csptr) struct */
  setupCS( mat, n, 1 );

  if( rsa ) { /* RSA HB matrix */
    for( j = 0; j < n; j++ ) {
      len = ia[j+1] - ia[j];
      mat->nzcount[j] = len;
    }
    for( j = 0; j < n; j++ ) {
      for( j1 = ia[j]-1; j1 < ia[j+1]-1; j1++ ) {
        col = ja[j1] - 1;
        if( col != j ) mat->nzcount[col]++;
      }
    }
    for( j = 0; j < n; j++ ) {
      nnz = mat->nzcount[j];
      mat->ja[j] = (int *)Malloc( nnz * sizeof(int), "CSRcs" );
      mat->ma[j] = (double *)Malloc( nnz * sizeof(double), "CSRcs" );
      mat->nzcount[j] = 0;
    }
    for( j = 0; j < n; j++ ) {
      for( j1 = ia[j]-1; j1 < ia[j+1]-1; j1++ ) {
        col = ja[j1] - 1;
        mat->ja[j][mat->nzcount[j]] = col;
        mat->ma[j][mat->nzcount[j]] = a[j1];
        mat->nzcount[j]++;
        if( col != j ) {
          mat->ja[col][mat->nzcount[col]] = j;
          mat->ma[col][mat->nzcount[col]] = a[j1];
          mat->nzcount[col]++;
        }
      }
    }
    return 0;
  }

  for (j=0; j<n; j++) {
    len = ia[j+1] - ia[j];
    mat->nzcount[j] = len;
    if (len > 0) {
      bja = (int *) Malloc( len*sizeof(int), "CSRcs" );
      bra = (double *) Malloc( len*sizeof(double), "CSRcs" );
      i = 0;
      for (j1=ia[j]-1; j1<ia[j+1]-1; j1++) {
        bja[i] = ja[j1] - 1;
        bra[i] = a[j1] ;
        i++;
      }
      mat->ja[j] = bja;
      mat->ma[j] = bra;
    }
  }    
  return 0;
}
/*---------------------------------------------------------------------
|     end of CSRcs
|--------------------------------------------------------------------*/


int COOcs(int n, int nnz,  double *a, int *ja, int *ia, csptr bmat)
{
/*----------------------------------------------------------------------
| Convert COO matrix to SpaFmt struct
|----------------------------------------------------------------------
| on entry:
|==========
| a, ja, ia  = Matrix stored in COO format  -- a =  entries
|                                             ja = column indices
|                                             ia = row indices 
| On return:
|===========
|
| ( bmat )  =  Matrix stored as SpaFmt struct.
|
|       integer value returned:
|             0   --> successful return.
|             1   --> memory allocation error.
|--------------------------------------------------------------------*/
  int i, k, k1, l, job = 1;
  int *len;
  int setupCS(csptr, int, int);
  /*-------------------- setup data structure for bmat (csptr) struct */
  if (setupCS(bmat, n, job)) {
    printf(" ERROR SETTING UP bmat IN SETUPCS \n") ;
    exit(0);
  }
/*-------------------- determine lengths */
   len =  (int*) Malloc(n*sizeof(int), "COOcs:0" );
   for (k=0; k<n; k++) 
     len[k] = 0;
   for (k=0; k<nnz; k++) 
     ++len[ia[k]]; 
/*-------------------- allocate          */
   for (k=0; k<n; k++) {
     l = len[k];
     bmat->nzcount[k] = l;
     if (l > 0) {
       bmat->ja[k] = (int *) Malloc(l*sizeof(int), "COOcs:1" );
       bmat->ma[k] = ( double *) Malloc(l*sizeof( double),"COOcs:2" );
     }
     len[k] = 0;
   }
/*-------------------- Fill actual entries */
   for (k=0; k<nnz; k++) {
     i  = ia[k];
     k1 = len[i];
     (bmat->ja[i])[k1] = ja[k] ; 
     (bmat->ma[i])[k1] =  a[k] ;  
     len[i]++; 
   }
   free(len);
   return 0;
}
/*---------------------------------------------------------------------
|     end of COOcs
|--------------------------------------------------------------------*/


void coocsc(int n, int nnz, double *val, int *col, int *row, 
            double **a, int **ja, int **ia, int job)
/* job == 0, input coo in 0-indexing
   job == 1, input coo in 1-indexing
*/
{
  *a = (double *)Malloc(nnz*sizeof(double), "coocsc");
  *ja = (int *)Malloc(nnz*sizeof(int), "coocsc");
  *ia = (int *)Malloc((n+1)*sizeof(int), "coocsc");

  int i, *ir, *jc;

  if (job == 0) {
    ir = (int *)Malloc(nnz*sizeof(int), "coocsc");
    jc = (int *)Malloc(nnz*sizeof(int), "coocsc");
    for (i=0; i<nnz; i++) {
      ir[i] = row[i] + 1;
      jc[i] = col[i] + 1;
    }
  }
  else {
    ir = row;
    jc = col;
  }

  coocsr_(&n, &nnz, val, jc, ir, *a, *ja, *ia);

  if (job == 0) {
    free(ir);
    free(jc);
  }
}


int csrvbsrC( int job, int nBlk, int *nB, csptr csmat, vbsptr vbmat )
{
/*----------------------------------------------------------------------
 *  Compressed C-style Sparse Row to C-style Various Block Sparse Row
 *----------------------------------------------------------------------
 *
 * This  subroutine converts a matrix stored  in a C-style SpaFmt format
 * into a C-style various block SpaFmt format
 *
 * NOTE: the initial matrix does not have to have a block structure. 
 * zero padding is done for general sparse matrices. 
 *
 *----------------------------------------------------------------------
 * on entry:
 *----------
 * job   = if job == 0 on entry, pattern only is generated. 
 *
 * nBlk  = integer equal to the dimension of block matrix.
 * 
 * nB    = integer array of diagonals' block size
 *
 * csmat = Sparse Row format Matrix
 *
 * on return:
 *-----------
 * 
 * vbmat = Various Block Sparse Row format Matrix
 *
 * ierr  = integer, error code. 
 *              0  -- normal termination
 *             -1  -- error occur
 *
 *---------------------------------------------------------------------*/
    int n, i, j, k;
    int nnz, szofBlock, ipos, b_row, b_col, br, bc;
    int *iw = NULL;

    n  = csmat->n;            /* size of the original matrix          */
    setupVBMat( vbmat, nBlk, nB );
    iw = (int *)Malloc( sizeof(int)*nBlk, "csrvbsrC_1" );
    for( i = 0; i < nBlk; i++ ) iw[i] = 0;
    b_row = -1;
    for( i = 0; i < n; i += nB[b_row] ) {
        vbmat->nzcount[++b_row] = 0;

        /* calculate nzcount of the (b_row)-th row of the block matrix */
        for( j = i; j < i+nB[b_row]; j++ ) {
            int nnz_j = csmat->nzcount[j];
            for( k = 0; k < nnz_j; k++ ) {
                /* get the column ID of block matrix by giving the column ID
                   of the original matrix */
                b_col = col2vbcol( csmat->ja[j][k], vbmat );
                if( iw[b_col] == 0 ) {
                    iw[b_col] = 1;
                    vbmat->nzcount[b_row]++;
                }
            }
        }
        if( 0 == ( nnz = vbmat->nzcount[b_row] ) ) continue;
        vbmat->ja[b_row] = (int *)Malloc( sizeof(int)*nnz, "csrvbsrC_2" );

        /* calculate the pattern of the (b_row)-th row of the block matrix */
        for( j = 0, ipos = 0; j < nBlk; j++ ) {
            if( iw[j] != 0 ) {
                vbmat->ja[b_row][ipos] = j;
                iw[j] = ipos;
                ipos++;
            }
        }
        if( job == 0 ) goto NEXT_ROW;  /* stop here if patterns only */

        /* copy data to the (b_row)-th row of the block matrix from the
           original matrix */
        vbmat->ba[b_row] = (BData *)Malloc( sizeof(BData)*nnz, "csrvbsrC_3" );
        for( j = 0; j < nnz; j++ ) {
            szofBlock = sizeof(double)*nB[b_row]*nB[vbmat->ja[b_row][j]];
            vbmat->ba[b_row][j] = (BData)Malloc( szofBlock, "csrvbsrC_4" );
            memset( vbmat->ba[b_row][j], 0, szofBlock );
        }
        for( j = i; j < i+nB[b_row]; j++ ) {
            for( k = 0; k < csmat->nzcount[j]; k++ ) {
                /* get the column ID of block matrix by giving the column ID
                   of the original matrix */
                b_col = col2vbcol( csmat->ja[j][k], vbmat );
                ipos = iw[b_col];
                br = j - i;
                bc = csmat->ja[j][k] - vbmat->bsz[b_col];
                DATA(vbmat->ba[b_row][ipos],nB[b_row],br,bc) = csmat->ma[j][k];
            }
        }
NEXT_ROW:
        /* reset iw */
        for( j = 0; j < nnz; j++ ) iw[vbmat->ja[b_row][j]] = 0;
    }
    
    free( iw );
    return 0;
}

int col2vbcol( int col, vbsptr vbmat )
{
/*---------------------------------------------------------------------
 * get the column ID of block matrix by giving the column ID of the original
 * matrix
 *--------------------------------------------------------------------*/
    int *bsz = vbmat->bsz, n = vbmat->n;
    int begin = 0, mid, end = n-1;
    while( end - begin > 1 ) {
        mid = (begin+end)/2;
        if( col < bsz[mid] ) {
            end = mid;
        } else if( col >= bsz[mid+1] ) {
            begin = mid;
        } else {
            return mid;
        }
    }
    if( col >= bsz[end] ) {
        return end;
    }
    return begin;
}

int nnz_vbilu(vbiluptr lu )
{  
  /*-------------------- counts number of nonzero entries in preconditioner */
    int *bsz = lu->bsz;
    int nzcount, nnz = 0, i, j, col;
    for( i = 0; i < lu->n; i++ ) {
        nzcount = 0;
        for( j = 0; j < lu->L->nzcount[i]; j++ ) {
            col = lu->L->ja[i][j];
            nzcount += B_DIM(bsz,col);
        }
        for( j = 0; j < lu->U->nzcount[i]; j++ ) {
            col = lu->U->ja[i][j];
            nzcount += B_DIM(bsz,col);
        }
        nzcount += B_DIM(bsz,i);  /* diagonal */
        nzcount *= B_DIM(bsz,i);
        nnz += nzcount;
    }
    return nnz;
}


int nnz_ilu( iluptr lu )
{
  int nnz = 0, i;
  for( i = 0; i < lu->n; i++ ) {
    nnz += lu->L->nzcount[i];
    nnz += lu->U->nzcount[i];
    nnz++;
  }
  return nnz;
}


int nnz_lev4(p4ptr levmat, int *lev, FILE *ft) 
{
/* counts all nonzero elements in levmat struct  -- 
     recursive */
  int nnzT, nnzL, nnzU, nnzF, nnzE, nnzDown=0;
  p4ptr nextmat; 
  nnzL = nnz_cs(levmat->L); 
  nnzU = nnz_cs(levmat->U); 
  nnzF = nnz_cs(levmat->F); 
  nnzE = nnz_cs(levmat->E); 
  nnzT = nnzL+nnzU+nnzF+nnzE;
/*-------------------- print */
  if (ft) {
  if (*lev == 0) 
 fprintf(ft,"\nnnz/lev used:      L        U        F        E    subtot\n");  
 fprintf(ft,"    Level %2d %8d %8d %8d %8d %8d\n",*lev, nnzL,
	 nnzU,nnzF,nnzE,nnzT);
  }
  (*lev)++;
  nextmat = levmat->next; 
  if (nextmat != NULL) 
   nnzDown = nnz_lev4(nextmat, lev, ft);
   return (nnzT+nnzDown); 
}
  
int nnz_cs (csptr A) {
/*-------------------- counts number of nonzeros in CSR matrix A */
  int i, n=A->n, nnz=0; 
  for (i=0; i<n; i++) 
    nnz +=A->nzcount[i];
  return nnz;
  }

int nnz_arms (arms PreSt,  FILE *ft)
{ 
/*-------------------------------------------------------
| computes and prints out total number of nonzero elements
| used in ARMS factorization 
+--------------------------------------------------------*/
  p4ptr levmat   = PreSt->levmat; 
  ilutptr ilschu = PreSt->ilus; 
  int nlev       = PreSt->nlev;
  int ilev=0,nnz_lev,nnz_sch,nnz_tot; 
  nnz_lev = 0; 
  if (nlev) nnz_lev+= nnz_lev4(levmat, &ilev, ft);
  nnz_sch = nnz_cs(ilschu->L)+nnz_cs(ilschu->U);
  if (nlev) nnz_sch += nnz_cs(ilschu->C);
  nnz_tot = nnz_lev+nnz_sch; 
  if (ft){
    fprintf(ft,"\n"); 
    fprintf(ft,"Total nonzeros for interm. blocks.... =  %10d\n",nnz_lev);
    fprintf(ft,"Total nonzeros for last level ....... =  %10d\n",nnz_sch);
    fprintf(ft,"Grand total.......................... =  %10d\n",nnz_tot);
  }
  return nnz_tot;
}

int CSClum( int n, double *a, int *ja, int *ia, iluptr mat, int rsa )
{
/*----------------------------------------------------------------------
| Convert CSC matrix to LUSparMat struct
|----------------------------------------------------------------------
| on entry:
|==========
| a, ja, ia  = Matrix stored in CSC format (with FORTRAN indexing).
| rsa        = 0: rua matrix
|              1: source file is symmetric HB matrix 
|              2: pattern symmetrization (add zeros)
|
| On return:
|===========
|
| ( mat )  =  Matrix stored as LUSparMat struct.
|
|       integer value returned:
|             0   --> successful return.
|             1   --> memory allocation error.
|--------------------------------------------------------------------*/
    int row, col, i, k, j1, j2, nnz, id;
    double val, *D;
    csptr L, U;
    setupILU(mat, n); 

    L = mat->L;
    U = mat->U;
    D = mat->D;

    for( i = 0; i < n; i++ ) {
        L->nzcount[i] = 0;
        U->nzcount[i] = 0;
    }
    for( col = 0; col < n; col++ ) {
        j1 = ia[col] - 1;
        j2 = ia[col+1] - 1;
        for( k = j1; k < j2; k++ ) {
            row = ja[k] - 1;
            if( row > col ) {
                L->nzcount[col]++;
                if( rsa == 1 ) U->nzcount[col]++;
            } else if ( row == col ) {
                U->nzcount[row]++;
            } else {
                U->nzcount[row]++;
                if( rsa == 1 ) L->nzcount[row]++;
            }
        }
    }
    for( i = 0; i < n; i++ ) {
        if( rsa == 2 ) {
            /* nzcount(A + A^T) <= nzcount(A) + nnzcol(A^T) */
            nnz = L->nzcount[i] + U->nzcount[i];
            L->ja[i] = (int *)Malloc( nnz*sizeof(int), "CSClum1" );
            L->ma[i] = (double *)Malloc( nnz*sizeof(double), "CSClum2" );
            L->nzcount[i] = 0;
            U->ja[i] = (int *)Malloc( nnz*sizeof(int), "CSClum3" );
            U->ma[i] = (double *)Malloc( nnz*sizeof(double), "CSClum4" );
            U->nzcount[i] = 0;
        } else {
            nnz = L->nzcount[i];
            L->ja[i] = (int *)Malloc( nnz*sizeof(int), "CSClum4" );
            L->ma[i] = (double *)Malloc( nnz*sizeof(double), "CSClum5" );
            L->nzcount[i] = 0;
            nnz = U->nzcount[i];
            U->ja[i] = (int *)Malloc( nnz*sizeof(int), "CSClum6" );
            U->ma[i] = (double *)Malloc( nnz*sizeof(double), "CSClum7" );
            U->nzcount[i] = 0;
        }
    }
    for( col = 0; col < n; col++ ) {
        j1 = ia[col] - 1;
        j2 = ia[col+1] - 1;
	/*-------------------- */
        for( k = j1; k < j2; k++ ) {
            row = ja[k] - 1;
            val = a[k];
            if( row > col ) {
                id = L->nzcount[col];
                L->ja[col][id] = row;
                L->ma[col][id] = val;
                L->nzcount[col]++;
                if( rsa == 1 ) {
                    id = U->nzcount[col];
                    U->ja[col][id] = row;
                    U->ma[col][id] = val;
                    U->nzcount[col]++;
                }
            } else if ( row == col ) {
                id = U->nzcount[row];
                U->ja[row][id] = col;
                U->ma[row][id] = val;
                U->nzcount[row]++;
                D[col] = val;
            } else {
                id = U->nzcount[row];
                U->ja[row][id] = col;
                U->ma[row][id] = val;
                U->nzcount[row]++;
                if( rsa == 1 ) {
                    id = L->nzcount[row];
                    L->ja[row][id] = col;
                    L->ma[row][id] = val;
                    L->nzcount[row]++;
                }
            }
        }
    }
    if( rsa == 2 ) {
/* add zeros to make pattern symmetrization */
        int *idU, *idL;
        int nzcount, nnzcol, nnzU, nnzL, j;
        idU = (int *)Malloc( n*sizeof(int), "CSClum6" );
        idL = (int *)Malloc( n*sizeof(int), "CSClum7" );
        for( i = 0; i < n; i++ ) {
            idU[i] = 0;
            idL[i] = 0;
        }
        for( i = 0; i < n; i++ ) {
            nzcount = U->nzcount[i];
            for( j = 0; j < nzcount; j++ ) idU[U->ja[i][j]] = 1;
            nnzcol = L->nzcount[i];
            for( j = 0; j < nnzcol; j++ ) idL[L->ja[i][j]] = 1;
            nnzU = nzcount;
            for( j = 0; j < nnzcol; j++ ) {
                row = L->ja[i][j];
                if( idU[row] == 0 ) {
                    U->ja[i][nnzU] = row;
                    U->ma[i][nnzU] = 0.0;
                    nnzU++;
                }
            }
            nnzL = nnzcol;
            for( j = 0; j < nzcount; j++ ) {
                col = U->ja[i][j];
                if( col == i ) continue;
                if( idL[col] == 0 ) {
                    L->ja[i][nnzL] = col;
                    L->ma[i][nnzL] = 0.0;
                    nnzL++;
                }
            }
            for( j = 0; j < nzcount; j++ ) idU[U->ja[i][j]] = 0;
            for( j = 0; j < nnzcol; j++ ) idL[L->ja[i][j]] = 0;
            if( nnzU - nnzL != 1 ) {
                fprintf( stderr, "error in pattern symmetrization...\n" );
                exit(-1);
            }
            L->nzcount[i] = nnzL;
            U->nzcount[i] = nnzU;
        }
        free( idU );
        free( idL );
    }

    return 0;
}
/*---------------------------------------------------------------------
|     end of CSClum
|--------------------------------------------------------------------*/

int CSClumC(csptr amat, iluptr mat, int rsa )
{
/*----------------------------------------------------------------------
| Convert cs matrix to LUSparMat struct
|----------------------------------------------------------------------
| on entry:  
|==========
| amat      = matrix in cs format 
|
| On return:
|===========
|
| ( mat )  =  Matrix stored as LUSparMat struct.
|
|       integer value returned:
|             0   --> successful return.
|             1   --> memory allocation error.
|--------------------------------------------------------------------*/
    int row, col, i, k, n, nz, nnz, id, *ja;
    double val, *D, *ma;
    csptr L, U;
/*-------------------- begin */
    n = amat->n;
    setupILU(mat, n); 
    L = mat->L;
    U = mat->U;
    D = mat->D;

    for( i = 0; i < n; i++ ) {
        L->nzcount[i] = 0;
        U->nzcount[i] = 0;
    }
    for( col = 0; col < n; col++ ) {
      nz = amat->nzcount[col];
      ja = amat->ja[col];
      for( k = 0; k<nz; k++) { 
	row = ja[k];
	if( row > col ) {
	  L->nzcount[col]++;
	  if( rsa == 1 ) U->nzcount[col]++;
	} else if ( row == col ) {
	  U->nzcount[row]++;
	} else {
	  U->nzcount[row]++;
	  if( rsa == 1 ) L->nzcount[row]++;
	}
      }
    }
    for( i = 0; i < n; i++ ) {
        if( rsa == 2 ) {
            /* nzcount(A + A^T) <= nzcount(A) + nnzcol(A^T) */
            nnz = L->nzcount[i] + U->nzcount[i];
            L->ja[i] = (int *)Malloc( nnz*sizeof(int), "CSClum1" );
            L->ma[i] = (double *)Malloc( nnz*sizeof(double), "CSClum2" );
            L->nzcount[i] = 0;
            U->ja[i] = (int *)Malloc( nnz*sizeof(int), "CSClum3" );
            U->ma[i] = (double *)Malloc( nnz*sizeof(double), "CSClum4" );
            U->nzcount[i] = 0;
        } else {
	  nnz = L->nzcount[i];
	  L->ja[i] = (int *)Malloc( nnz*sizeof(int), "CSClum4" );
	  L->ma[i] = (double *)Malloc( nnz*sizeof(double), "CSClum5" );
	  L->nzcount[i] = 0;
	  nnz = U->nzcount[i];
	  U->ja[i] = (int *)Malloc( nnz*sizeof(int), "CSClum6" );
	  U->ma[i] = (double *)Malloc( nnz*sizeof(double), "CSClum7" );
	  U->nzcount[i] = 0;
        }
    }
    for( col = 0; col < n; col++ ) {
      nz = amat->nzcount[col];
      ja = amat->ja[col];
      ma = amat->ma[col];
      for( k = 0;  k < nz; k++ ) {
            row = ja[k];
            val = ma[k];
            if( row > col ) {
	      id = L->nzcount[col];
	      L->ja[col][id] = row;
	      L->ma[col][id] = val;
	      L->nzcount[col]++;
	      if( rsa == 1 ) {
		id = U->nzcount[col];
		U->ja[col][id] = row;
		U->ma[col][id] = val;
		U->nzcount[col]++;
	      }
            } else if ( row == col ) {
	      id = U->nzcount[row];
	      U->ja[row][id] = col;
	      U->ma[row][id] = val;
	      U->nzcount[row]++;
	      D[col] = val;
            } else {
	      id = U->nzcount[row];
	      U->ja[row][id] = col;
	      U->ma[row][id] = val;
	      U->nzcount[row]++;
	      if( rsa == 1 ) {
		id = L->nzcount[row];
		L->ja[row][id] = col;
		L->ma[row][id] = val;
		L->nzcount[row]++;
	      }
            }
      }
    }
    if( rsa == 2 ) {
      /* add zeros to make pattern symmetrization */
      int *idU, *idL;
      int nzcount, nnzcol, nnzU, nnzL, j;
      idU = (int *)Malloc( n*sizeof(int), "CSClum6" );
      idL = (int *)Malloc( n*sizeof(int), "CSClum7" );
      for( i = 0; i < n; i++ ) {
	idU[i] = 0;
	idL[i] = 0;
      }
      for( i = 0; i < n; i++ ) {
	nzcount = U->nzcount[i];
	for( j = 0; j < nzcount; j++ ) idU[U->ja[i][j]] = 1;
	nnzcol = L->nzcount[i];
	for( j = 0; j < nnzcol; j++ ) idL[L->ja[i][j]] = 1;
	nnzU = nzcount;
	for( j = 0; j < nnzcol; j++ ) {
	  row = L->ja[i][j];
	  if( idU[row] == 0 ) {
	    U->ja[i][nnzU] = row;
	    U->ma[i][nnzU] = 0.0;
	    nnzU++;
	  }
	}
	nnzL = nnzcol;
	for( j = 0; j < nzcount; j++ ) {
	  col = U->ja[i][j];
	  if( col == i ) continue;
	  if( idL[col] == 0 ) {
	    L->ja[i][nnzL] = col;
	    L->ma[i][nnzL] = 0.0;
	    nnzL++;
	  }
	}
	for( j = 0; j < nzcount; j++ ) idU[U->ja[i][j]] = 0;
	for( j = 0; j < nnzcol; j++ ) idL[L->ja[i][j]] = 0;
	if( nnzU - nnzL != 1 ) {
	  fprintf( stderr, "error in pattern symmetrization...\n" );
	  exit(-1);
	}
	L->nzcount[i] = nnzL;
	U->nzcount[i] = nnzU;
      }
      free( idU );
      free( idL );
    }
    return 0;
}
/*---------------------------------------------------------------------
|     end of CSClum
|--------------------------------------------------------------------*/
