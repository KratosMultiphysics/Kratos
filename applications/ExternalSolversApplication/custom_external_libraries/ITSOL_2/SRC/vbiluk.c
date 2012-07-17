#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "globheads.h"

#ifndef min
#define min(a,b) (((a)>(b))?(b):(a))
#endif
#ifndef max
#define max(a,b) (((a)>(b))?(a):(b))
#endif
#define SVD 1
#define dgemm dgemm_

/*-------------------- protos */
void *Malloc(int nbytes, char *msg); 
int vblusolC(double *y, double *x, vbiluptr lu); 
int invGauss(int nn, double *A); 
int invSVD(int nn, double *A) ;
void dgemm(char*, char*, int*, int*, int*, double*, double*, int*, 
	   double*, int*, double*, double*, int*) ; 
int setupVBMat(vbsptr vbmat, int n, int *nB);
int mallocVBRow(vbiluptr lu, int nrow); 
void zrmC(int m, int n, BData data); 
int lofC( int lofM, vbsptr vbmat, vbiluptr lu, FILE *fp ); 
int setupVBILU(vbiluptr lu, int n, int *bsz);
void copyBData(int m, int n, BData dst, BData src, int isig);
int lofC( int lofM, vbsptr vbmat, vbiluptr lu, FILE *fp ); 
/*-------------------- END of protos */


int vbilukC( int lofM, vbsptr vbmat, vbiluptr lu, FILE *fp )
{
/*----------------------------------------------------------------------------
 * Block ILUK preconditioner
 * Block incomplete LU factorization with level of fill dropping
 * This version uses svd to invert diagonal blocks
 *----------------------------------------------------------------------------
 * Parameters
 *----------------------------------------------------------------------------
 * on entry:
 * =========
 * lofM     = level of fill: all entries with level of fill > lofM are
 *            dropped. Setting lofM = 0 gives BILU(0).
 * vbmat    = block matrix stored in VBSpaFmt format -- see globheads.h for
 *            details on format, the block sizes might be different
 * lu       = pointer to a VBILUSpar struct -- see globheads.h for details
 *            on format
 * fp       = file pointer for error log ( might be stderr )
 *
 * on return:
 * ==========
 * ierr     = return value.
 *            ierr  = 0   --> successful return.
 *            ierr  = -1  --> error in lofC
 *            ierr  = -2  --> singular diagonal block
 * lu->n    = dimension of the block matrix
 *   ->bsz  = the row/col of the first element of each diagonal block
 *            the size of the i-th row block should be bsz[i+1] - bsz[i]
 *   ->L    = L part -- stored in VBSpaFmt format
 *   ->D    = Diagonals
 *   ->U    = U part -- stored in VBSpaFmt format
 *----------------------------------------------------------------------------
 * Notes:
 * ======
 * All the diagonal blocks of the input block matrix must not be singular
 *--------------------------------------------------------------------------*/
    int ierr;
    int n = vbmat->n, *bsz = vbmat->bsz;
    int *jw, i, j, k, col, jpos, jrow, dim, sz;
    int mm, nn, kk;
    double alpha1 = 1.0, beta1 = 0.0, alpha2 = -1.0, beta2 = 1.0;
    vbsptr L, U;

    setupVBILU( lu, n, bsz );
    L = lu->L;
    U = lu->U;

    /* symbolic factorization to calculate level of fill index arrays */
    if( ( ierr = lofC( lofM, vbmat, lu, fp ) ) != 0 ) {
      fprintf( fp, "Error: lofC\n" );
      return -1;
    }

    jw = lu->work;
    /* set indicator array jw to -1 */
    for( j = 0; j < n; j++ ) jw[j] = -1;

    /* beginning of main loop */
    for( i = 0; i < n; i++ ) {
        dim = B_DIM(bsz,i);  /* number of rows of blocks in i-th row */
        /* set up the i-th row accroding to the nonzero information from
           symbolic factorization */
        mallocVBRow( lu, i );

        /* setup array jw[], and initial i-th row */
        for( j = 0; j < L->nzcount[i]; j++ ) {  /* initialize L part   */
            col = L->ja[i][j];
            sz = B_DIM(bsz,col);
            jw[col] = j;
            zrmC( dim, sz, L->ba[i][j] );
        }
        jw[i] = i;
        zrmC( dim, dim, lu->D[i] );            /* initialize diagonal */
        for( j = 0; j < U->nzcount[i]; j++ ) {  /* initialize U part   */
            col = U->ja[i][j];
            sz = B_DIM(bsz,col);
            jw[col] = j;
            zrmC( dim, sz, U->ba[i][j] );
        }

        /* copy row from vbmat into lu */
        for( j = 0; j < vbmat->nzcount[i]; j++ ) {
            col = vbmat->ja[i][j];
            sz = B_DIM(bsz,col);  /* number of columns of current block */
            jpos = jw[col];
            if( col < i ) {
                copyBData( dim, sz, L->ba[i][jpos], vbmat->ba[i][j], 0 );
            } else if( col == i ) {
                copyBData( dim, sz, lu->D[i], vbmat->ba[i][j], 0 );
            } else {
                copyBData( dim, sz, U->ba[i][jpos], vbmat->ba[i][j], 0 );
            }
        }

        /* eliminate previous rows */
        for( j = 0; j < L->nzcount[i]; j++ ) {
            jrow = L->ja[i][j];
            mm = dim;              /* number of rows of current block */
            nn = B_DIM(bsz,jrow);  /* number of cols of current block */
            /* get the multiplier for row to be eliminated (jrow) */
            dgemm( "n", "n", &mm, &nn, &nn, &alpha1, L->ba[i][j], &mm,
                    lu->D[jrow], &nn, &beta1, lu->bf, &mm );
            copyBData( mm, nn, L->ba[i][j], lu->bf, 0 );

            /* combine current row and row jrow */
            for( k = 0; k < U->nzcount[jrow]; k++ ) {
                col = U->ja[jrow][k];
                jpos = jw[col];
                if( jpos == -1 ) continue;
                if( col < i ) {
                    kk = B_DIM(bsz,col);
                    dgemm( "n", "n", &mm, &kk, &nn, &alpha2, L->ba[i][j],
                           &mm, U->ba[jrow][k], &nn, &beta2,
                           L->ba[i][jpos], &mm );
                } else if( col == i ) {
                    dgemm( "n", "n", &mm, &mm, &nn, &alpha2, L->ba[i][j],
                           &mm, U->ba[jrow][k], &nn, &beta2,
                           lu->D[i], &mm );
                } else {
                    kk = B_DIM(bsz,col);
                    dgemm( "n", "n", &mm, &kk, &nn, &alpha2, L->ba[i][j],
                           &mm, U->ba[jrow][k], &nn, &beta2,
                           U->ba[i][jpos], &mm );
                }
            }
        }

/*-------------------- reset double-pointer to -1 ( U-part) */
        for( j = 0; j < L->nzcount[i]; j++ )        {
	  col = L->ja[i][j];
	  jw[col] = -1;
        }
        jw[i] = -1;
        for( j = 0; j < U->nzcount[i]; j++ ) {
	  col = U->ja[i][j];
	  jw[col] = -1;
        }

/*-------------------- calculate truncated inverse of diagonal element of U */
	if (SVD)
	  ierr = invSVD(dim,lu->D[i]);
	else
	  ierr = invGauss(dim,lu->D[i]);
        if( ierr != 0 ) {
            for( j = i+1; j < n; j++ ) {
                lu->D[j] = NULL;
                L->ba[j] = NULL;
                U->ba[j] = NULL;
            }
            fprintf( fp, "fatal error: Singular diagonal block...\n" );
            return -2;
        }
    }
    lu->DiagOpt = 2;
    return 0;
}

int lofC( int lofM, vbsptr vbmat, vbiluptr lu, FILE *fp )
{
/*--------------------------------------------------------------------
 * symbolic ilu factorization to calculate structure of ilu matrix
 * for specified level of fill
 *--------------------------------------------------------------------
 * on entry:
 * =========
 * lofM     = level of fill, lofM >= 0
 * vbmat    = block matrix stored in VBSpaFmt format -- see globheads.h for
 *            details on format, size of blocks might be different
 * lu       = pointer to a VBILUSpar struct -- see globheads.h for details
 *            on format
 * fp       = file pointer for error log ( might be stderr )
 *--------------------------------------------------------------------
 * on return:
 * ==========
 * ierr     = return value.
 *            ierr  = 0   --> successful return.
 *            ierr != 0   --> error
 * lu->n    = dimension of the block matrix
 *   ->L    = L part -- stored in BSpaFmt format, patterns only in lofC
 *   ->U    = U part -- stored in BSpaFmt format, patterns only in lofC
 *------------------------------------------------------------------*/
    int n = vbmat->n;
    int *levls = NULL, *jbuf = NULL, *iw = lu->work;
    int **ulvl;  /*  stores lev-fils for U part of ILU factorization*/
    vbsptr L = lu->L, U = lu->U;
/*--------------------------------------------------------------------
 * n        = number of rows or columns in matrix
 * inc      = integer, count of nonzero(fillin) element of each row
 *            after symbolic factorization
 * ju       = entry of U part of each row
 * lvl      = buffer to store levels of each row
 * jbuf     = buffer to store column index of each row
 * iw       = work array
 *------------------------------------------------------------------*/
    int i, j, k, col, ip, it, jpiv;
    int incl, incu, jmin, kmin;

    levls  = (int *)Malloc( n*sizeof(int), "lofC" );
    jbuf = (int *)Malloc( n*sizeof(int), "lofC" );
    ulvl = (int **)Malloc( n*sizeof(int *), "lofC" );

    /* initilize iw */
    for( j = 0; j < n; j++ ) iw[j] = -1;
    for( i = 0; i < n; i++ ) {
        incl = 0;
        incu = i;
/*-------------------- assign lof = 0 for matrix elements */
        for( j = 0; j < vbmat->nzcount[i]; j++ ) {
            col = vbmat->ja[i][j];
            if( col < i ) {
/*-------------------- L-part  */
	        jbuf[incl] = col;
	        levls[incl] = 0;
	        iw[col] = incl++;
            }
            else if (col > i) {
/*-------------------- U-part  */
	        jbuf[incu] = col;
	        levls[incu] = 0;
	        iw[col] = incu++;
            }
        }
/*-------------------- symbolic k,i,j Gaussian elimination  */
        jpiv = -1;
        while (++jpiv < incl) {
            k = jbuf[jpiv] ;
/*-------------------- select leftmost pivot */
            kmin = k;
            jmin = jpiv;
            for( j = jpiv + 1; j< incl; j++) {
	        if( jbuf[j] < kmin ) {
	            kmin = jbuf[j];
	            jmin = j;
	        }
            }
/*-------------------- swap  */
            if( jmin != jpiv ) {
	        jbuf[jpiv] = kmin;
	        jbuf[jmin] = k;
	        iw[kmin] = jpiv;
	        iw[k] = jmin;
	        j = levls[jpiv] ;
	        levls[jpiv] = levls[jmin];
	        levls[jmin] = j;
	        k = kmin;
            }
/*-------------------- symbolic linear combinaiton of rows  */
            for( j = 0; j < U->nzcount[k]; j++ ) {
	        col = U->ja[k][j];
	        it = ulvl[k][j]+levls[jpiv]+1 ;
	        if( it > lofM ) continue;
	        ip = iw[col];
	        if( ip == -1 ) {
	            if( col < i) {
	                jbuf[incl] = col;
	                levls[incl] = it;
	                iw[col] = incl++;
                    }
	            else if( col > i ) {
	                jbuf[incu] = col;
	                levls[incu] = it;
	                iw[col] = incu++;
	            }
                }
                else
	            levls[ip] = min(levls[ip], it);
            }
        }   /* end - while loop */
/*-------------------- reset iw */
        for( j = 0; j < incl; j++ ) iw[jbuf[j]] = -1;
        for( j = i; j < incu; j++ ) iw[jbuf[j]] = -1;
/*-------------------- copy L-part */
        L->nzcount[i] = incl;
        if(incl > 0 ) {
            L->ja[i] = (int *)Malloc( incl*sizeof(int), "lofC" );
            memcpy( L->ja[i], jbuf, sizeof(int)*incl);
        }
/*-------------------- copy U - part        */
        k = incu-i;
        U->nzcount[i] = k;
        if( k > 0 ) {
            U->ja[i] = (int *)Malloc( sizeof(int)*k, "lofC" );
            memcpy(U->ja[i], jbuf+i, sizeof(int)*k );
/*-------------------- update matrix of levels */
            ulvl[i] = (int *)Malloc( k*sizeof(int), "lofC" );
            memcpy( ulvl[i], levls+i, k*sizeof(int) );
        }
    }

/*-------------------- free temp space and leave --*/
    free(levls);
    free(jbuf);
    for(i = 0; i < n-1; i++ ) {
        if (U->nzcount[i]) free(ulvl[i]) ;
    }
    free(ulvl);

    return 0;
}
