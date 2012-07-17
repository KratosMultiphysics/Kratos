#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "globheads.h"
#include "protos.h"

#ifndef min
#define min(a,b) (((a)>(b))?(b):(a))
#endif
#ifndef max
#define max(a,b) (((a)>(b))?(a):(b))
#endif

/*--------------------protos */
int lofC( int lofM, csptr csmat, iluptr lu, FILE *fp ); 
/*--------------------end protos */

int ilukC( int lofM, csptr csmat, iluptr lu, FILE *fp )
{
/*----------------------------------------------------------------------------
 * ILUK preconditioner
 * incomplete LU factorization with level of fill dropping
 *----------------------------------------------------------------------------
 * Parameters
 *----------------------------------------------------------------------------
 * on entry:
 * =========
 * lofM     = level of fill: all entries with level of fill > lofM are
 *            dropped. Setting lofM = 0 gives BILU(0).
 * csmat    = matrix stored in SpaFmt format -- see heads.h for details
 *            on format
 * lu       = pointer to a ILUKSpar struct -- see heads.h for details
 *            on format
 * fp       = file pointer for error log ( might be stderr )
 *
 * on return:
 * ==========
 * ierr     = return value.
 *            ierr  = 0   --> successful return.
 *            ierr  = -1  --> error in lofC
 *            ierr  = -2  --> zero diagonal found
 * lu->n    = dimension of the matrix
 *   ->L    = L part -- stored in SpaFmt format
 *   ->D    = Diagonals
 *   ->U    = U part -- stored in SpaFmt format
 *----------------------------------------------------------------------------
 * Notes:
 * ======
 * All the diagonals of the input matrix must not be zero
 *--------------------------------------------------------------------------*/
    int ierr;
    int n = csmat->n;
    int *jw, i, j, k, col, jpos, jrow;
    csptr L, U;
    double *D;

    setupILU( lu, n );
    L = lu->L;
    U = lu->U;
    D = lu->D;

    /* symbolic factorization to calculate level of fill index arrays */
    if( ( ierr = lofC( lofM, csmat, lu, fp ) ) != 0 ) {
      fprintf( fp, "Error: lofC\n" );
      return -1;
    }

    jw = lu->work;
    /* set indicator array jw to -1 */
    for( j = 0; j < n; j++ ) jw[j] = -1;

    /* beginning of main loop */
    for( i = 0; i < n; i++ ) {
        /* set up the i-th row accroding to the nonzero information from
           symbolic factorization */
        mallocRow( lu, i );

        /* setup array jw[], and initial i-th row */
        for( j = 0; j < L->nzcount[i]; j++ ) {  /* initialize L part   */
            col = L->ja[i][j];
            jw[col] = j;
            L->ma[i][j] = 0;
        }
        jw[i] = i;
        D[i] = 0; /* initialize diagonal */
        for( j = 0; j < U->nzcount[i]; j++ ) {  /* initialize U part   */
            col = U->ja[i][j];
            jw[col] = j;
            U->ma[i][j] = 0;
        }

        /* copy row from csmat into lu */
        for( j = 0; j < csmat->nzcount[i]; j++ ) {
            col = csmat->ja[i][j];
            jpos = jw[col];
            if( col < i )
                L->ma[i][jpos] = csmat->ma[i][j];
            else if( col == i )
                D[i] = csmat->ma[i][j];
            else
                U->ma[i][jpos] = csmat->ma[i][j];
        }

        /* eliminate previous rows */
        for( j = 0; j < L->nzcount[i]; j++ ) {
            jrow = L->ja[i][j];
            /* get the multiplier for row to be eliminated (jrow) */
            L->ma[i][j] *= D[jrow];

            /* combine current row and row jrow */
            for( k = 0; k < U->nzcount[jrow]; k++ ) {
                col = U->ja[jrow][k];
                jpos = jw[col];
                if( jpos == -1 ) continue;
                if( col < i )
                    L->ma[i][jpos] -= L->ma[i][j] * U->ma[jrow][k];
                else if( col == i )
                    D[i] -= L->ma[i][j] * U->ma[jrow][k];
                else
                    U->ma[i][jpos] -= L->ma[i][j] * U->ma[jrow][k];
            }
        }

        /* reset double-pointer to -1 ( U-part) */
        for( j = 0; j < L->nzcount[i]; j++ )
        {
            col = L->ja[i][j];
            jw[col] = -1;
        }
        jw[i] = -1;
        for( j = 0; j < U->nzcount[i]; j++ )
        {
            col = U->ja[i][j];
            jw[col] = -1;
        }

        if( D[i] == 0 ) {
            for( j = i+1; j < n; j++ ) {
                L->ma[j] = NULL;
                U->ma[j] = NULL;
            }
            fprintf( fp, "fatal error: Zero diagonal found...\n" );
            return -2;
        }
        D[i] = 1.0 / D[i];
    }

    return 0;
}

int lofC( int lofM, csptr csmat, iluptr lu, FILE *fp )
{
/*--------------------------------------------------------------------
 * symbolic ilu factorization to calculate structure of ilu matrix
 * for specified level of fill
 *--------------------------------------------------------------------
 * on entry:
 * =========
 * lofM     = level of fill, lofM >= 0
 * csmat    = matrix stored in SpaFmt format -- see heads.h for details
 *            on format
 * lu       = pointer to a ILUSpar struct -- see heads.h for details
 *            on format
 * fp       = file pointer for error log ( might be stderr )
 *--------------------------------------------------------------------
 * on return:
 * ==========
 * ierr     = return value.
 *            ierr  = 0   --> successful return.
 *            ierr != 0   --> error
 * lu->n    = dimension of the block matrix
 *   ->L    = L part -- stored in SpaFmt format, patterns only in lofC
 *   ->U    = U part -- stored in SpaFmt format, patterns only in lofC
 *------------------------------------------------------------------*/
    int n = csmat->n;
    int *levls = NULL, *jbuf = NULL, *iw = lu->work;
    int **ulvl;  /*  stores lev-fils for U part of ILU factorization*/
    csptr L = lu->L, U = lu->U;
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
        for( j = 0; j < csmat->nzcount[i]; j++ ) {
            col = csmat->ja[i][j];
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

