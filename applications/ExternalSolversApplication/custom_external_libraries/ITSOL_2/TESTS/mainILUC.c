/*-------------------------------------------------------------------*
 * main test driver for Crout version of ILU
 *-------------------------------------------------------------------*
 * Na Li, Mar 25, 2002                                               *
 *                                                                   *
 * Report bugs / send comments to: saad@cs.umn.edu, nli@cs.umn.edu   *
 *-------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "globheads.h"
#include "protos.h"  
#include "ios.h" 
#include <time.h>

/*-------------------- protos */
void output_header(io_t *pio );
void output_result( int lfil, io_t *pio, int iparam );
int read_coo(double **AA, int **JA, int **IA, io_t *pio,
	     double **hs, double **ol, int job); 
int readhb_2(int*, double**, int**, int**, 
             io_t*, double**, double**, int*, int);
int read_inputs( char *in_file, io_t *pio );
int get_matrix_info( FILE *fmat, io_t *pio );
void randvec (double *v, int n);
int CondestC(iluptr lu, FILE *fp );
void matvecLDU(SMatptr mat, double *x, double *y); 
void matvecCSC(SMatptr mat, double *x, double *y); 
int preconLDU(double *x, double *y, SPreptr mat);
void coocsc(int,int,double*,int*,int*,double**,int**,int**,int);
/*-------------------- end protos */

#define DRP_MTH 0          /* drop method see ilutc code */
int main() { 
  int ierr = 0;
/*-------------------------------------------------------------------
 * options
 *-----------------------------------------------------------------*/
  int dropmthd = DRP_MTH, plotting = 0, nnz, output_lu = 0; 
  int pattern_symm = 0;
  char pltfile[256];
  FILE *fits = NULL;
/*-------------------- main structs and wraper structs.     */
  csptr csmat = NULL;    /* matrix in csr formt             */
  SMatptr MAT = NULL;    /* Matrix structure for matvecs    */
  SPreptr PRE = NULL;    /* General precond structure       */
  iluptr lumat = NULL;   /* ilu preconditioner structure    */
  iluptr lu = NULL;      /* a temporary lu matrix           */
  double *sol = NULL, *x = NULL, *rhs = NULL;
/*-------------------- method for incrementing lfil is set here */
  int lfil; 
  double tol;
/*-------------------- harwell boeing temp. arrays */
  double *AA;
  int *IA, *JA;
  int rsa=0; 
  int n; 
/*-------------------- IO-related */  
  FILE *flog = stdout, *fmat = NULL;
  io_t io;
  double tm1, tm2;
  int mat, numat, iparam, i;
  double terr;
  char line[MAX_LINE];
  MAT = (SMatptr)Malloc( sizeof(SMat), "main:MAT" );
  PRE = (SPreptr)Malloc( sizeof(SPre), "main:PRE" );
/*------------------ read and set parameters and other inputs  */
  memset( &io, 0, sizeof(io) );
  if( read_inputs( "inputs", &io ) != 0 ) {
    fprintf( flog, "Invalid inputs file...\n" );
    exit(1); 
  }
/*-----------------------------------------------------------------*/
  if( NULL == ( fmat = fopen( "matfile", "r" ) ) ) {
    fprintf( flog, "Can't open matfile...\n" );
    exit(2) ; 
  }
  memset( line, 0, MAX_LINE );
  fgets( line, MAX_LINE, fmat );
  if( ( numat = atoi( line ) ) <= 0 ) {
    fprintf( flog, "Invalid count of matrices...\n" );
    exit(3); 
  }
/*-------------------- open file OUT/ILUC.out for all performance
                       results of this run (all matrices and params) 
                       also set io->PrecMeth */
  strcpy(io.outfile,"OUT/ILUC.out");
  strcpy(io.PrecMeth,"ILUC");
  if( NULL == ( io.fout = fopen(io.outfile,"w" ) ) ) {
    fprintf(flog,"Can't open output file %s...\n", io.outfile);
    exit(4);
  }  
/*-------------------- LOOP THROUGH MATRICES */
  for( mat = 1; mat <= numat; mat++ ) {
    if( get_matrix_info( fmat, &io ) != 0 ) {
      fprintf( flog, "Invalid format in matfile...\n" );
      exit(5); 
    }
    fprintf( flog, "MATRIX: %s...\n", io.MatNam);
/*-------------------- Read  matrix */
    lumat = (iluptr)Malloc(sizeof(LDUmat),"main:lumat");
    csmat = (csptr)Malloc( sizeof(SparMat), "main:csmat");    
/*-------------------- case: COO formats (0-indexing) */
    if (io.Fmt > HB) {      
      ierr = read_coo(&AA, &JA, &IA, &io, &rhs, &sol,0);
      if (ierr == 0) 
        fprintf(flog,"matrix read successfully\n");
      else {
	fprintf(flog, "read_coo error = %d\n", ierr);
	exit(6);
      }
      n = io.ndim;
      nnz = io.nnz;
/*-------------------- conversion from COO to CS format *
                       NOTE: csmat is in colum format   */
      if((ierr = COOcs(n, nnz, AA, IA, JA, csmat)) != 0) {
	fprintf(stderr, "mainILUC: COOcs error\n");
	return ierr;
      }
    }
    else if (io.Fmt == HB) {
/*-------------------- NOTE : (AA,JA,IA) is in CSC format */
      ierr = readhb_2(&n, &AA, &JA, &IA, &io, &rhs, &sol, &rsa, 0);
      if(ierr != 0) {
	fprintf(flog, "readhb_c error = %d\n", ierr);
	exit(7);
      }
      nnz = io.nnz;
/*-------------------- convert matrix to cs format for matvecs */
      if((ierr = CSRcs(n, AA, JA, IA, csmat, rsa)) != 0) {
	fprintf(flog, "CSRcs error\n");
	return (ierr);
      }
    }
/*-------------------- convert to lum format for iluc + symmetrize */
    if( rsa == 0 && pattern_symm)
      rsa = 2;  
    if( ( ierr = CSClumC(csmat, lumat, rsa ) ) != 0 ) {
     fprintf( stderr, " error: CSClum error\n" );
     return (ierr);
    }
/*-------------------- free CSR/CSC arrays no longer needed  */
    free(AA); AA = NULL;
    free(IA); IA = NULL;
    free(JA); JA = NULL;
/*-------- ##FIXME: diag scaling removed - see old version */
/*---------------------------------------------------------*/
    x = (double *)Malloc( io.ndim * sizeof(double), "main" );
    output_header( &io );
/*-------------------- set initial lfil and tol */ 
   lfil = io.lfil0;
   tol  = io.tol0;
/*-------------------- LOOP through parameters */
    for( iparam = 1; iparam <= io.nparam; iparam++ ) {
      fprintf( flog, "iparam = %d\n", iparam );
      lu = (iluptr)Malloc( sizeof(ILUSpar), "main" );
      tm1 = sys_timer();
/*-------------------- call ILUC preconditioner set-up  */
      ierr = ilutc( lumat, lu, lfil, tol, dropmthd, flog);
      tm2 = sys_timer();
      if( ierr != 0 ) {
	fprintf( io.fout, " *** ILUC error - code %d \n", ierr);
	io.its = -1;
	io.tm_i = -1;
	io.enorm = -1;
	io.rnorm = -1;
	goto NEXT_PARA;
      }
      if( output_lu ){
	char matdata[MAX_LINE];
	sprintf( matdata, "OUT/%s.dat",io.MatNam );
	outputLU( lu, matdata );
      }
      io.tm_p = tm2 - tm1;
      io.fillfact = (double)nnz_ilu( lu )/(double)(io.nnz + 1);
      fprintf( flog, "ilutc ends, fill factor (mem used) = %f\n",\
	       io.fillfact); 
      
      if( CondestC( lu, flog ) != 0 ) {
	fprintf( flog, "Not attempting iterative solution.\n" );
	fprintf( io.fout, "Not attempting iterative solution.\n" );
	io.its = -1;
	io.tm_i = -1;
	io.enorm = -1;
	io.rnorm = -1;
	goto NEXT_PARA;
      }      
/*-------------------- initial guess */     
/* */  for( i = 0; i < n; i++ )  x[i] = 0.0; 
//      randvec(x, n);          
/*-------------------- create a file for printing
                       'its -- time -- res' info from fgmres */
      if(plotting ) { 
	sprintf( pltfile, "OUT/%s_ILUC_F%05d_T%08.6f",\
        io.MatNam, lfil,tol);
	if( NULL == ( fits = fopen( pltfile, "w" ) ) ) {
	  fprintf( flog, "Can't open output file %s...\n",\
          pltfile );
	  exit(8);
	}
      } else 
	fits = NULL;
/*-------------------- set up the structs before calling fgmr */
      MAT->n = n;
      MAT->CS = csmat; /* in column format */
      MAT->matvec = matvecCSC; /* column matvec */
      PRE->ILU = lu;
      PRE->precon = preconLDU;
/*-------------------- call fgmr */
      io.its = io.maxits;
      tm1 = sys_timer();
      fgmr(MAT, PRE, rhs, x, io.tol, io.im, &io.its, fits );
      tm2 = sys_timer();
      io.tm_i = tm2 - tm1;
      if( io.its < io.maxits ) 
        fprintf( flog, "param %03d OK: converged in %d steps...\n\n",\
        iparam, io.its );
      else 
        fprintf( flog, "not converged in %d steps...\n\n",\
        io.maxits );
      if( fits ) 
        fclose( fits );
/*-------------------- calculate error norm */
      terr = 0.0;
      for( i = 0; i < io.ndim; i++ )
	terr += ( x[i] - sol[i] ) * ( x[i] - sol[i] );
      io.enorm = sqrt(terr);
/*-------------------- calculate res norm */
      matvecCSC( MAT, x, sol );
      terr = 0.0;
      for( i = 0; i < io.ndim; i++ )
	terr += ( rhs[i] - sol[i] ) * ( rhs[i] - sol[i] );
      io.rnorm = sqrt(terr);
/*-------------------- Test with next params   */     
NEXT_PARA:
      output_result( lfil, &io, iparam );
      lfil += io.lfilInc;
      tol  *= io.tolMul;
      cleanILU(lu);
    }
/*-------------------- NEXT_MAT: */
    cleanCS(csmat);
    cleanILU(lumat); 
    free( sol );
    free( x );
    free( rhs );
  }
  fclose( io.fout );
  if( flog != stdout ) fclose ( flog );
  fclose( fmat );
  free(MAT) ; 
  free (PRE); 
  return 0;
}
