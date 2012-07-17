/*-------------------------------------------------------------------*
 * main test driver for ILUK                                         *
 *-------------------------------------------------------------------*
 * Na Li, Oct 31, 2001                                               *
 *                                                                   *
 * Report bugs / send comments to: saad@cs.umn.edu, nli@cs.umn.edu   *
 *-------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h> 
#include <math.h>
#include "globheads.h"
#include "defs.h" 
#include "protos.h"
#include "ios.h"  
#include <sys/time.h>

/*-------------------- protos */
void output_header( io_t *pio );
void output_result(int lfil, io_t *pio, int iparam );
int read_coo(double **AA, int **JA, int **IA, io_t *pio,
             double **hs, double **ol, int job);
int readhb_c(int *NN, double **AA, int **JA, int **IA, io_t *pio, 
             double **rhs, double **guess, int *rsa);
int read_inputs( char *in_file, io_t *pio );
int get_matrix_info( FILE *fmat, io_t *pio );
void randvec (double *v, int n);
/*-------------------- end protos */

/*------------------------------------------------------------*/
int main(){
  int ierr = 0;
/* --- plotting is for writing stats into OUT/... in raw form */
  int plotting = 0, skip_its = 0;
  char pltfile[256];
  FILE *fits = NULL;
/*-------------------- main structs and wraper structs.   */
  csptr csmat = NULL;  /* matrix in csr formt             */
  SMatptr MAT;         /* Matrix structure for matvecs    */
  SPreptr PRE;         /* general precond structure       */
  iluptr lu = NULL;    /* ilu preconditioner structure    */
  double *sol = NULL, *x = NULL, *rhs = NULL;
/*-------------------- temp Harwell Boeing arrays */
  double *AA;
  int *JA, *IA;
  int n, nnz, lfil, rsa;
/*-------------------- IO */
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
/*------------------ file "matfile" contains paths to matrices */
  if( NULL == ( fmat = fopen( "matfile", "r" ) ) ) {
    fprintf( flog, "Can't open matfile...\n" );
    exit(2);
  }
  memset( line, 0, MAX_LINE );
  fgets( line, MAX_LINE, fmat );
  if( ( numat = atoi( line ) ) <= 0 ) {
    fprintf( flog, "Invalid count of matrices...\n" );
    exit(3);
  }
/*-------------------- open file OUT/ILUK.out for all performance
                       results of this run (all matrices and params) 
                       also set io->PrecMeth */
  strcpy(io.outfile,"OUT/ILUK.out");
  strcpy(io.PrecMeth,"ILUK");
  if( NULL == ( io.fout = fopen( io.outfile, "w" ) ) ) {
    fprintf(flog,"Can't open output file %s...\n", io.outfile);
    exit(4);
  }
/* ------------------- LOOP THROUGH MATRICES -------------------*/
  for( mat = 1; mat <= numat; mat++ ) {
    if( get_matrix_info( fmat, &io ) != 0 ) {
      fprintf( flog, "Invalid format in matfile_hb...\n" );
      exit(5);
    }
    fprintf( flog, "MATRIX: %s...\n", io.MatNam );
/*------------------------- Read matrix */
    csmat = (csptr)Malloc( sizeof(SparMat), "main" );
/*-------------------- case: COO formats */
    if (io.Fmt > HB) {      
      ierr = read_coo(&AA,&JA, &IA, &io, &rhs, &sol,0);
      if (ierr == 0) 
        fprintf(flog,"matrix read successfully\n");
      else {
	fprintf(flog, "read_coo error = %d\n", ierr);
	exit(6);
      }
      n = io.ndim;
      nnz = io.nnz;
/*-------------------- conversion from COO to CSR format */
      if((ierr = COOcs(n, nnz, AA, JA, IA, csmat)) != 0) {
	fprintf(stderr, "mainARMS: COOcs error\n");
	return ierr;
      }
    }
    else if (io.Fmt == HB) {
/*-------------------- NOTE: (AA,JA,IA) is in CSR format */
      ierr = readhb_c(&n, &AA, &JA, &IA, &io, &rhs, &sol, &rsa);
      if(ierr != 0) {
	fprintf(flog, "readhb_c error = %d\n", ierr);
	exit(7);
      }
      nnz = io.nnz;
      if((ierr = CSRcs(n, AA, JA, IA, csmat, rsa)) != 0) {
	fprintf(flog, "readhb_c: CSRcs error\n");
	return ierr;
      }
    }
/*-------------------- free COO/HB arrays  */
    free( AA);    AA = NULL;
    free( JA );   JA = NULL;
    free( IA );   IA = NULL; 
/*---------------------------------------------------------*/
    x = (double *)Malloc( io.ndim * sizeof(double), "main" );
    output_header( &io );
/*-------------------- set initial lfil and tol */ 
    lfil = io.fill_lev;
    io.tol0 = 0.0;
/*--------------------------LOOP THROUGH PARAMETERS */
    for( iparam = 1; iparam <= io.nparam; iparam++ ) {
      fprintf( flog, "iparam = %d\n", iparam );
      lu = (iluptr)Malloc( sizeof(ILUSpar), "main" );
      fprintf( flog, "begin iluk(%d)\n",lfil );
      tm1 = sys_timer();
/*-------------------- call ILUK preconditioner set-up  */
      ierr = ilukC(lfil, csmat, lu, flog );
/*----------------------------------------------------- */
      tm2 = sys_timer();
      if( ierr == -2 ) {
	fprintf( io.fout, "zero diagonal element found...\n" );
	cleanILU( lu );
	goto NEXT_MAT;
      } else if( ierr != 0 ) {
	fprintf( flog, "*** iluk error, ierr != 0 ***\n" );
	exit(-1);
      }
      io.tm_p = tm2 - tm1;
      io.fillfact = (double)nnz_ilu( lu )/(double)(io.nnz + 1);
      fprintf( flog, "iluk ends, fill factor (mem used) = %f\n",\
      io.fillfact );
      if( skip_its ) {
	io.its = -1;
	io.tm_i = -1;
	io.enorm = -1;
	io.rnorm = -1;
	goto NEXT_PARA;
      }
/*------------------------ condition number estimation */
      if( condestLU( lu, flog ) != 0 ) {
	fprintf( flog, "Not attempting iterative solution.\n" );
	fprintf( io.fout, "Not attempting iterative solution.\n" );
	io.its = -1;
	io.tm_i = -1;
	io.enorm = -1;
	io.rnorm = -1;
	goto NEXT_PARA;
      }
/* ----------------------- initial guess */
      /*for( i = 0; i < io.ndim; i++ ) x[i] = 0.0;*/             
      randvec(x, n);
/*-------------------- create a file for printing
                       'its -- time -- res' info from fgmres */
      if (plotting ) { 
	sprintf( pltfile, "OUT/%s_ILUK_F%05d", io.MatNam, lfil);
	if( NULL == ( fits = fopen( pltfile, "w" ) ) ) {
	  fprintf( flog, "Can't open output file %s...\n", pltfile );
	  exit(8);
	}
      } else 
	fits = NULL;      
/*-------------------- set up the structs before calling fgmr */
      MAT->n = n;
      MAT->CS = csmat;
      MAT->matvec = matvecCSR; 
      PRE->ILU = lu;
      PRE->precon = preconILU;
/*-------------------- call fgmr */
      io.its = io.maxits;
      tm1 = sys_timer();
      fgmr(MAT, PRE, rhs, x, io.tol, io.im, &io.its, fits);
      tm2 = sys_timer();
      io.tm_i = tm2 - tm1;
      if( io.its < io.maxits ) 
        fprintf( flog, "param %03d OK: converged in %d steps...\n\n",\
        iparam, io.its );
      else 
        fprintf( flog, "not converged in %d steps...\n\n", io.maxits );
      if( fits )
        fclose( fits );
/*-------------------- calculate error norm */
      terr = 0.0;
      for( i = 0; i < io.ndim; i++ )
	terr += ( sol[i] - x[i] ) * ( sol[i] - x[i] );
      io.enorm = sqrt(terr);
/*-------------------- calculate residual norm */
      matvec( csmat, x, sol );
      terr = 0.0;
      for( i = 0; i < io.ndim; i++ )
	terr += ( rhs[i] - sol[i] ) * ( rhs[i] - sol[i] );
      io.rnorm = sqrt(terr);
/*-------------------- Test with next params   */           
NEXT_PARA:
      output_result( lfil, &io, iparam );
      lfil += io.fill_lev_inc;
      cleanILU(lu);
    }
/*-------------------- Test with next matrix   */              
NEXT_MAT:
    cleanCS( csmat );
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

