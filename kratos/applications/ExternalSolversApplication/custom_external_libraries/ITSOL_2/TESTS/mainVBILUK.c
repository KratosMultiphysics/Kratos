/*-----------------------------------------------------------------------*
 * main test driver for VBILUK                                           *
 *-----------------------------------------------------------------------*
 * Na Li, Aug 26, 2001  -- Y.Saad 07/05                                  *
 *                                                                       *
 * Report bugs / send comments to: saad@cs.umn.edu, nli@cs.umn.edu       *
 *-----------------------------------------------------------------------*/
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
void output_header_vb( io_t *pio );
void output_result(int lfil, io_t *pio, int iparam );
int read_coo(double **AA, int **JA, int **IA, io_t *pio,
             double **hs, double **ol,int job);
int readhb_c(int *NN, double **AA, int **JA, int **IA, io_t *pio, 
             double **rhs, double **guess, int *rsa);
int read_inputs( char *in_file, io_t *pio );
int get_matrix_info( FILE *fmat, io_t *pio );
void randvec (double *v, int n);
void output_blocks( int nBlock, int *nB, FILE *f );
void output_perm( int n, int *perm, FILE *f );
/*-------------------- end protos */

int main () { 
  int ierr = 0;
/*-------------------------------------------------------------------
 * options
 *-----------------------------------------------------------------*/
  int plotting = 0, output_mat = 0, diagscal = 0;
  char pltfile[256];
  FILE *fits = NULL;
  int lfil, nBlock, *nB = NULL, *perm = NULL;
/*-------------------- main structs and wraper structs.   */
  csptr csmat = NULL;  /* matrix in csr formt             */
  vbsptr vbmat = NULL; 
  vbiluptr lu = NULL;  /* vbilu preconditioner structure  */
  SMatptr MAT;         /* Matrix structure for matvecs    */
  SPreptr PRE;         /* general precond structure       */
  double *sol = NULL, *x = NULL, *prhs = NULL, *rhs = NULL;
/*-------------------- temp COO/HB arrays - fortran style */
  int *IA,  *JA;
  double *AA;
  int rsa, n, nnz; 
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
/*------------------ set any parameters manually */
/* io.eps is the angle tolerance for grouping two columns in same
   supernode. This is a cosine and should be <= 1.  */
  io.eps = 0.8;
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
/*-------------------- open file OUT/VBILUK.out for all performance
                       results of this run (all matrices and params) 
                       also set io->PrecMeth */
  strcpy(io.outfile,"OUT/VBILUK.out");
  strcpy(io.PrecMeth,"Variable Block ILUK (VBILUK)");
  if( NULL == ( io.fout = fopen( io.outfile, "w" ) ) ) {
    fprintf(flog,"Can't open output file %s...\n", io.outfile);
    exit(4);
  }
/*-------------------- LOOP through MATRICES */
  for( mat = 1; mat <= numat; mat++ ) {
    if( get_matrix_info( fmat, &io ) != 0 ) {
      fprintf( flog, "Invalid format in matfile_hb...\n" );
      exit(5);
    }
    fprintf( flog, "MATRIX: %s...\n", io.MatNam );
/* ------------------- Read in matrix and allocate memory */
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
/*-------------- Free memory */    
    free( AA);    AA = NULL;
    free( JA );   JA = NULL;
    free( IA );   IA = NULL; 
/*-------------------- diag scaling */
    if (diagscal == 1) {  
      int nrm=1;
      double *diag;
      diag = (double *)Malloc( sizeof(double)*n, "mainILUC:diag" );
      ierr = roscalC( csmat, diag, nrm);
      if( ierr != 0 ) {
	fprintf( stderr, "main-vbiluk: roscal: a zero row...\n" );
	return ierr;
      }
      ierr = coscalC( csmat, diag, nrm );
      if( ierr != 0 ) {
	fprintf( stderr, "main-vbiluk: roscal: a zero col...\n" );
	return ierr;
      }
      free(diag);
    }
/*------------------------------------------------------------*/
    x = (double *)Malloc( io.ndim * sizeof(double), "main" );
    ierr = init_blocks( csmat, &nBlock, &nB, &perm, 
			io.eps, &io.tm_h, &io.tm_a );
    io.tm_b = io.tm_h + io.tm_a;
    if( ierr != 0 ) {
      fprintf( flog, "*** in init_blocks ierr != 0 ***\n" );
      exit(8);
    }
/*--------------- permutes the rows and columns of the matrix */
    if( dpermC( csmat, perm ) != 0 ) {
      fprintf( flog, "*** dpermC error ***\n" );
      exit(9);
    }
/*-------------------- permute the right-hand-side  */
    prhs = (double *)Malloc(n*sizeof(double), "main" );
    for( i = 0; i < n; i++ ) 
      prhs[perm[i]] = rhs[i];
/*-------------------- convert to block matrix. */
    vbmat = (vbsptr)Malloc( sizeof(VBSparMat), "main" );
    ierr = csrvbsrC( 1, nBlock, nB, csmat, vbmat );
    if( ierr != 0 ) {
      fprintf( flog, "*** in csrvbsr ierr != 0 ***\n" );
      exit(10);
    }
/*------------------ OUTPUT MATRIX */
    if( output_mat ) {
      char matdata[MAX_LINE];
      FILE *fmatlab;
      int ii, jj;
      sprintf( matdata, "OUT/%s.dat", io.MatNam );
      if( NULL != ( fmatlab = fopen( matdata, "w" ) ) ) {
	fprintf( fmatlab, "%d %d 0\n", csmat->n, csmat->n );
	for( ii = 0; ii < csmat->n; ii++ )
	  for( jj = 0; jj < csmat->nzcount[ii]; jj++ )
	    fprintf( fmatlab, "%d %d 1\n", ii+1, csmat->ja[ii][jj]+1 );
	fclose( fmatlab );
      }
    }
/*------------- info of Block matrix */
    io.rt_v = (double)csmat->n / (double)vbmat->n;
    io.rt_e = (double)nnz_cs( csmat ) / (double)nnzVBMat( vbmat );
    io.ceff = (double)nnz_cs( csmat ) / (double)memVBMat( vbmat ) * 100;
/*---------------------------*/    
    output_header_vb( &io );
    lfil = io.fill_lev;
    io.tol0 = 0.0; /* make sure this is set to zero*/
/*---------------------- LOOP through parameters */
    for( iparam = 1; iparam <= io.nparam; iparam++ ) {
      fprintf( flog, "iparam = %d\n", iparam );
      lu = (vbiluptr)Malloc( sizeof(VBILUSpar), "main" );
      fprintf( flog, "begin vbiluk\n" );
      tm1 = sys_timer();
/*-------------------- call VBILUK preconditioner set-up  */
      ierr = vbilukC( lfil, vbmat, lu, flog );
/*----------------------------------------------------- */
      tm2 = sys_timer();
      if( ierr == -2 ) {
	fprintf( io.fout, "Singular diagonal block...\n" );
	cleanVBILU( lu );
	goto NEXT_MAT;
      } else if( ierr != 0 ) {
	fprintf( flog, "*** vbilu error, ierr != 0 ***\n" );
	exit(11);
      }
      io.tm_p = tm2 - tm1;
      io.fillfact = (double)nnz_vbilu( lu )/(double)(io.nnz + 1);
      fprintf( flog, "vbiluk ends, fill factor (mem used) = %f\n",\
      io.fillfact );
/*------------- get rough idea of cond number - exit if too big */     
      if( VBcondestC( lu, flog ) != 0 ) {
	fprintf( flog, "Not attempting iterative solution.\n" );
	fprintf( io.fout, "Not attempting iterative solution.\n" );
	io.its = -1;
	io.tm_i = -1;
	io.enorm = -1;
	io.rnorm = -1;
	goto NEXT_PARA;
      }
/*-------------------- initial guess */      
/* for( i = 0; i < io.ndim; i++ ) x[i] = 0.0; */
      randvec(x, n);
/*-------------------- create a file for printing
                       'its -- time -- res' info from fgmres */
      if (plotting ) {
	sprintf( pltfile, "OUT/%s_VBILUK_F%05d", io.MatNam, lfil);
	if( NULL == ( fits = fopen( pltfile, "w" ) ) ) {
	  fprintf( flog, "Can't open output file %s...\n", pltfile );
	  exit(12);
	}
      } else 
	fits = NULL;
/*-------------------- set up the structs before calling fgmr */
      MAT->n = n;
      MAT->CS = csmat;
      MAT->matvec = matvecCSR; 
      PRE->VBILU = lu; 
      PRE->precon = preconVBR;
/*-------------------- call fgmr */
      io.its = io.maxits;
      tm1 = sys_timer();
      fgmr(MAT, PRE, prhs, x, io.tol, io.im, &io.its, fits);
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
/*-------------------- calculate the error norm */
/*-------------------- P*sol ?= x */
      terr = 0.0;
      for( i = 0; i < io.ndim; i++ )
	terr += (sol[i] - x[perm[i]]) * (sol[i] - x[perm[i]]);
      io.enorm = sqrt(terr);
/*-------------------- calculate the residual norm */
      vbmatvec( vbmat, x, sol );
      terr = 0.0;
      for( i = 0; i < io.ndim; i++ )
	terr += (prhs[i] - sol[i]) * (prhs[i] - sol[i]);
      io.rnorm = sqrt(terr);
/*-------------------- next params */      
    NEXT_PARA:
      output_result(lfil, &io, iparam );
      lfil += io.fill_lev_inc;
      cleanVBILU( lu );
    }
/*-------------------- next matrix */          
  NEXT_MAT:
    /*  output_blocks( nBlock, nB, io.fout );*/
    cleanCS( csmat );
    cleanVBMat( vbmat );
    free( nB );
    free( perm );
    free( sol );
    free( x );
    free( prhs );
    free( rhs );
  }

  fclose( io.fout );  
  if( flog != stdout ) fclose ( flog );
  fclose( fmat );
  free(MAT);
  free(PRE);
  return 0;
}

