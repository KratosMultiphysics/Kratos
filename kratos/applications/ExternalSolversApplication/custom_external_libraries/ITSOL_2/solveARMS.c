/*-----------------------------------------------------------------*
 * main test driver for the ARMS2 preconditioner for
 * Matrices in the COO/Harwell Boeing format
 *-----------------------------------------------------------------*
 * Yousef Saad - Aug. 2005.                                        *
 *                                                                 *
 * Report bugs / send comments to: saad@cs.umn.edu                 *
 *-----------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "globheads.h"
#include "defs.h"
#include "protos.h"
#include "ios.h"
#include <time.h>
#include <assert.h>


#define TOL_DD 0.7   /* diagonal dominance tolerance for */
/* independent sets                 */

void output_header(io_t *pio);
void output_result(int lfil, io_t *pio, int iparam);
void set_arms_pars(io_t* io,  int Dscale, int *ipar,
                   double *tolcoef, int *lfil);
int read_inputs(char *in_file, io_t *pio);
int get_matrix_info(FILE *fmat, io_t *pio);
int read_coo(double **AA, int **JA, int **IA, io_t *pio,
             double **hs, double **ol, int job);
int readhb_c(int *NN, double **AA, int **JA, int **IA, io_t *pio,
             double **rhs, double **guess, int *rsa);
void randvec (double *v, int n);
int dumpArmsMat(arms PreSt, FILE *ft);

int solveARMS(io_t* io , int echo_level, double* AA, int * IA, int* JA, double *x, double* rhs)
{
    int n = io->ndim; 
    
    int nnz = io->nnz;
/*----------------------------------------------------------
     * solves a system using ARMS preconditioned fgmres
     *----------------------------------------------------------*/
    FILE *flog = NULL; //stdout;
    if(echo_level >= 1) flog = stdout;
    int ierr = 0;
    /*-------------------- options*/
    int diagscal = 1;
    //char pltfile[256];
     FILE *fits = NULL; //stdout;
    if(echo_level >= 2) fits = stdout;
    
    double tol, tolind = TOL_DD;
    int j,  lfil;
    /*-------------------- main structs and wrapper structs.     */
    csptr csmat = NULL;    /* matrix in csr formt             */
    arms ArmsSt = NULL;    /* arms preconditioner structure   */
    SMatptr MAT = NULL;    /* Matrix structure for matvecs    */
    SPreptr PRE = NULL;    /* general precond structure       */
    //double *sol = NULL;
    /*---------------- method for incrementing lfil is set here */
    int lfil_arr[7];
    double droptol[7], dropcoef[7];
    int ipar[18];
    
    
    /*-------------------- harwell boeing temporary arrays */
    //int rsa = 0;

    double tm1, tm2;
    int  iparam, i;
    //double terr;
    //char line[MAX_LINE];
    MAT = (SMatptr)Malloc(sizeof(SMat), "main:MAT");
    PRE = (SPreptr)Malloc(sizeof(SPre), "main:PRE");
    
    /*------------------ read and set parameters and other inputs */
//     memset(&io, 0, sizeof(io));
// // HERE SET THE INPUTS HARD CODED
//     io.nparam = 1; //     1. nparam  = number of tests for the preconditioner (see below)
//     io.ndim = n;  //     
//     io.nnz = nnz;
//     io.im = 200;  //     2. dim     = dimension of Krylov subspace in (outer) FGMRES
//     io.maxits = 200; //             3. maxits  = maxits in outer fgmres.
//     io.tol = 1.0e-8; //         4. tol     = tolerance for stopping iteration
//     io.lfil0 = 50;  //            5. lfil0   = initial lfil
//     io.lfilInc = 1;  //            6. lfilInc = increment for lfil
//     io.tol0 = 0.001; //          7. tol0    = initial tol
//     io.tolMul = 0.01;  //          8. tolMul  = multiple increment for tol0
//     io.fill_lev = 1;  //            9. USED BY ILUK ONLY: fill_lev = fill level
//     io.perm_type = 2;  //           10. ARMS ONLY: PQ perms or Ind. Sets.
//     io.Bsize  = 30; //            11. ARMS ONLY: Block-size for independent sets/last block

    /*------------------ file "matfile" contains paths to matrices*/

    /*-------------------- set parameters for arms */
    set_arms_pars(io, diagscal, ipar, dropcoef, lfil_arr);
    
    ipar[0] = 5; //setting the maximum number of levels.
    
    if(echo_level < 1) ipar[3] = 0;



    //do the reading
     csmat = (csptr)Malloc( sizeof(SparMat), "main:csmat" );

     //copy to ITSOL matrix structure
    int  j1, len;
  double *bra;
  int *bja;
  /*    setup data structure for mat (csptr) struct */
  setupCS( csmat, n, 1 );

  for (j=0; j<n; j++) {
    len = IA[j+1] - IA[j];
    csmat->nzcount[j] = len;
    if (len > 0) {
      bja = (int *) Malloc( len*sizeof(int), "CSRcs" );
      bra = (double *) Malloc( len*sizeof(double), "CSRcs" );
      i = 0;
      for (j1=IA[j]; j1<IA[j+1]; j1++) {
        bja[i] = JA[j1] ;
        bra[i] = AA[j1] ;
        i++;
      }
      csmat->ja[j] = bja;
      csmat->ma[j] = bra;
    }
  }    

    /*----------------------------------------------------------*/
    n = csmat->n;
    /*-------------------- set initial lfil and tol */
    lfil = io->lfil0;
    tol  = io->tol0;

    /*-------------------- LOOP THROUGH PARAMETERS */
    for (iparam = 1; iparam <= io->nparam; iparam++)
    {
        if(echo_level >= 1) fprintf(flog, "Parameter case = %d\n", iparam);
        for (j=0; j<7; j++)
        {
            lfil_arr[j] = lfil*((int) nnz/n);
            droptol[j] = tol*dropcoef[j];
        }

        ArmsSt = (arms) Malloc(sizeof(armsMat),"main:ArmsSt");
        setup_arms(ArmsSt);
	if(echo_level >= 1)   fprintf(flog, "begin arms\n");
        tm1 = sys_timer();
        /*-------------------- call ARMS preconditioner set-up  */
        ierr =
            arms2(csmat, ipar, droptol, lfil_arr, tolind, ArmsSt, flog);
        /*----------------------------------------------------- */
        tm2 = sys_timer();
        if (ierr != 0)
        {
            //fprintf(io->fout, " ** ARMS2 error - code %d...\n",ierr);
            io->its = -1;
            io->tm_i = -1;
            io->enorm = -1;
            io->rnorm = -1;
	    return 1;
//             goto NEXT_PARA;
        }
        io->tm_p = tm2 - tm1;
        io->fillfact = (double)nnz_arms(ArmsSt, flog)/(double)(nnz + 1);
	
	if(echo_level >= 1) fprintf(flog, "ARMS ends, fill factor (mem used) = %f\n",  io->fillfact);
	
        /*---------------- get rough idea of cond number - exit if too big */
        if (condestArms(ArmsSt, x, flog) != 0)
        {
	    if(echo_level >= 1)  fprintf(flog, "Not attempting iterative solution.\n");
            //fprintf(io->fout, "Not attempting iterative solution.\n");
            io->its = -1;
            io->tm_i = -1;
            io->enorm = -1;
            io->rnorm = -1;
	    return 1;
//             goto NEXT_PARA;
        }

/*-------------------- set up the structs before calling fgmr */
        MAT->n = n;
        MAT->CS = csmat;
        MAT->matvec = matvecCSR;
        PRE->ARMS = ArmsSt;
        PRE->precon = preconARMS;
        /*-------------------- call fgmr */

	io->its = io->maxits;
        tm1 = sys_timer();
        fgmr(MAT, PRE, rhs, x, io->tol, io->im, &io->its, fits);
        tm2 = sys_timer();
        io->tm_i = tm2 - tm1;
        if (io->its < io->maxits)
	  {
	  if(echo_level >= 1) fprintf(flog, "param %03d OK: converged in %d steps...\n", iparam, io->its);
	  }
        else
	  {
	  if(echo_level >= 1)  fprintf(flog, "not converged in %d steps...\n", io->maxits);
	  }
//         if (fits)
//             fclose(fits);
        /*-------------------- actual error norm */
        /*      terr = 0.0;
              for( i = 0; i < n; i++)
                terr += ( x[i] - sol[i]) * ( x[i] - sol[i]);
        	io->enorm = sqrt(terr);*/
        /*-------------------- calculate residual norm from generated rhs */
        /*      matvec(csmat, x, sol);
              terr = 0.0;
              for(i = 0; i < io->ndim; i++)
        	terr += (rhs[i] - sol[i]) * (rhs[i] - sol[i]);
              io->rnorm = sqrt(terr);*/
        /*-------------------- go to next param case */
// NEXT_PARA:
//         output_result(lfil, &io, iparam);
//         lfil += io->lfilInc;
//         tol  *= io->tolMul;
        cleanARMS(ArmsSt);

    }
    /*-------------------- NEXT_MATRIX: */
    cleanCS(csmat);
    //free(sol);

    free(MAT);
    free(PRE);
    
    return 0;
}
