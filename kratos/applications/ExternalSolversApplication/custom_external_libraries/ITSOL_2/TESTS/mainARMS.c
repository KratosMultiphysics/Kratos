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

/*-------------------- protos */
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
/*-------------------- end protos */

int main ()
{
    /*----------------------------------------------------------
     * solves a system using ARMS preconditioned fgmres
     *----------------------------------------------------------*/
    int ierr = 0;
    /*-------------------- options*/
    int plotting=0,  diagscal = 1;
    char pltfile[256];
    FILE *fits = NULL;
    double tol, tolind = TOL_DD;
    int j, nnz=0, lfil;
    /*-------------------- main structs and wraper structs.     */
    csptr csmat = NULL;    /* matrix in csr formt             */
    arms ArmsSt = NULL;    /* arms preconditioner structure   */
    SMatptr MAT = NULL;    /* Matrix structure for matvecs    */
    SPreptr PRE = NULL;    /* general precond structure       */
    double *sol = NULL, *x = NULL, *rhs = NULL;
    /*---------------- method for incrementing lfil is set here */
    int lfil_arr[7];
    double droptol[7], dropcoef[7];
    int ipar[18];
    /*-------------------- harwell boeing temporary arrays */
    double *AA;
    int *IA, *JA;
    int rsa;
    int n;
    /*-------------------- IO-related           */
    FILE *flog = stdout;   /* to output stats */
    FILE *fmat = NULL;     /* matrix file     */
    io_t io;               /* structure for handling I/O
                            functions + a few other things */
    double tm1, tm2;
    int mat, numat, iparam, i;
    double terr;
    char line[MAX_LINE];
    MAT = (SMatptr)Malloc(sizeof(SMat), "main:MAT");
    PRE = (SPreptr)Malloc(sizeof(SPre), "main:PRE");
    /*------------------ read and set parameters and other inputs */
    memset(&io, 0, sizeof(io));
    if (read_inputs("inputs", &io) != 0)
    {
        fprintf(flog, "ERROR reading inputs from file...\n");
        exit(1);
    }
    /*------------------ file "matfile" contains paths to matrices*/
    if (NULL == (fmat = fopen("matfile", "r")))
    {
        fprintf(flog, "Can't open matfile...\n");
        exit(2);
    }
    memset(line, 0, MAX_LINE);
    fscanf(fmat,"%d",&numat);
    if (numat <= 0)
    {
        fprintf(flog, "Invalid count of matrices...\n");
        exit(3);
    }
    /*-------------------- set parameters for arms */
    set_arms_pars(&io, diagscal, ipar, dropcoef, lfil_arr);
    /*-------------------- open file OUT/ARMS.out for all performance
                           results of this run (all matrices and params)
                           also set io->PrecMeth */
    /*-------------------- open file OUT/ARMS.out for all performance
      --- results of this run (all matrices) also set io->PrecMeth */
    strcpy(io.outfile,"OUT/ARMS.out");
    strcpy(io.PrecMeth,"ARMS");
    if (NULL == (io.fout = fopen(io.outfile, "w")))
    {
        fprintf(flog,"Can't open output file %s...\n", io.outfile);
        exit(4);
    }
    /*-------------------- LOOP through matrices -*/
    for (mat = 1; mat <= numat; mat++)
    {
        if (get_matrix_info(fmat, &io) != 0)
        {
            fprintf(flog, "Invalid format in matfile ...\n");
            exit(5);
        }
        fprintf(flog, "MATRIX: %s...\n", io.MatNam);
        /*-------------------- Read matrix - case: COO formats */
        csmat = (csptr)Malloc( sizeof(SparMat), "main:csmat" );
        if (io.Fmt > HB)
        {
            ierr = read_coo(&AA,&JA, &IA, &io, &rhs, &sol,0);
            if (ierr == 0)
                fprintf(flog,"matrix read successfully\n");
            else
            {
                fprintf(flog, "read_coo error = %d\n", ierr);
                exit(6);
            }
            n = io.ndim;
            nnz = io.nnz;
            /*-------------------- conversion from COO to CSR format */
            if ((ierr = COOcs(n, nnz, AA, JA, IA, csmat)) != 0)
            {
                fprintf(stderr, "mainARMS: COOcs error\n");
                return ierr;
            }
        }
        /*-------------------- Read matrix - case: HB formats */
        if (io.Fmt == HB)
        {
            /* NOTE: (AA,JA,IA) is in CSR format */
            ierr = readhb_c(&n, &AA, &JA, &IA, &io, &rhs, &sol, &rsa);
            if (ierr != 0)
            {
                fprintf(flog, "readhb_c error = %d\n", ierr);
                exit(7);
            }
            nnz = io.nnz;
            if ((ierr = CSRcs(n, AA, JA, IA, csmat, rsa)) != 0)
            {
                fprintf(flog, "readhb_c: CSRcs error\n");
                return ierr;
            }
        }
        /*---------------- COO / HB arrays no longer needed -- free */
        free(IA);
        IA = NULL;
        free(AA);
        AA = NULL;
        free(JA);
        JA = NULL;
        /*----------------------------------------------------------*/
        n = csmat->n;
        x = (double *)Malloc(n* sizeof(double), "main:x");
        output_header(&io);
        /*-------------------- set initial lfil and tol */
        lfil = io.lfil0;
        tol  = io.tol0;
        /*-------------------- LOOP THROUGH PARAMETERS */
        for (iparam = 1; iparam <= io.nparam; iparam++)
        {
            fprintf(flog, "Parameter case = %d\n", iparam);
            for (j=0; j<7; j++)
            {
                lfil_arr[j] = lfil*((int) nnz/n);
                droptol[j] = tol*dropcoef[j];
            }

            ArmsSt = (arms) Malloc(sizeof(armsMat),"main:ArmsSt");
            setup_arms(ArmsSt);
            fprintf(flog, "begin arms\n");
            tm1 = sys_timer();
            /*-------------------- call ARMS preconditioner set-up  */
            ierr =
                arms2(csmat, ipar, droptol, lfil_arr, tolind, ArmsSt, flog);
            /*----------------------------------------------------- */
            tm2 = sys_timer();
            if (ierr != 0)
            {
                fprintf(io.fout, " ** ARMS2 error - code %d...\n",ierr);
                io.its = -1;
                io.tm_i = -1;
                io.enorm = -1;
                io.rnorm = -1;
                goto NEXT_PARA;
            }
            io.tm_p = tm2 - tm1;
            io.fillfact = (double)nnz_arms(ArmsSt, flog)/(double)(nnz + 1);
            fprintf(flog, "ARMS ends, fill factor (mem used) = %f\n",\
                    io.fillfact);
            /*---------------- get rough idea of cond number - exit if too big */
            if (condestArms(ArmsSt, x, flog) != 0)
            {
                fprintf(flog, "Not attempting iterative solution.\n");
                fprintf(io.fout, "Not attempting iterative solution.\n");
                io.its = -1;
                io.tm_i = -1;
                io.enorm = -1;
                io.rnorm = -1;
                goto NEXT_PARA;
            }
            /*-------------------- initial guess */
            /*    for(i=0; i < io.ndim; i++) x[i] = 0.0 */
            randvec(x, n);
            /*-------------------- create a file for printing
                                   'its -- time -- res' info from fgmres */
            if (plotting)
            {
                sprintf(pltfile, "OUT/%s_ARMS_F%05d_T%08.6f", \
                        io.MatNam, lfil,tol);
                if (NULL == (fits = fopen(pltfile, "w")))
                {
                    fprintf(flog, "Can't open output file %s...\n", pltfile);
                    exit(8);
                }
            }
            else
                fits = NULL;
            /*-------------------- set up the structs before calling fgmr */
            MAT->n = n;
            MAT->CS = csmat;
            MAT->matvec = matvecCSR;
            PRE->ARMS = ArmsSt;
            PRE->precon = preconARMS;
            /*-------------------- call fgmr */
            io.its = io.maxits;
            tm1 = sys_timer();
            fgmr(MAT, PRE, rhs, x, io.tol, io.im, &io.its, fits);
            tm2 = sys_timer();
            io.tm_i = tm2 - tm1;
            if (io.its < io.maxits)
                fprintf(flog, "param %03d OK: converged in %d steps...\n", \
                        iparam, io.its);
            else
                fprintf(flog, "not converged in %d steps...\n", io.maxits);

            if (fits)
                fclose(fits);
            /*-------------------- actual error norm */
            terr = 0.0;
            for ( i = 0; i < n; i++)
                terr += ( x[i] - sol[i]) * ( x[i] - sol[i]);
            io.enorm = sqrt(terr);
            /*-------------------- calculate residual norm from generated rhs */
            matvec(csmat, x, sol);
            terr = 0.0;
            for (i = 0; i < io.ndim; i++)
                terr += (rhs[i] - sol[i]) * (rhs[i] - sol[i]);
            io.rnorm = sqrt(terr);
            /*-------------------- go to next param case */
NEXT_PARA:
            output_result(lfil, &io, iparam);
            lfil += io.lfilInc;
            tol  *= io.tolMul;
            cleanARMS(ArmsSt);
        }
        /*-------------------- NEXT_MATRIX: */
        cleanCS(csmat);
        free(sol);
        free(x);
        free(rhs);
    }

    fclose(io.fout);
    if (flog != stdout) fclose (flog);
    fclose(fmat);
    free(MAT);
    free(PRE);
    return 0;
}
