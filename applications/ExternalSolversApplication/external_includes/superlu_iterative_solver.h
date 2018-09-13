/*
* =======================================================================*
* kkkk   kkkk  kkkkkkkkkk   kkkkk    kkkkkkkkkk kkkkkkkkkk kkkkkkkkkK    *
* kkkk  kkkk   kkkk   kkkk  kkkkkk   kkkkkkkkkk kkkkkkkkkk kkkkkkkkkK    *
* kkkkkkkkk    kkkk   kkkk  kkkkkkk     kkkk    kkk    kkk  kkkk         *
* kkkkkkkkk    kkkkkkkkkkk  kkkk kkk	kkkk    kkk    kkk    kkkk       *
* kkkk  kkkk   kkkk  kkkk   kkkk kkkk   kkkk    kkk    kkk      kkkk     *
* kkkk   kkkk  kkkk   kkkk  kkkk  kkkk  kkkk    kkkkkkkkkk  kkkkkkkkkk   *
* kkkk    kkkk kkkk    kkkk kkkk   kkkk kkkk    kkkkkkkkkk  kkkkkkkkkk 	 *
*                                                                        *
* krATos: a fREe opEN sOURce CoDE for mULti-pHysIC aDaptIVe SoLVErS,     *
* aN extEnsIBLe OBjeCt oRiEnTEd SOlutION fOR fInITe ELemEnt fORmULatIONs *
* Copyleft by 2003 ciMNe                                                 *
* Copyleft by 2003 originary authors Copyleft by 2003 your name          *
* This library is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License as         *
* published by the Free Software Foundation; either version 2.1 of       *
* the License, or any later version.                                     *
*                                                                        *
* This library is distributed in the hope that it will be useful, but    *
* WITHOUT ANY WARRANTY; without even the implied warranty of             *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                   *
* See the GNU Lesser General Public License for more details.            *
*                                                                        *
* You should have received a copy of the GNU Lesser General Public       *
* License along with this library; if not, write to International Centre *
* for Numerical Methods in Engineering (CIMNE),                          *
* Edifici C1 - Campus Nord UPC, Gran Capità s/n, 08034 Barcelona.        *
*                                                                        *
* You can also contact us to the following email address:                *
* kratos@cimne.upc.es                                                    *
* or fax number: +34 93 401 65 17                                        *
*                                                                        *
* Created at Institute for Structural Mechanics                          *
* Ruhr-University Bochum, Germany                                        *
* Last modified by:    $Author: rrossi $  				 *
* Date:                $Date: 2009-01-15 11:11:35 $			 *
* Revision:            $Revision: 1.5 $ 				 *
*========================================================================*
* International Center of Numerical Methods in Engineering - CIMNE	 *
* Barcelona - Spain 							 *
*========================================================================*
*/

#if !defined(KRATOS_SUPERLU_ITERATIVE_SOLVER_H_INCLUDED )
#define  KRATOS_SUPERLU_ITERATIVE_SOLVER_H_INCLUDED

// External includes

#include "includes/ublas_interface.h"
// #include "external_includes/superlu/superlu.hpp"

#include "SRC/slu_ddefs.h"

// Project includes
#include "includes/define.h"
#include "linear_solvers/direct_solver.h"

namespace ublas = boost::numeric::ublas;

namespace Kratos
{

extern "C"
{
    void dcopy_(int *, double [], int *, double [], int *);
}


//some definitions for SuperLU iterative solver
char *GLOBAL_EQUED;
superlu_options_t *GLOBAL_OPTIONS;
double *GLOBAL_R, *GLOBAL_C;
int *GLOBAL_PERM_C, *GLOBAL_PERM_R;
SuperMatrix *GLOBAL_A, *GLOBAL_L, *GLOBAL_U;
// SuperMatrix *GLOBAL_A, *GLOBAL_A_ORIG, *GLOBAL_L, *GLOBAL_U;
SuperLUStat_t *GLOBAL_STAT;
mem_usage_t   *GLOBAL_MEM_USAGE;

void dmatvec_mult(double alpha, double x[], double beta, double y[])
{
    SuperMatrix *A = GLOBAL_A;

    char flag[] = "N";
    sp_dgemv(flag, alpha, A, x, 1, beta, y, 1);
}

void dpsolve(int n, double x[], double y[])
{

    int i_1 = 1;
    SuperMatrix *L = GLOBAL_L, *U = GLOBAL_U;
    SuperLUStat_t *stat = GLOBAL_STAT;
    int *perm_c = GLOBAL_PERM_C, *perm_r = GLOBAL_PERM_R;
    int info;
    static DNformat X;
    static SuperMatrix XX = {SLU_DN, SLU_D, SLU_GE, 1, 1, &X};

    dcopy_(&n, y, &i_1, x, &i_1);
    XX.nrow = n;
    X.lda = n;
    X.nzval = x;
    dgstrs(NOTRANS, L, U, perm_c, perm_r, &XX, stat, &info);
}


//some initial definitions
extern "C"
{
    int dfgmr( int n,
    void (*matvec_mult)(double, double [], double, double []),
    void (*psolve)(int n, double [], double[]),
    double *rhs, double *sol, double tol, int restrt, int *itmax,
    FILE *fits);
}
extern "C"
{
    int dfill_diag(int n, NCformat *Astore);
}

extern "C"
{
    double dnrm2_(int *, double [], int *);
}
extern "C"
{
    void daxpy_(int *, double *, double [], int *, double [], int *);
}

extern "C"
{
  void dCompRow_to_CompCol(int , int , int , double *, int *, int *, double **, int **, int **);
}















template< class TSparseSpaceType, class TDenseSpaceType,
class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType> >
class SuperLUIterativeSolver : public DirectSolver< TSparseSpaceType,
        TDenseSpaceType, TReordererType>
{
public:
    /**
     * Counted pointer of SuperLUIterativeSolver
     */
    KRATOS_CLASS_POINTER_DEFINITION(  SuperLUIterativeSolver );

    typedef LinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType> BaseType;

    typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

    typedef typename TSparseSpaceType::VectorType VectorType;

    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;
    
    SuperLUIterativeSolver(Parameters settings)
    {
        KRATOS_TRY
        
        
        Parameters default_parameters( R"(
        {
        "solver_type": "SuperLUIterativeSolver",
        "tolerance" : 1.0e-6,
        "max_iteration" : 200,
        "gmres_krylov_space_dimension": 100,
        "drop_tolerance":  1e-4  ,
        "fill_tolerance":  1e-2  ,
        "fill_factor"   :  10  ,
        "scaling":false
        }  )" );

        //now validate agains defaults -- this also ensures no type mismatch
        settings.ValidateAndAssignDefaults(default_parameters);
        
        mTol = settings["tolerance"].GetDouble();
        mmax_it = settings["max_iteration"].GetInt();
        mrestart = settings["gmres_krylov_space_dimension"].GetInt();
        mDropTol = settings["drop_tolerance"].GetDouble();
        mFillTol = settings["fill_tolerance"].GetDouble();
        mFillFactor = settings["fill_factor"].GetInt();
        
        KRATOS_CATCH("")
    }

    /**
     * Default constructor
     */
    SuperLUIterativeSolver(double NewMaxTolerance,
                           int NewMaxIterationsNumber,
                           int restart,
                           double DropTol=1e-4,
                           double FillTol = 1e-2,
                           double FillFactor=10)
    {
        mTol = NewMaxTolerance;
        mmax_it = NewMaxIterationsNumber;
        mrestart = restart;
        mDropTol = DropTol;
        mFillTol = FillTol;
        mFillFactor = FillFactor;
    }

    SuperLUIterativeSolver()
    {
        mTol = 1e-6;
        mrestart = 150;
        mmax_it = mrestart*3;
        mDropTol = 1e-4;
        mFillTol = 1e-2;
        mFillFactor = 10;
    }

    /**
     * Destructor
     */
    ~SuperLUIterativeSolver() override {};

    /**
     * Normal solve method.
     * Solves the linear system Ax=b and puts the result on SystemVector& rX.
     * rX is also th initial guess for iterative methods.
     * @param rA. System matrix
     * @param rX. Solution vector.
     * @param rB. Right hand side vector.
     */
    bool Solve(SparseMatrixType& rA, VectorType& rX, VectorType& rB) override
    {
        //void dmatvec_mult(double alpha, double x[], double beta, double y[]);
        //void dpsolve(int n, double x[], double y[]);
        //////////////////////////////////////////////////////////////////////////////extern int dfgmr( int n,
        //	void (*matvec_mult)(double, double [], double, double []),
        //	void (*psolve)(int n, double [], double[]),
        //	double *rhs, double *sol, double tol, int restrt, int *itmax,
        //	FILE *fits);
        //    extern int dfill_diag(int n, NCformat *Astore);

        char     equed[1] = {'B'};
        //yes_no_t equil;
//        trans_t  trans;

//         SuperMatrix A, AA, L, U;
        SuperMatrix A, L, U;
        SuperMatrix B, X;

        
	//double   *a_orig;
        //int       *asub_orig, *xa_orig;
        int      *etree;
        int      *perm_c; /* column permutation vector */
        int      *perm_r; /* row permutations from partial pivoting */
        int      /*nrhs,*/ lwork, info, m, n;
        double   *work = NULL;
        double   *R, *C;
        double   rpg, rcond;
        double zero = 0.0;
        //double one = 1.0;
        mem_usage_t   mem_usage;
        superlu_options_t options;
        SuperLUStat_t stat;

        int i;
        //int restrt, iter, maxit, i;
        double resid;
        double *x, *b;

#ifdef DEBUG
        extern int num_drop_L, num_drop_U;
#endif

#if ( DEBUGlevel>=1 )
        CHECK_MALLOC("Enter main()");
#endif

        /* Defaults */
        lwork = 0;
//        nrhs  = 1;
//        trans = NOTRANS;

        /* Set the default input options:
        options.Fact = DOFACT;
        options.Equil = YES;
        options.ColPerm = COLAMD;
        options.DiagPivotThresh = 0.1; //different from complete LU
        options.Trans = NOTRANS;
        options.IterRefine = NOREFINE;
        options.SymmetricMode = NO;
        options.PivotGrowth = NO;
        options.ConditionNumber = NO;
        options.PrintStat = YES;
        options.RowPerm = LargeDiag;
        options.ILU_DropTol = 1e-4;
        options.ILU_FillTol = 1e-2;
        options.ILU_FillFactor = 10.0;
        options.ILU_DropRule = DROP_BASIC | DROP_AREA;
        options.ILU_Norm = INF_NORM;
        options.ILU_MILU = SILU;
         */

	ilu_set_default_options(&options);
	options.ILU_MILU = SILU; //SMILU_3;
	 
	options.PrintStat = NO;
	options.Trans = NOTRANS;
	
// 	options.RowPerm = NO;


        options.ILU_DropTol = mDropTol;
        options.ILU_FillTol = mFillTol;
        options.ILU_FillFactor = mFillFactor; 
	
// 	options.SymmetricMode = YES;
//  	options.ILU_DropRule = DROP_DYNAMIC;
//         options.ILU_Norm = ONE_NORM;

// 	options.IterRefine = SLU_DOUBLE;

        /* Modify the defaults. */
         options.PivotGrowth = YES;	  /* Compute reciprocal pivot growth */
         options.ConditionNumber = YES;/* Compute reciprocal condition number */

        if ( lwork > 0 )
        {
            work =(double*) SUPERLU_MALLOC(lwork);
            if ( !work ) ABORT("Malloc fails for work[].");
        }


        //copy matrix from kratos
        //Fill the SuperLU matrices
        //NRformat *Astore;
        //NCformat *Ustore;
        //SCformat *Lstore;
        m = rA.size1();
        n = rA.size2();

        //create a copy of the matrix
        int *index1_vector = new (std::nothrow) int[rA.index1_data().size()];
        int *index2_vector = new (std::nothrow) int[rA.index2_data().size()];

        for ( int unsigned i = 0; i < rA.index1_data().size(); i++ )
            index1_vector[i] = (int)rA.index1_data()[i];

        for ( unsigned int i = 0; i < rA.index2_data().size(); i++ )
            index2_vector[i] = (int)rA.index2_data()[i];

        //works also with dCreate_CompCol_Matrix
//         dCreate_CompCol_Matrix (&A, rA.size1(), rA.size2(),
//                                 rA.nnz(),
//                                 rA.value_data().begin(),
//                                 index2_vector, //can not avoid a copy as ublas uses unsigned int internally
//                                 index1_vector, //can not avoid a copy as ublas uses unsigned int internally
//                                 SLU_NC, SLU_D, SLU_GE //ARGH...maybe this should be SLU_NC
//                                );

// 	SuperMatrix Akratos;
// 	dCreate_CompRow_Matrix (&A, rA.size1(), rA.size2(),
//                                 rA.nnz(),
//                                 rA.value_data().begin(),
//                                 index2_vector, //can not avoid a copy as ublas uses unsigned int internally
//                                 index1_vector, //can not avoid a copy as ublas uses unsigned int internally
//                                 SLU_NR, SLU_D, SLU_GE
//                                );
	
	int_t *asubt, *xat;
	double  *at;
	dCompRow_to_CompCol(rA.size1(), rA.size2(), rA.nnz(),rA.value_data().begin(), index2_vector, index1_vector,&at, &asubt, &xat);
				
	dCreate_CompCol_Matrix(&A, rA.size1(), rA.size2(),rA.nnz(),  at, asubt, xat, SLU_NC, SLU_D, SLU_GE);
	
        delete [] index1_vector;
        delete [] index2_vector;

	

	
// 	double **values;;
// 	int **i1;
//         int **i2;	       
// 	dCompRow_to_CompCol(rA.size1(), rA.size2(), rA.nnz(), 
// 		    rA.value_data().begin(), index2_vector, index1_vector,
// 		    values, i1, i2);

	//jmc
	//NCformat *Astore;
        //Astore = (NCformat*) A.Store;
        dfill_diag(n, (NCformat*) A.Store);

	//jmc
        //printf("Dimension %dx%d; # nonzeros %d\n", A.nrow, A.ncol, ((NCformat*) A.Store)->nnz);

        fflush(stdout);

        /* Make a copy of the original matrix. */
//         nnz = Astore->nnz;
//         a_orig = doubleMalloc(nnz);
//         asub_orig = intMalloc(nnz);
//         xa_orig = intMalloc(n+1);
//         for (i = 0; i < nnz; ++i)
//         {
//             a_orig[i] = ((double *)Astore->nzval)[i];
//             asub_orig[i] = Astore->rowind[i];
//         }
//         for (i = 0; i <= n; ++i) xa_orig[i] = Astore->colptr[i];
//         dCreate_CompCol_Matrix(&AA, m, n, nnz, a_orig, asub_orig, xa_orig,
//                                SLU_NR, SLU_D, SLU_GE);

        /* Generate the right-hand side */
        // if ( !(rhsb = doubleMalloc(m * nrhs)) ) ABORT("Malloc fails for rhsb[].");
        // if ( !(rhsx = doubleMalloc(m * nrhs)) ) ABORT("Malloc fails for rhsx[].");
        dCreate_Dense_Matrix (&B, rB.size(), 1,&rB[0],rB.size(),SLU_DN, SLU_D, SLU_GE);
        dCreate_Dense_Matrix (&X, rX.size(), 1,&rX[0],rX.size(),SLU_DN, SLU_D, SLU_GE);
        //dCreate_Dense_Matrix(&B, m, nrhs, rhsb, m, SLU_DN, SLU_D, SLU_GE);
        //dCreate_Dense_Matrix(&X, m, nrhs, rhsx, m, SLU_DN, SLU_D, SLU_GE);
        //xact = doubleMalloc(n * nrhs);
        //ldx = n;
        //dGenXtrue(n, nrhs, xact, ldx);
        //dFillRHS(trans, nrhs, xact, ldx, &A, &B);

        if ( !(etree = intMalloc(n)) ) ABORT("Malloc fails for etree[].");
        if ( !(perm_r = intMalloc(m)) ) ABORT("Malloc fails for perm_r[].");
        if ( !(perm_c = intMalloc(n)) ) ABORT("Malloc fails for perm_c[].");
        if ( !(R = (double *) SUPERLU_MALLOC(A.nrow * sizeof(double))) )
            ABORT("SUPERLU_MALLOC fails for R[].");
        if ( !(C = (double *) SUPERLU_MALLOC(A.ncol * sizeof(double))) )
            ABORT("SUPERLU_MALLOC fails for C[].");

        info = 0;
#ifdef DEBUG
        num_drop_L = 0;
        num_drop_U = 0;
#endif

        /* Initialize the statistics variables. */
        StatInit(&stat);

        /* Compute the incomplete factorization and compute the condition number
           and pivot growth using dgsisx. */
        B.ncol = 0;  /* not to perform triangular solution */
        dgsisx(&options, &A, perm_c, perm_r, etree, equed, R, C, &L, &U, work,
               lwork, &B, &X, &rpg, &rcond, &mem_usage, &stat, &info);

        /* Set RHS for GMRES. */
        if (!(b = doubleMalloc(m))) ABORT("Malloc fails for b[].");
        for (i = 0; i < m; i++) b[i] = rB[i];
	
	//jmc
        //printf("dgsisx(): info %d, equed %c\n", info, equed[0]);
        //if (info > 0 || rcond < 1e-8 || rpg > 1e8)
        //    printf("WARNING: This preconditioner might be unstable.\n");

        //if ( info == 0 || info == n+1 )
        //{
        //    if ( options.PivotGrowth == YES )
        //        printf("Recip. pivot growth = %e\n", rpg);
        //    if ( options.ConditionNumber == YES )
        //        printf("Recip. condition number = %e\n", rcond);
        //}
        //else if ( info > 0 && lwork == -1 )
        //{
        //    printf("** Estimated memory: %d bytes\n", info - n);
	//}

	//jmc
        //NCformat *Ustore;
        //SCformat *Lstore;
        //Lstore = (SCformat *) L.Store;
        //Ustore = (NCformat *) U.Store;
        //printf("n(A) = %d, nnz(A) = %d\n", n, Astore->nnz);
        //printf("No of nonzeros in factor L = %d\n", Lstore->nnz);
        //printf("No of nonzeros in factor U = %d\n", Ustore->nnz);
        //printf("No of nonzeros in L+U = %d\n", Lstore->nnz + Ustore->nnz - n);
        //printf("Fill ratio: nnz(F)/nnz(A) = %.3f\n",
        //       ((double)(Lstore->nnz) + (double)(Ustore->nnz) - (double)n)
        //       / (double)Astore->nnz);
        //printf("L\\U MB %.3f\ttotal MB needed %.3f\n",
        //       mem_usage.for_lu/1e6, mem_usage.total_needed/1e6);
        //fflush(stdout);

        /* Set the global variables. */
        GLOBAL_A = &A;
//         GLOBAL_A_ORIG = &AA;
        GLOBAL_L = &L;
        GLOBAL_U = &U;
        GLOBAL_STAT = &stat;
        GLOBAL_PERM_C = perm_c;
        GLOBAL_PERM_R = perm_r;
        GLOBAL_OPTIONS = &options;
        GLOBAL_EQUED = equed;
        GLOBAL_R = R;
        GLOBAL_C = C;
        GLOBAL_MEM_USAGE = &mem_usage;

        /* Set the options to do solve-only. */
        options.Fact = FACTORED;
        options.PivotGrowth = NO;
        options.ConditionNumber = NO;

        /* Set the variables used by GMRES. */
        //restrt = SUPERLU_MIN(n / 3 + 1, 150);
        int restrt = SUPERLU_MIN(n / 3 + 1, mrestart);
        //jmc
	//int maxit = mmax_it; //1000;
        int iter = mmax_it;
        resid = mTol; // 1e-8;
        if (!(x = doubleMalloc(n))) ABORT("Malloc fails for x[].");

        if (info <= n + 1)
        {
            int i_1 = 1;

	    double nrmB, res, t;
            
            //extern double dnrm2_(int *, double [], int *);
            //extern void daxpy_(int *, double *, double [], int *, double [], int *);

            /* Initial guess */
            for (i = 0; i < n; i++) x[i] = zero;

            t = SuperLU_timer_();

            /* Call GMRES */
            //dfgmr(n, dmatvec_mult, dpsolve, b, x, resid, restrt, &iter, stdout);
            dfgmr(n, dmatvec_mult, dpsolve, b, x, resid, restrt, &iter, NULL);

            /* Scale the solution back if equilibration was performed. */
            if (*equed == 'C' || *equed == 'B')
            {
                #pragma omp parallel for
                for (int i = 0; i < n; i++) x[i] *= C[i];
            }
            
            t = SuperLU_timer_() - t;

            /* Output the result. */
	    //jmc
	    //double nrmA;
            //nrmA = dnrm2_(&(Astore->nnz), (double *)((DNformat *)A.Store)->nzval,
            //              &i_1);

            nrmB = dnrm2_(&m, b, &i_1);
            sp_dgemv( ( char *)"N", -1.0, &A, x, 1, 1.0, b, 1);
            res = dnrm2_(&m, b, &i_1);
            resid = res / nrmB;
	    //jmc
            //printf("||A||_F = %.1e, ||B||_2 = %.1e, ||B-A*X||_2 = %.1e, "
            //       "relres = %.1e\n", nrmA, nrmB, res, resid);

            //if (iter >= maxit)
            //{
            //    if (resid >= 1.0){
	    //      iter = -180;
	    //      std::cout << "final residual is VERY VERY LARGE" << std::endl;
	    //}
            //    else if (resid > 1e-8) iter = -111;
            //}
            //printf("iteration: %d\nresidual: %.1e\nGMRES time: %.2f seconds.\n",
            //       iter, resid, t);

            //for (i = 0; i < m; i++) {
            //  maxferr = SUPERLU_MAX(maxferr, fabs(x[i] - xact[i]));
            //}
            //printf("||X-X_true||_oo = %.1e\n", maxferr);
        }
        fflush(stdout);

#pragma omp parallel for
        for (int i = 0; i < n; i++) rX[i] = x[i];

        if ( options.PrintStat ) StatPrint(&stat);
        StatFree(&stat);

        //SUPERLU_FREE (rhsb);
        //SUPERLU_FREE (rhsx);
        //SUPERLU_FREE (xact);

        SUPERLU_FREE (etree);

        SUPERLU_FREE (perm_r);

        SUPERLU_FREE (perm_c);

        SUPERLU_FREE (R);

        SUPERLU_FREE (C);

        Destroy_SuperMatrix_Store(&A);

//         Destroy_CompCol_Matrix(&AA);

        Destroy_SuperMatrix_Store(&B);

        Destroy_SuperMatrix_Store(&X);

        if ( lwork >= 0 )
        {

            Destroy_SuperNode_Matrix(&L);

            Destroy_CompCol_Matrix(&U);
        }

        SUPERLU_FREE(b);

        SUPERLU_FREE(x);



	SUPERLU_FREE(asubt);
	SUPERLU_FREE(xat);
	SUPERLU_FREE(at);
#if ( DEBUGlevel>=1 )
        CHECK_MALLOC("Exit main()");
#endif


        return true;
    }

    /**
     * Multi solve method for solving a set of linear systems with same coefficient matrix.
     * Solves the linear system Ax=b and puts the result on SystemVector& rX.
     * rX is also th initial guess for iterative methods.
     * @param rA. System matrix
     * @param rX. Solution vector.
     * @param rB. Right hand side vector.
     */
    bool Solve(SparseMatrixType& rA, DenseMatrixType& rX, DenseMatrixType& rB) override
    {
        /**
         * TODO:
             * translate SparseMatrixType into SuperMatrix
             * call solving routine from SuperLU
         */
//                 slu::gssv ( rA, rB, slu::atpla_min_degree);

//                 std::cout<<"Matrix Test:"<<std::endl;
//                 std::cout<<"boost matrix:"<<std::endl;
//                 KRATOS_WATCH( rA );
        //             const int size1 = TDenseSpaceType::Size1(rX);
//             const int size2 = TDenseSpaceType::Size2(rX);

        bool is_solved = true;

//             VectorType x(size1);
//             VectorType b(size1);

        // define an object to store skyline matrix and factorization
//             LUSkylineFactorization<TSparseSpaceType, TDenseSpaceType> myFactorization;
        // copy myMatrix into skyline format
//             myFactorization.copyFromCSRMatrix(rA);
        // factorize it
//             myFactorization.factorize();

//             for(int i = 0 ; i < size2 ; i++)
//             {
//                 TDenseSpaceType::GetColumn(i,rX, x);
//                 TDenseSpaceType::GetColumn(i,rB, b);

        // and back solve
//                 myFactorization.backForwardSolve(size1, b, x);

//                 TDenseSpaceType::SetColumn(i,rX, x);
//                 TDenseSpaceType::SetColumn(i,rB, b);
//             }

        return is_solved;
    }

    /**
     * Print information about this object.
     */
    void  PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "SuperLU solver finished.";
    }

    /**
     * Print object's data.
     */
    void  PrintData(std::ostream& rOStream) const override
    {
    }

private:

    double mTol; 
    int mmax_it;
    int mrestart;
    double mDropTol;
    double mFillTol;
    double mFillFactor;

    /**
     * Assignment operator.
     */
    SuperLUIterativeSolver& operator=(const SuperLUIterativeSolver& Other);

    /**
     * Copy constructor.
     */
    SuperLUIterativeSolver(const SuperLUIterativeSolver& Other);

}; // Class SkylineLUFactorizationSolver


/**
 * input stream function
 */
template<class TSparseSpaceType, class TDenseSpaceType,class TReordererType>
inline std::istream& operator >> (std::istream& rIStream, SuperLUIterativeSolver< TSparseSpaceType,
                                  TDenseSpaceType, TReordererType>& rThis)
{
    return rIStream;
}

/**
 * output stream function
 */
template<class TSparseSpaceType, class TDenseSpaceType, class TReordererType>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const SuperLUIterativeSolver<TSparseSpaceType,
                                  TDenseSpaceType, TReordererType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}


}  // namespace Kratos.

#endif // KRATOS_SUPERLU_ITERATIVE_SOLVER_H_INCLUDED  defined 


