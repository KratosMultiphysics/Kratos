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
#include "boost/smart_ptr.hpp"

#include "includes/ublas_interface.h"
// #include "external_includes/superlu/superlu.hpp"

#include "SRC/slu_ddefs.h"

// Project includes
#include "includes/define.h"
#include "linear_solvers/direct_solver.h"

namespace ublas = boost::numeric::ublas;

namespace Kratos
{
 
  extern "C" {
      void dcopy_(int *, double [], int *, double [], int *);
  }
  
  
    //some definitions for SuperLU iterative solver
    int *GLOBAL_PERM_C, *GLOBAL_PERM_R;
    SuperMatrix *GLOBAL_A, *GLOBAL_L, *GLOBAL_U;
    SuperLUStat_t *GLOBAL_STAT;

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
    extern "C" {int dfgmr( int n,
	void (*matvec_mult)(double, double [], double, double []),
	void (*psolve)(int n, double [], double[]),
	double *rhs, double *sol, double tol, int restrt, int *itmax,
	FILE *fits);}
    extern "C" { int dfill_diag(int n, NCformat *Astore); }

    extern "C" {double dnrm2_(int *, double [], int *); }
    extern "C" {void daxpy_(int *, double *, double [], int *, double [], int *);}
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
    template< class TSparseSpaceType, class TDenseSpaceType, 
              class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType> >
        class SuperLUIterativeSolver : public DirectSolver< TSparseSpaceType, 
              TDenseSpaceType, TReordererType>
    {
        public:
            /**
             * Counted pointer of SuperLUIterativeSolver
             */
            typedef boost::shared_ptr<SuperLUIterativeSolver> Pointer;
            
            typedef LinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType> BaseType; 
            
            typedef typename TSparseSpaceType::MatrixType SparseMatrixType;
            
            typedef typename TSparseSpaceType::VectorType VectorType;
            
            typedef typename TDenseSpaceType::MatrixType DenseMatrixType;
            
            /**
             * Default constructor
             */
            SuperLUIterativeSolver(double NewMaxTolerance, 
				   int NewMaxIterationsNumber, 
				   int restart, 
				   double DropTol=1e-6, 
				   double FillTol = 1e-4, 
				   double FillFactor=20)
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
	      mrestart = 50;
	      mmax_it = mrestart*3;
	      mDropTol = 1e-6;
	      mFillTol = 1e-4;
	      mFillFactor = 20;
	    }
            
            /**
             * Destructor
             */
	virtual ~SuperLUIterativeSolver(){};
            
            /** 
             * Normal solve method.
             * Solves the linear system Ax=b and puts the result on SystemVector& rX.
             * rX is also th initial guess for iterative methods.
             * @param rA. System matrix
             * @param rX. Solution vector.
             * @param rB. Right hand side vector.
             */
            bool Solve(SparseMatrixType& rA, VectorType& rX, VectorType& rB)
            {
		std::cout << "matrix size in solver: " << rA.size1() << std::endl;
		std::cout << "RHS size in solver: " << rB.size() << std::endl;
                               
                if(this->IsNotConsistent(rA, rX, rB))
                    return false;
		
		int m = rA.size1();
		int n = rA.size2();
		double zero = 0.0;
		//double one = 1.0;
		

		char     equed[1] = {'B'};
		/* Defaults */
		int lwork = 0;
		//int nrhs  = 1;
		//yes_no_t equil = YES;
		//double u	  = 0.1; /* u=1.0 for complete factorization */
		//trans_t trans = NOTRANS;
		    
                superlu_options_t options;
		SuperLUStat_t stat;

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
		//set options for ILU
		ilu_set_default_options(&options);
        options.ColPerm = MMD_ATA;
		options.PivotGrowth = YES;	  /* Compute reciprocal pivot growth */
		options.ConditionNumber = YES;/* Compute reciprocal condition number */
		options.ILU_DropTol = mDropTol;
 		options.ILU_FillTol = mFillTol;
		options.ILU_FillFactor = mFillFactor;
		
		double   *work = NULL;
		if ( lwork > 0 ) {
		  work = (double*)SUPERLU_MALLOC(lwork);
		    if ( !work ) ABORT("Malloc fails for work[].");
		}
		
                
                //Fill the SuperLU matrices
                SuperMatrix Aslu, B, L, U, X;
		//NRformat *Astore;
		NCformat *Ustore;
		SCformat *Lstore;
		
		//create a copy of the matrix
		int *index1_vector = new (std::nothrow) int[rA.index1_data().size()];
		int *index2_vector = new (std::nothrow) int[rA.index2_data().size()];
		
		for( int unsigned i = 0; i < rA.index1_data().size(); i++ )
		  index1_vector[i] = (int)rA.index1_data()[i];

		for( unsigned int i = 0; i < rA.index2_data().size(); i++ )
		    index2_vector[i] = (int)rA.index2_data()[i];
		
		//works also with dCreate_CompCol_Matrix
		dCreate_CompCol_Matrix (&Aslu, rA.size1(), rA.size2(), 
				rA.nnz(),
			      rA.value_data().begin(), 
			      index2_vector, //can not avoid a copy as ublas uses unsigned int internally
			      index1_vector, //can not avoid a copy as ublas uses unsigned int internally
			      SLU_NR, SLU_D, SLU_GE
			      );
		dfill_diag(rA.size1(),(NCformat*) Aslu.Store);
				
		
		dCreate_Dense_Matrix (&B, rB.size(), 1,&rB[0],rB.size(),SLU_DN, SLU_D, SLU_GE);   
		dCreate_Dense_Matrix (&X, rX.size(), 1,&rX[0],rX.size(),SLU_DN, SLU_D, SLU_GE);   
		
		//allocate memory for permutation arrays
		int *perm_c;
		int *perm_r;
		int *etree;
		double   *R, *C;
		if ( !(perm_c = intMalloc(rA.size1())) ) ABORT("Malloc fails for perm_c[].");
		if ( !(perm_r = intMalloc(rA.size2())) ) ABORT("Malloc fails for perm_r[].");
		if ( !(etree = intMalloc(n)) ) ABORT("Malloc fails for etree[].");
		if ( !(R = (double *) SUPERLU_MALLOC(Aslu.nrow * sizeof(double))) )
		    ABORT("SUPERLU_MALLOC fails for R[].");
		if ( !(C = (double *) SUPERLU_MALLOC(Aslu.ncol * sizeof(double))) )
		    ABORT("SUPERLU_MALLOC fails for C[].");		

		//initialize container for statistical data
		StatInit(&stat);

		//call solver routine
		int info = 0;
		double   rpg, rcond;
		mem_usage_t   mem_usage;
		/* Compute the incomplete factorization and compute the condition number
		  and pivot growth using dgsisx. */
		dgsisx(&options, &Aslu, perm_c, perm_r, etree, equed, R, C, &L, &U, work,
		      lwork, &B, &X, &rpg, &rcond, &mem_usage, &stat, &info);

		Lstore = (SCformat *) L.Store;
		Ustore = (NCformat *) U.Store;
		printf("dgsisx(): info %d\n", info);
		if (info > 0 || rcond < 1e-8 || rpg > 1e8)
		    printf("WARNING: This preconditioner might be unstable.\n");

		if ( info == 0 || info == n+1 ) {

		    if ( options.PivotGrowth == YES )
			printf("Recip. pivot growth = %e\n", rpg);
		    if ( options.ConditionNumber == YES )
			printf("Recip. condition number = %e\n", rcond);

		} else if ( info > 0 && lwork == -1 ) {
		    printf("** Estimated memory: %d bytes\n", info - n);
		}
 		printf("n(A) = %d, nnz(A) = %d\n", n, (((NRformat *)Aslu.Store)->nnz));
		printf("No of nonzeros in factor L = %d\n", Lstore->nnz);
		printf("No of nonzeros in factor U = %d\n", Ustore->nnz);
		printf("No of nonzeros in L+U = %d\n", Lstore->nnz + Ustore->nnz - n);
 		printf("Fill ratio: nnz(F)/nnz(A) = %.3f\n",
 			((double)(Lstore->nnz) + (double)(Ustore->nnz) - (double)n)
 			/ (double)(((NRformat *)Aslu.Store)->nnz) );
		printf("L\\U MB %.3f\ttotal MB needed %.3f\n",
		      mem_usage.for_lu/1e6, mem_usage.total_needed/1e6);
		fflush(stdout);

		/* Set the global variables. */
		GLOBAL_A = &Aslu;
		GLOBAL_L = &L;
		GLOBAL_U = &U;
		GLOBAL_STAT = &stat;
		GLOBAL_PERM_C = perm_c;
		GLOBAL_PERM_R = perm_r;
		


		/* Set the variables used by GMRES. */
		int restrt = SUPERLU_MIN(n / 3 + 1, mrestart);
		int maxit = mmax_it; //1000;
		int iter = mmax_it;
		double resid = mTol; // 1e-8;
		double *b;
		double *x;
		if (!(b = doubleMalloc(m))) ABORT("Malloc fails for b[].");
		if (!(x = doubleMalloc(n))) ABORT("Malloc fails for x[].");

		if (info <= n + 1)
		{
		    int i_1 = 1;
		    double nrmA, nrmB, res, t;
		    //double temp;

		    /* Call GMRES. */
		    #pragma omp parallel for
		    for (int i = 0; i < n; i++) b[i] = rB[i]; //it was rhsb
		    #pragma omp parallel for
		    for (int i = 0; i < n; i++) x[i] = zero;

		    t = SuperLU_timer_();

		    dfgmr(n, dmatvec_mult, dpsolve, b, x, resid, restrt, &iter, stdout);

		    t = SuperLU_timer_() - t;

		    /* Output the result. */
		    nrmA = dnrm2_(&(((NRformat *)Aslu.Store)->nnz), (double *)((DNformat *)Aslu.Store)->nzval,
			    &i_1);
		    nrmB = dnrm2_(&m, b, &i_1);

                    char flag[] = "N";
		    sp_dgemv(flag, -1.0, &Aslu, x, 1, 1.0, b, 1);
		    res = dnrm2_(&m, b, &i_1);
		    resid = res / nrmB;
		    printf("||A||_F = %.1e, ||B||_2 = %.1e, ||B-A*X||_2 = %.1e, "
			    "relres = %.1e\n", nrmA, nrmB, res, resid);

		    if (iter >= maxit)
		    {
			if (resid >= 1.0) iter = -180;
			else if (resid > 1e-8) iter = -111;
		    }
		    printf("iteration: %d\nresidual: %.1e\nGMRES time: %.2f seconds.\n",
			    iter, resid, t);

		    /* Scale the solution back if equilibration was performed. */
		    if (*equed == 'C' || *equed == 'B') 
			#pragma omp parallel for
			for (int i = 0; i < n; i++) x[i] *= C[i];

		}
	    #ifdef DEBUG
		printf("%d entries in L and %d entries in U dropped.\n",
			num_drop_L, num_drop_U);
	    #endif
		fflush(stdout);

		if ( options.PrintStat ) StatPrint(&stat);
		StatFree(&stat);
        		
		//deallocate memory used
		SUPERLU_FREE (perm_r);
		SUPERLU_FREE (perm_c);
		Destroy_SuperMatrix_Store(&Aslu); //note that by using the "store" function we will take care of deallocation ourselves
		Destroy_SuperMatrix_Store(&B);
		Destroy_SuperMatrix_Store(&X);
// 		Destroy_SuperNode_Matrix(&L);
// 		Destroy_CompCol_Matrix(&U);
		//	SUPERLU_FREE (rhsb);
		//SUPERLU_FREE (rhsx);
		SUPERLU_FREE (etree);
		SUPERLU_FREE (R);
		SUPERLU_FREE (C);
		if ( lwork >= 0 ) {
		    Destroy_SuperNode_Matrix(&L);
		    Destroy_CompCol_Matrix(&U);
		}
		SUPERLU_FREE(b);
		SUPERLU_FREE(x);
		
		delete [] index1_vector;
		delete [] index2_vector;
// 		delete [] b_vector;

		//CHECK WITH VALGRIND IF THIS IS NEEDED ...or if it is done by the lines above
                //deallocate tempory storage used for the matrix
//                 if(b_vector!=NULL) delete [] index1_vector;
// //   		if(b_vector!=NULL) delete [] index2_vector;
//   		if(b_vector!=NULL) delete [] values_vector;
// 		if(b_vector!=NULL) delete [] b_vector;
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
            bool Solve(SparseMatrixType& rA, DenseMatrixType& rX, DenseMatrixType& rB)
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
            void  PrintInfo(std::ostream& rOStream) const
            {
                rOStream << "SuperLU solver finished.";
            }
            
            /**
             * Print object's data.
             */
            void  PrintData(std::ostream& rOStream) const 
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


