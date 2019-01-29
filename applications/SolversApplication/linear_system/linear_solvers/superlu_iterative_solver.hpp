//
//   Project Name:        KratosSolversApplication $
//   Created by:          $Author:     JMCarbonell $
//   Last modified by:    $Co-Author:              $
//   Date:                $Date:      January 2019 $
//   Revision:            $Revision:           0.0 $
//
//

#if !defined(KRATOS_SUPERLU_ITERATIVE_GMRES_SOLVER_H_INCLUDED )
#define  KRATOS_SUPERLU_ITERATIVE_GMRES_SOLVER_H_INCLUDED

// System includes

// External includes
#include "includes/ublas_interface.h"

extern "C"{
   #include "slu_ddefs.h"
}

// Project includes
#include "linear_solvers/linear_solver.h"
#include "utilities/openmp_utils.h"

namespace Kratos
{

///@name Kratos Classes
///@{

//some definitions for SuperLU iterative solver
char *GLOBAL_EQUED;
superlu_options_t *GLOBAL_OPTIONS;
double *GLOBAL_R, *GLOBAL_C;
int *GLOBAL_PERM_C, *GLOBAL_PERM_R;
SuperMatrix *GLOBAL_A, *GLOBAL_L, *GLOBAL_U;
SuperLUStat_t *GLOBAL_STAT;
mem_usage_t   *GLOBAL_MEM_USAGE;


//some initial definitions
extern "C"
{
  void dcopy_(int *, double [], int *, double [], int *);

  int dfgmr(int n,
            void (*matvec_mult)(double, double [], double, double []),
            void (*psolve)(int n, double [], double[]),
            double *rhs, double *sol, double tol, int restrt, int *itmax,
            FILE *fits);

  int dfill_diag(int n, NCformat *Astore);

  double dnrm2_(int *, double [], int *);

  void daxpy_(int *, double *, double [], int *, double [], int *);
}

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

///@name Kratos Classes
///@{

template< class TSparseSpaceType, class TDenseSpaceType, class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType> >
class SuperLUIterativeGMRESSolver : public LinearSolver< TSparseSpaceType, TDenseSpaceType, TReordererType>
{
 public:
  ///@name Kratos Globals
  ///@{

  ///@}
  ///@name Type Definitions
  ///@{

  /// Pointer definition for the Solver
  KRATOS_CLASS_POINTER_DEFINITION( SuperLUIterativeGMRESSolver );

  typedef LinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType> BaseType;

  typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

  typedef typename TSparseSpaceType::VectorType VectorType;

  typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

  ///@}
  ///@name  Enum's
  ///@{

  ///@}
  ///@name  Functions
  ///@{

  ///@}
  ///@name Kratos Classes
  ///@{

  /// Empty constructor.
  SuperLUIterativeGMRESSolver() {}

  /// Default constructor.
  SuperLUIterativeGMRESSolver(Parameters settings)
  {
    KRATOS_TRY

    Parameters default_parameters( R"(
    {
        "solver_type"     : "superlu_iterative",
        "krylov_dimension": 100,
        "drop_tolerance"  : 1e-4,
        "fill_tolerance"  : 1e-2,
        "fill_factor"     : 10,
        "tolerance"       : 1e-8,
        "max_iteration"   : 1000,
        "scaling"         : false,
        "verbosity"       : 1
    }  )" );

    //now validate agains defaults -- this also ensures no type mismatch
    settings.ValidateAndAssignDefaults(default_parameters);

    mTolerance = settings["tolerance"].GetDouble();
    mMaxIteration = settings["max_iteration"].GetInt();
    mRestart = settings["krylov_dimension"].GetInt();
    mDropTolerance = settings["drop_tolerance"].GetDouble();
    mFillTolerance = settings["fill_tolerance"].GetDouble();
    mFillFactor = settings["fill_factor"].GetInt();
    mEchoLevel = settings["verbosity"].GetInt();

    KRATOS_CATCH("")
  }

  /// Copy constructor.
  SuperLUIterativeGMRESSolver(const SuperLUIterativeGMRESSolver& Other) {}

  /// Destructor.
  ~SuperLUIterativeGMRESSolver() override {};

  ///@}
  ///@name Operators
  ///@{

  /// Assignment operator.
  SuperLUIterativeGMRESSolver& operator=(const SuperLUIterativeGMRESSolver& Other) {}

  ///@}
  ///@name Operations
  ///@{

  /// Solve
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
    //KRATOS_TRY

    // double start_define = OpenMPUtils::GetCurrentTime();

    char equed[1] = {'B'};
    SuperMatrix lhs, L, U;
    SuperMatrix rhs, X;

    GlobalLU_t Glu; /* facilitate multiple factorizations with
                       SamePattern_SameRowPerm */
    double rpg, rcond;
    mem_usage_t mem_usage;
    superlu_options_t options;

    double* work = NULL;

#ifdef DEBUG
    extern int num_drop_L, num_drop_U;
#endif

    /* Defaults */
    int lwork = 0;

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

    /* Modify the defaults. */
    options.ILU_MILU = SILU;
    options.PrintStat = NO;
    options.Trans = NOTRANS;

    options.ILU_DropTol = mDropTolerance;
    options.ILU_FillTol = mFillTolerance;
    options.ILU_FillFactor = mFillFactor;
    options.PivotGrowth = YES;	  /* Compute reciprocal pivot growth */
    options.ConditionNumber = YES;/* Compute reciprocal condition number */
    //options.SymmetricMode = YES;
    //options.IterRefine = SLU_DOUBLE;

    /* Make a copy of the original matrix. */

    //copy matrix from kratos and fill the SuperLU matrices
    int n =rA.size1();
    int m =rA.size2();

    //create a copy of the matrix
    int* index1_vector = new (std::nothrow) int[rA.index1_data().size()];
    for( int unsigned i = 0; i < rA.index1_data().size(); i++ )
      index1_vector[i] = (int)rA.index1_data()[i];

    //can not avoid a copy as ublas uses unsigned int internally
    int* index2_vector = new (std::nothrow) int[rA.index2_data().size()];
    for( unsigned int i = 0; i < rA.index2_data().size(); i++ )
      index2_vector[i] = (int)rA.index2_data()[i];

    //double* values_vector = new (std::nothrow) double[rA.value_data().size()];
    //for( unsigned int i = 0; i < rA.value_data().size(); i++ )
    //  values_vector[i] = (double)rA.value_data()[i];

    // Creation of a Column Matrix (option0)
    //dCreate_CompCol_Matrix (&lhs, n, m, rA.nnz(), values_vector, index2_vector, index1_vector, SLU_NR, SLU_D, SLU_GE );

    // Creation of a Column Matrix (option1)
    //dCreate_CompCol_Matrix (&lhs, n, m, rA.nnz(), rA.value_data().begin(), index2_vector, index1_vector, SLU_NR, SLU_D, SLU_GE );

    // Creation of a Row Matrix (option2)
    dCreate_CompRow_Matrix (&lhs, n, m, rA.nnz(), rA.value_data().begin(), index2_vector, index1_vector, SLU_NR, SLU_D, SLU_GE);

    // Creation of a Row Matrix (option3)
    // int_t *asubt, *xat;
    // double  *at;
    // dCompRow_to_CompCol(n, m, rA.nnz(), rA.value_data().begin(), index2_vector, index1_vector, &at, &asubt, &xat);

    // dCreate_CompCol_Matrix(&lhs, n, m, rA.nnz(), at, asubt, xat, SLU_NC, SLU_D, SLU_GE);

    dfill_diag(n, (NCformat*) lhs.Store);

    /* Generate the right-hand side */
    dCreate_Dense_Matrix (&rhs, m, 1, &rB[0], rB.size(), SLU_DN, SLU_D, SLU_GE);
    dCreate_Dense_Matrix (&X, m, 1, &rX[0], rX.size(), SLU_DN, SLU_D, SLU_GE);

    //allocate memory for permutation arrays
    int* etree;
    int* perm_c; /* column permutation vector */
    int* perm_r; /* row permutations from partial pivoting */
    if ( !(etree = intMalloc(n)) ) ABORT("Malloc fails for etree[].");
    if ( !(perm_c = intMalloc(n)) ) ABORT("Malloc fails for perm_c[].");
    if ( !(perm_r = intMalloc(m)) ) ABORT("Malloc fails for perm_r[].");
    //allocate memory for other arrays
    double* R;
    double* C;
    if ( !(R = (double *) SUPERLU_MALLOC(lhs.nrow * sizeof(double))) )
      ABORT("SUPERLU_MALLOC fails for R[].");
    if ( !(C = (double *) SUPERLU_MALLOC(lhs.ncol * sizeof(double))) )
      ABORT("SUPERLU_MALLOC fails for C[].");

    // double stop_define = OpenMPUtils::GetCurrentTime();
    // if(mEchoLevel > 2)
    //   KRATOS_INFO("superlu_define_time") << stop_define - start_define << std::endl;

    //initialize container for statistical data
    SuperLUStat_t stat;
    int info = 0;

    // double start_factorization = OpenMPUtils::GetCurrentTime();

    // Resolution of the linear system:

    /* Initialize the statistics variables. */
    StatInit(&stat);

    /* Compute the incomplete factorization and compute the condition number
       and pivot growth using dgsisx. */
    rhs.ncol = 0;  /* not to perform triangular solution */
    dgsisx(&options, &lhs, perm_c, perm_r, etree, equed, R, C, &L, &U, work,
           lwork, &rhs, &X, &rpg, &rcond, &Glu, &mem_usage, &stat, &info);


    // double stop_factorization = OpenMPUtils::GetCurrentTime();
    // if(mEchoLevel > 2)
    //   KRATOS_INFO("superlu_factorization_time") << stop_factorization - start_factorization << std::endl;

    // double start_solve = OpenMPUtils::GetCurrentTime();

    /* Set RHS for GMRES. */
    double* b;
    if (!(b = doubleMalloc(m))) ABORT("Malloc fails for b[].");
    for (int i=0; i<m; ++i)
      b[i] = rB[i];

    /* Set the global variables. */
    GLOBAL_A = &lhs;
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
    int restrt, iter;
    double resid;
    // int maxit, i;
    // restrt = SUPERLU_MIN(n / 3 + 1, 50);
    // maxit = 1000;
    // iter = maxit;
    // resid = 1e-8;

    //allocate memory for the solution
    restrt = SUPERLU_MIN(n / 3 + 1, mRestart);
    double* x;
    if (!(x = doubleMalloc(n))) ABORT("Malloc fails for x[].");

    /* Initial guess */
    for (int i=0; i<m; ++i)
      x[i] = 0.0;

    iter  = mMaxIteration;
    resid = mTolerance;

    /* Call GMRES */
    dfgmr(n, dmatvec_mult, dpsolve, b, x, resid, restrt, &iter, NULL);

    /* Scale the solution back if equilibration was performed. */
    if (*equed == 'C' || *equed == 'B')
    {
#pragma omp parallel for
      for (int i=0; i<m; ++i)
        x[i] *= C[i];
    }

    //resubstitution of results

    /* Output the result. */
    sp_dgemv("N", -1.0, &lhs, x, 1, 1.0, b, 1);

#pragma omp parallel for
    for (int i=0; i<m; ++i)
      rX[i] = x[i];

    if(options.PrintStat)
    {
      StatPrint(&stat);
    }

    // Deallocate memory used
    StatFree(&stat);

    // Free superlu allocated space
    SUPERLU_FREE (etree);
    SUPERLU_FREE (perm_r);
    SUPERLU_FREE (perm_c);
    SUPERLU_FREE (R);
    SUPERLU_FREE (C);

    // Destroy the matrices and vectors that are not used anymore
    Destroy_SuperMatrix_Store(&lhs);
    Destroy_SuperMatrix_Store(&rhs);
    Destroy_SuperMatrix_Store(&X);
    Destroy_SuperNode_Matrix(&L);
    Destroy_CompCol_Matrix(&U);

    SUPERLU_FREE(b);
    SUPERLU_FREE(x);

    // SUPERLU_FREE(asubt);
    // SUPERLU_FREE(xat);
    // SUPERLU_FREE(at);

    delete [] index1_vector;
    delete [] index2_vector;

    // double stop_solve = OpenMPUtils::GetCurrentTime();
    // if(mEchoLevel > 2)
    //   KRATOS_INFO("superlu_solve_time") << stop_solve - start_solve << std::endl;

    return true;

    //KRATOS_CATCH("")
  }

  ///@}
  ///@name Access
  ///@{
  ///@}
  ///@name Inquiry
  ///@{
  ///@}
  ///@name Input and output
  ///@{

  /// Turn back information as a string.
  std::string Info() const override
  {
    return "SuperLU iterative solver";
  }

  /// Print information about this object.
  void  PrintInfo(std::ostream& rOStream) const override
  {
    rOStream << "SuperLU iterative solver";
  }

  /// Print object's data.
  void  PrintData(std::ostream& rOStream) const override
  {
  }

  ///@}
  ///@name Friends
  ///@{
  ///@}

 protected:

  ///@name Protected static Member Variables
  ///@{
  ///@}
  ///@name Protected member Variables
  ///@{

  double mTolerance;
  double mDropTolerance;
  double mFillTolerance;
  double mFillFactor;

  int mEchoLevel;
  int mMaxIteration;
  int mRestart;

  ///@}
  ///@name Protected Operators
  ///@{
  ///@}
  ///@name Protected Operations
  ///@{
  ///@}
  ///@name Protected  Access
  ///@{
  ///@}
  ///@name Protected Inquiry
  ///@{
  ///@}
  ///@name Protected LifeCycle
  ///@{
  ///@}

 private:
  ///@name Static Member Variables
  ///@{
  ///@}
  ///@name Member Variables
  ///@{
  ///@}
  ///@name Private Operators
  ///@{
  ///@}
  ///@name Private Operations
  ///@{
  ///@}
  ///@name Private  Access
  ///@{
  ///@}
  ///@name Private Inquiry
  ///@{
  ///@}
  ///@name Un accessible methods
  ///@{
  ///@}

}; // Class SuperLUIterativeGMRESSolver

///@}

///@name Input and output
///@{

/// input stream function
template<class TSparseSpaceType, class TDenseSpaceType,class TReordererType>
inline std::istream& operator >> (std::istream& rIStream, SuperLUIterativeGMRESSolver< TSparseSpaceType,
                                  TDenseSpaceType, TReordererType>& rThis)
{
  return rIStream;
}

/// output stream function
template<class TSparseSpaceType, class TDenseSpaceType, class TReordererType>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const SuperLUIterativeGMRESSolver<TSparseSpaceType,
                                  TDenseSpaceType, TReordererType>& rThis)
{
  rThis.PrintInfo(rOStream);
  rOStream << std::endl;
  rThis.PrintData(rOStream);

  return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_SUPERLU_ITERATIVE_GMRES_SOLVER_H_INCLUDED  defined
