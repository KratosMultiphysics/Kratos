//
//   Project Name:        KratosSolversApplication $
//   Created by:          $Author:     JMCarbonell $
//   Last modified by:    $Co-Author:              $
//   Date:                $Date:      January 2019 $
//   Revision:            $Revision:           0.0 $
//
//

#if !defined(KRATOS_SUPERLU_DIRECT_SOLVER_H_INCLUDED)
#define  KRATOS_SUPERLU_DIRECT_SOLVER_H_INCLUDED

// System includes

// External includes
#include "includes/ublas_interface.h"

extern "C"{
   #include "slu_ddefs.h"
}

// Project includes
#include "linear_solvers/direct_solver.h"
#include "utilities/openmp_utils.h"


namespace Kratos
{

///@name Kratos Classes
///@{

template< class TSparseSpaceType, class TDenseSpaceType, class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType> >
class SuperLUDirectSolver : public DirectSolver< TSparseSpaceType, TDenseSpaceType, TReordererType>
{
 public:
  ///@name Kratos Globals
  ///@{

  ///@}
  ///@name Type Definitions
  ///@{

  /// Pointer definition for the Solver
  KRATOS_CLASS_POINTER_DEFINITION( SuperLUDirectSolver );

  typedef DirectSolver<TSparseSpaceType, TDenseSpaceType, TReordererType> BaseType;

  typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

  typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

  typedef typename TSparseSpaceType::VectorType VectorType;

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
  SuperLUDirectSolver() {}

  /// Default constructor.
  SuperLUDirectSolver(Parameters settings) : BaseType(settings)
  {
    KRATOS_TRY

    Parameters default_parameters( R"(
    {
        "solver_type"         : "superlu_direct",
        "tolerance"           : 1e-7,
	"max_iteration"       : 5000,
        "scaling"             : false,
        "verbosity"           : 1
    }  )" );

    //now validate agains defaults -- this also ensures no type mismatch
    settings.ValidateAndAssignDefaults(default_parameters);

    mEchoLevel = settings["verbosity"].GetInt();

    KRATOS_CATCH("")
  }

  /// Copy constructor.
  SuperLUDirectSolver(const SuperLUDirectSolver& Other) : BaseType(Other) {}

  /// Destructor.
  ~SuperLUDirectSolver() override {}

  ///@}
  ///@name Operators
  ///@{

  /// Assignment operator.
  SuperLUDirectSolver& operator=(const SuperLUDirectSolver& Other)
  {
    this->mpReorderer = Other.mpReorderer;

    return *this;
  }

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
    // double start_define = OpenMPUtils::GetCurrentTime();

    if(this->IsNotConsistent(rA, rX, rB))
      return false;

    superlu_options_t options;

    /* Set the default input options: (constructor settings)
       options.Fact = DOFACT;
       options.Equil = YES;
       options.ColPerm = COLAMD;
       options.DiagPivotThresh = 1.0;
       options.Trans = NOTRANS;
       options.IterRefine = NOREFINE;
       options.SymmetricMode = NO;
       options.PivotGrowth = NO;
       options.ConditionNumber = NO;
       options.PrintStat = YES;
    */

    set_default_options(&options);
    options.PrintStat = NO;
    options.IterRefine = SLU_DOUBLE;

    //Fill the SuperLU matrices
    SuperMatrix lhs, rhs, L, U;


    //can not avoid a copy as ublas uses unsigned int internally
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


    //create a copy of the rhs vector (it will be overwritten with the solution)
    int size = rB.size();
    VectorType b_rhs(size);
    noalias(b_rhs) = rB;
    // double *b_rhs = new (std::nothrow) double[rB.size()];
    // for( unsigned int i = 0; i < rB.size(); i++ )
    //   b_rhs[i] = rB[i];

    // Creation of a Column Matrix (option0)
    //dCreate_CompCol_Matrix (&lhs, rA.size1(), rA.size2(), rA.nnz(), values_vector, index2_vector, index1_vector, SLU_NR, SLU_D, SLU_GE );

    // Creation of a Column Matrix (option1)
    //dCreate_CompCol_Matrix (&lhs, rA.size1(), rA.size2(), rA.nnz(), rA.value_data().begin(), index2_vector, index1_vector, SLU_NR, SLU_D, SLU_GE );

    // Creation of a Row Matrix (option2)
    dCreate_CompRow_Matrix (&lhs, rA.size1(), rA.size2(), rA.nnz(), rA.value_data().begin(), index2_vector, index1_vector, SLU_NR, SLU_D, SLU_GE);

    dCreate_Dense_Matrix (&rhs, size, 1, &b_rhs[0], size, SLU_DN, SLU_D, SLU_GE);

    //allocate memory for permutation arrays
    int* perm_c;
    int* perm_r;
    if ( !(perm_c = intMalloc(rA.size1())) ) ABORT("Malloc fails for perm_c[].");
    if ( !(perm_r = intMalloc(rA.size2())) ) ABORT("Malloc fails for perm_r[].");

    // double stop_define = OpenMPUtils::GetCurrentTime();
    // if(mEchoLevel > 2)
    //   KRATOS_INFO("superlu_define_time") << stop_define - start_define << std::endl;


    //initialize container for statistical data
    SuperLUStat_t stat;
    int info;

    // double start_solve = OpenMPUtils::GetCurrentTime();

    // Resolution of the linear system:
    StatInit(&stat);

    dgssv(&options, &lhs, perm_c, perm_r, &L, &U, &rhs, &stat, &info);

    if(options.PrintStat)
    {
      StatPrint(&stat);
    }

    //resubstitution of results
#pragma omp parallel for
    for(int i=0; i<size; ++i)
      rX[i] = b_rhs[i]; // rhs(i,0);

    // Deallocate memory used
    StatFree(&stat);

    // Free superlu allocated space
    SUPERLU_FREE (perm_r);
    SUPERLU_FREE (perm_c);

    // Destroy the matrices and vectors that are not used anymore
    Destroy_SuperMatrix_Store(&lhs);
    Destroy_SuperMatrix_Store(&rhs);
    Destroy_SuperNode_Matrix(&L);
    Destroy_CompCol_Matrix(&U);

    delete [] index1_vector;
    delete [] index2_vector;

    // double stop_solve = OpenMPUtils::GetCurrentTime();
    // if(mEchoLevel > 2)
    //   KRATOS_INFO("superlu_solve_time") << stop_solve - start_solve << std::endl;

    return true;
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
    return "SuperLU direct solver";
  }

  /// Print information about this object.
  void PrintInfo(std::ostream& rOStream) const override
  {
    rOStream << "SuperLU direct solver";
  }

  /// Print object's data.
  void PrintData(std::ostream& rOStream) const override
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

  int mEchoLevel;

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


}; // Class SuperLUDirectSolver

///@}

///@name Input and output
///@{

/// input stream function
template<class TSparseSpaceType, class TDenseSpaceType,class TReordererType>
inline std::istream& operator >> (std::istream& rIStream, SuperLUDirectSolver< TSparseSpaceType,
                                  TDenseSpaceType, TReordererType>& rThis)
{
  return rIStream;
}

/// output stream function
template<class TSparseSpaceType, class TDenseSpaceType, class TReordererType>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const SuperLUDirectSolver<TSparseSpaceType,
                                  TDenseSpaceType, TReordererType>& rThis)
{
  rThis.PrintInfo(rOStream);
  rOStream << std::endl;
  rThis.PrintData(rOStream);

  return rOStream;
}
///@}

}  // namespace Kratos.

#endif // KRATOS_SUPERLU_DIRECT_SOLVER_H_INCLUDED  defined
