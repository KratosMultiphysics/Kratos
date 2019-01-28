//
//   Project Name:        KratosSolversApplication $
//   Created by:          $Author:     JMCarbonell $
//   Last modified by:    $Co-Author:              $
//   Date:                $Date:      January 2019 $
//   Revision:            $Revision:           0.0 $
//
//

#if !defined(KRATOS_SUPERLU_MT_DIRECT_SOLVER_H_INCLUDED )
#define  KRATOS_SUPERLU_MT_DIRECT_SOLVER_H_INCLUDED

// System includes

// External includes
#include "includes/ublas_interface.h"

extern "C"{
   #include "slu_mt_ddefs.h"
}

// Project includes
#include "linear_solvers/direct_solver.h"
#include "utilities/openmp_utils.h"


namespace Kratos
{

template< class TSparseSpaceType, class TDenseSpaceType, class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType> >
class SuperLUmtDirectSolver : public DirectSolver< TSparseSpaceType, TDenseSpaceType, TReordererType>
{
 public:
  ///@name Kratos Globals
  ///@{

  ///@}
  ///@name Type Definitions
  ///@{

  /// Pointer definition for the Solver
  KRATOS_CLASS_POINTER_DEFINITION( SuperLUmtDirectSolver );

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
  SuperLUmtDirectSolver() {}

  /// Default constructor.
  SuperLUmtDirectSolver(Parameters settings) : BaseType(settings)
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

    mEchoLevel    = settings["verbosity"].GetInt();

    KRATOS_CATCH("")
  }

  /// Copy constructor.
  SuperLUmtDirectSolver(const SuperLUmtDirectSolver& Other) : BaseType(Other) {}

  /// Destructor.
  ~SuperLUmtDirectSolver() override {}

  ///@}
  ///@name Operators
  ///@{

  /// Assignment operator.
  SuperLUmtDirectSolver& operator=(const SuperLUmtDirectSolver& Other)
  {
    this->mpReorderer = Other.mpReorderer;

    return *this;
  }

  ///@}
  ///@name Operations
  ///@{


  bool Solve(SparseMatrixType& rA, VectorType& rX, VectorType& rB) override
  {
    // double start_define = OpenMPUtils::GetCurrentTime();

    if(this->IsNotConsistent(rA, rX, rB))
      return false;

    //Fill the SuperLUmt matrices
    SuperMatrix lhs, rhs, L, U;

    // number of threads
    int_t nprocs = OpenMPUtils::GetNumThreads();

    //can not avoid a copy as ublas uses unsigned int internally
    int_t *index1_vector = new (std::nothrow) int[rA.index1_data().size()];
    for( int unsigned i = 0; i < rA.index1_data().size(); i++ )
      index1_vector[i] = (int_t)rA.index1_data()[i];

    //can not avoid a copy as ublas uses unsigned int internally
    int_t *index2_vector = new (std::nothrow) int[rA.index2_data().size()];
    for( unsigned int i = 0; i < rA.index2_data().size(); i++ )
      index2_vector[i] = (int_t)rA.index2_data()[i];

    // double *values_vector = new (std::nothrow) double[rA.value_data().size()];
    // for( unsigned int i = 0; i < rA.value_data().size(); i++ )
    //   values_vector[i] = (double)rA.value_data()[i];

    //create a copy of the rhs vector (it will be overwritten with the solution)
    int size = rB.size();
    VectorType b_rhs(size);
    noalias(b_rhs)= rB;

    // Creation of a Column Matrix (option0)
    //dCreate_CompCol_Matrix (&lhs, rA.size1(), rA.size2(), rA.nnz(), values_vector, index2_vector, index1_vector, SLU_NR, SLU_D, SLU_GE);

    // Creation of a Column Matrix (option1)
    // dCreate_CompCol_Matrix (&lhs, rA.size1(), rA.size2(), rA.nnz(), rA.value_data().begin(), index2_vector, index1_vector, SLU_NR, SLU_D, SLU_GE);

    // Creation of a Row Matrix (option2)
    dCreate_CompRow_Matrix (&lhs, rA.size1(), rA.size2(), rA.nnz(), rA.value_data().begin(), index2_vector, index1_vector, SLU_NR, SLU_D, SLU_GE);

    dCreate_Dense_Matrix (&rhs, size, 1, &b_rhs[0], size, SLU_DN, SLU_D, SLU_GE);

    //allocate memory for permutation arrays
    int_t* perm_c;
    int_t* perm_r;
    if ( !(perm_c = intMalloc(rA.size1())) ) SUPERLU_ABORT("Malloc fails for perm_c[].");
    if ( !(perm_r = intMalloc(rA.size2())) ) SUPERLU_ABORT("Malloc fails for perm_r[].");

    // double stop_define = OpenMPUtils::GetCurrentTime();
    // if(mEchoLevel > 2)
    //   KRATOS_INFO("superlu_define_time") << stop_define - start_define << std::endl;

    //initialize container for information
    int_t info;

    //double start_solve = OpenMPUtils::GetCurrentTime();

    /*
     * Get column permutation vector perm_c[], according to permc_spec:
     *   permc_spec = 0: natural ordering
     *   permc_spec = 1: minimum degree ordering on structure of A'*A
     *   permc_spec = 2: minimum degree ordering on structure of A'+A
     *   permc_spec = 3: approximate minimum degree for unsymmetric matrices
     */
    int_t permc_spec = 0;
    get_perm_c(permc_spec, &lhs, perm_c);

    // Resolution of the linear system:
    pdgssv(nprocs, &lhs, perm_c, perm_r, &L, &U, &rhs, &info);

    //resubstitution of results
#pragma omp parallel for
    for(int i=0; i<size; ++i)
      rX[i] = b_rhs[i]; // rhs(i,0);

    // Free superlu allocated space
    SUPERLU_FREE (perm_r);
    SUPERLU_FREE (perm_c);

    // Destroy the matrices and vectors that are not used anymore
    //Destroy_CompCol_Matrix(&lhs);
    Destroy_SuperMatrix_Store(&lhs);
    Destroy_SuperMatrix_Store(&rhs);
    Destroy_SuperNode_SCP(&L);
    Destroy_CompCol_NCP(&U);

    delete [] index1_vector;
    delete [] index2_vector;
    //delete [] values_vector;

    std::cout.clear();

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
    return "SuperLUmtDirectSolver";
  }

  /// Print information about this object.
  void PrintInfo(std::ostream& rOStream) const override
  {
    rOStream << "SuperLUmtDirectSolver";
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


}; // Class SuperLUmtDirectSolver

///@}

///@name Input and output
///@{

/// input stream function
template<class TSparseSpaceType, class TDenseSpaceType,class TReordererType>
inline std::istream& operator >> (std::istream& rIStream, SuperLUmtDirectSolver< TSparseSpaceType,
                                  TDenseSpaceType, TReordererType>& rThis)
{
  return rIStream;
}

/// output stream function
template<class TSparseSpaceType, class TDenseSpaceType, class TReordererType>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const SuperLUmtDirectSolver<TSparseSpaceType,
                                  TDenseSpaceType, TReordererType>& rThis)
{
  rThis.PrintInfo(rOStream);
  rOStream << std::endl;
  rThis.PrintData(rOStream);

  return rOStream;
}
///@}

}  // namespace Kratos.

#endif // KRATOS_SUPERLU_MT_DIRECT_SOLVER_H_INCLUDED  defined
