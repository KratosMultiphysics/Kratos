//------------------------------------------------------------------
//           ___      _                                            .
//   KRATOS / __| ___| |_ _ ___ _ _ ___                            .
//          \__ \/ _ \ \ V / -_) '_|_-<                            .
//          |___/\___/_|\_/\___|_| /__/ APPLICATION                .
//                                                                 .
//   License:(BSD)	  SolversApplication/license.txt           .
//   Main authors:        Josep Maria Carbonell                    .
//                        ..                                       .
//------------------------------------------------------------------
//
//   Project Name:        KratosSolversApplication $
//   Created by:          $Author:     JMCarbonell $
//   Last modified by:    $Co-Author:              $
//   Date:                $Date:      January 2019 $
//   Revision:            $Revision:           0.0 $
//
//

#if !defined(KRATOS_SOLVERS_APPLICATION_H_INCLUDED )
#define  KRATOS_SOLVERS_APPLICATION_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/kratos_application.h"
#include "includes/variables.h"
#include "spaces/ublas_space.h"
#include "includes/standard_linear_solver_factory.h"

//linear solvers
#ifdef INCLUDE_SUPERLU_MT
  #include "linear_system/linear_solvers/superlu_mt_direct_solver.hpp"
#else
  #include "linear_system/linear_solvers/superlu_direct_solver.hpp"
#endif

#ifdef INCLUDE_FEAST
  #include "linear_system/linear_solvers/feast_solver.hpp"
#endif

namespace Kratos {

///@name Kratos Globals
///@{
///@}
///@name Type Definitions
///@{
///@}
///@name  Enum's
///@{
///@}
///@name  Functions
///@{
///@}
///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.
 */
class KratosSolversApplication : public KratosApplication {
 public:
  ///@name Type Definitions
  ///@{

  /// Pointer definition of KratosSolversApplication
  KRATOS_CLASS_POINTER_DEFINITION(KratosSolversApplication);

  typedef Kratos::Vector                                                     DenseVectorType;
  typedef Kratos::Matrix                                                     DenseMatrixType;
  typedef boost::numeric::ublas::vector<double>                             SparseVectorType;
  typedef UblasSpace<double, CompressedMatrix, SparseVectorType>             SparseSpaceType;
  typedef UblasSpace<double, DenseMatrixType, DenseVectorType>                LocalSpaceType;
  typedef LinearSolver<SparseSpaceType, LocalSpaceType>                     LinearSolverType;

  ///@}
  ///@name Life Cycle
  ///@{

  /// Default constructor.
  KratosSolversApplication();

  /// Destructor.
  ~KratosSolversApplication() override {}

  ///@}
  ///@name Operators
  ///@{
  ///@}
  ///@name Operations
  ///@{

  void Register() override;

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
    return "KratosSolversApplication";
  }

  /// Print information about this object.
  void PrintInfo(std::ostream& rOStream) const override
  {
    rOStream << Info();
    PrintData(rOStream);
  }

  ///// Print object's data.
  void PrintData(std::ostream& rOStream) const override
  {
    KRATOS_WATCH("in Solvers application");
    KRATOS_WATCH(KratosComponents<VariableData>::GetComponents().size() );
    rOStream << "Variables:" << std::endl;
    KratosComponents<VariableData>().PrintData(rOStream);
    rOStream << std::endl;
    rOStream << "Elements:" << std::endl;
    KratosComponents<Element>().PrintData(rOStream);
    rOStream << std::endl;
    rOStream << "Conditions:" << std::endl;
    KratosComponents<Condition>().PrintData(rOStream);
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

#ifdef INCLUDE_SUPERLU_MT
  typedef SuperLUmtDirectSolver<SparseSpaceType, LocalSpaceType>   SuperLUmtDirectSolverType;
  const StandardLinearSolverFactory<SparseSpaceType, LocalSpaceType, SuperLUmtDirectSolverType> mSuperLUmtDirectSolverFactory;
#else
  typedef SuperLUDirectSolver<SparseSpaceType, LocalSpaceType>       SuperLUDirectSolverType;
  const StandardLinearSolverFactory<SparseSpaceType, LocalSpaceType, SuperLUDirectSolverType> mSuperLUDirectSolverFactory;
  //typedef SuperLUIterativeSolver<SparseSpaceType, LocalSpaceType> SuperLUIterativeSolverType;
  //const StandardLinearSolverFactory<SparseSpaceType, LocalSpaceType, SuperLUIterativeSolverType> mSuperLUIterativeSolverFactory;
#endif

#ifdef INCLUDE_FEAST
  typedef FEASTEigenValueSolver<SparseSpaceType, LocalSpaceType> FEASTEigenValueSolverType;
  const StandardLinearSolverFactory<SparseSpaceType, LocalSpaceType, FEASTEigenValueSolverType> mFEASTEigenValueSolverFactory;
#endif

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

  /// Assignment operator.
  KratosSolversApplication& operator=(KratosSolversApplication const& rOther);

  /// Copy constructor.
  KratosSolversApplication(KratosSolversApplication const& rOther);


  ///@}

}; // Class KratosSolversApplication

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif // KRATOS_SOLVERS_APPLICATION_H_INCLUDED  defined
