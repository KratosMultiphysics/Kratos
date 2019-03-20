//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:   michael.andre@tum.de $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:           September 2016 $
//   Revision:            $Revision:                  0.0 $
//
//


#if !defined(KRATOS_EIGENSOLVER_SCHEME_H_INCLUDED)
#define  KRATOS_EIGENSOLVER_SCHEME_H_INCLUDED


// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/condition.h"
#include "includes/process_info.h"
#include "includes/ublas_interface.h"
#include "custom_solvers/solution_schemes/solution_scheme.hpp"

// Application includes
#include "solid_mechanics_application_variables.h"

namespace Kratos
{

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

/// An adapter scheme for obtaining mass and stiffness matrices for dynamic eigenvalue problems.
template<class TSparseSpace,
         class TDenseSpace
         >
class EigensolverScheme : public SolutionScheme<TSparseSpace,TDenseSpace>
{
 public:
  ///@name Type Definitions
  ///@{

  KRATOS_CLASS_POINTER_DEFINITION( EigensolverScheme );

  typedef SolutionScheme<TSparseSpace,TDenseSpace>                                    BaseType;
  typedef typename BaseType::SolutionSchemePointerType                         BasePointerType;

  typedef typename BaseType::LocalSystemVectorType                       LocalSystemVectorType;
  typedef typename BaseType::LocalSystemMatrixType                       LocalSystemMatrixType;

  ///@}
  ///@name Life Cycle
  ///@{

  /// Default Constructor.
  EigensolverScheme()
      :BaseType()
  {
  }

  /// Constructor.
  EigensolverScheme(Flags& rOptions)
      :BaseType(rOptions)
  {
  }

  /// Destructor.
  ~EigensolverScheme() override {}

  ///@}
  ///@name Operators
  ///@{

  ///@}
  ///@name Operations
  ///@{

  void CalculateSystemContributions(
      Element::Pointer pCurrentElement,
      LocalSystemMatrixType& rLHS_Contribution,
      LocalSystemVectorType& rRHS_Contribution,
      Element::EquationIdVectorType& rEquationId,
      ProcessInfo& rCurrentProcessInfo
                                    ) override
  {
    KRATOS_TRY

    if (rCurrentProcessInfo[BUILD_LEVEL] == 1)
    { // mass matrix
      pCurrentElement->CalculateMassMatrix(rLHS_Contribution,rCurrentProcessInfo);
      std::size_t LocalSize = rLHS_Contribution.size1();
      if (rRHS_Contribution.size() != LocalSize)
        rRHS_Contribution.resize(LocalSize,false);
      noalias(rRHS_Contribution) = ZeroVector(LocalSize);
    }
    else if (rCurrentProcessInfo[BUILD_LEVEL] == 2) // stiffness matrix
    {
      pCurrentElement->CalculateLocalSystem(rLHS_Contribution,rRHS_Contribution,rCurrentProcessInfo);
    }
    else
    {
      KRATOS_ERROR <<"Invalid BUILD_LEVEL" << std::endl;
    }

    pCurrentElement->EquationIdVector(rEquationId,rCurrentProcessInfo);

    KRATOS_CATCH("")
  }

  void Calculate_LHS_Contribution(
      Element::Pointer pCurrentElement,
      LocalSystemMatrixType& rLHS_Contribution,
      Element::EquationIdVectorType& rEquationId,
      ProcessInfo& rCurrentProcessInfo) override
  {
    KRATOS_TRY

    LocalSystemVectorType RHS_Contribution;
    RHS_Contribution.resize(rLHS_Contribution.size1(), false);
    CalculateSystemContributions(
        pCurrentElement,
        rLHS_Contribution,
        RHS_Contribution,
        rEquationId,
        rCurrentProcessInfo);

    KRATOS_CATCH("")
  }

  void Condition_CalculateSystemContributions(
      Condition::Pointer pCurrentCondition,
      LocalSystemMatrixType& rLHS_Contribution,
      LocalSystemVectorType& rRHS_Contribution,
      Condition::EquationIdVectorType& rEquationId,
      ProcessInfo& rCurrentProcessInfo) override
  {
    KRATOS_TRY

    if (rCurrentProcessInfo[BUILD_LEVEL] == 1)
    { // mass matrix
      pCurrentCondition->CalculateMassMatrix(rLHS_Contribution,rCurrentProcessInfo);
      std::size_t LocalSize = rLHS_Contribution.size1();
      if (rRHS_Contribution.size() != LocalSize)
      {
        rRHS_Contribution.resize(LocalSize,false);
      }
      noalias(rRHS_Contribution) = ZeroVector(LocalSize);
    }
    else if (rCurrentProcessInfo[BUILD_LEVEL] == 2) // stiffness matrix
    {
      pCurrentCondition->CalculateLocalSystem(rLHS_Contribution,rRHS_Contribution,rCurrentProcessInfo);
    }
    else
    {
      KRATOS_ERROR <<"Invalid BUILD_LEVEL" << std::endl;
    }

    pCurrentCondition->EquationIdVector(rEquationId,rCurrentProcessInfo);

    KRATOS_CATCH("")
  }

  void Condition_Calculate_LHS_Contribution(
      Condition::Pointer pCurrentCondition,
      LocalSystemMatrixType& rLHS_Contribution,
      Condition::EquationIdVectorType& rEquationId,
      ProcessInfo& rCurrentProcessInfo) override
  {
    KRATOS_TRY

    LocalSystemVectorType RHS_Contribution;
    RHS_Contribution.resize(rLHS_Contribution.size1(), false);
    Condition_CalculateSystemContributions(
        pCurrentCondition,
        rLHS_Contribution,
        RHS_Contribution,
        rEquationId,
        rCurrentProcessInfo);

    KRATOS_CATCH("")
  }

  ///@}
  ///@name Access
  ///@{

  ///@}
  ///@name Inquiry
  ///@{

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

}; // Class Eigen Solver Scheme

///@}

///@name Type Definitions
///@{

///@}

}  // namespace Kratos.

#endif // KRATOS_EIGENSOLVER_SCHEME_H_INCLUDED  defined

