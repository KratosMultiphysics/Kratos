//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//   License:        BSD License
//   Kratos default license: kratos/license.txt

#if !defined(KRATOS_EIGENSOLVER_NITSCHE_DYNAMIC_SCHEME )
#define  KRATOS_EIGENSOLVER_NITSCHE_DYNAMIC_SCHEME

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/condition.h"
#include "includes/process_info.h"
#include "includes/ublas_interface.h"
#include "solving_strategies/schemes/scheme.h"

// Application includes
#include "iga_application_variables.h"

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
class EigensolverNitscheDynamicScheme : public Scheme<TSparseSpace,TDenseSpace>
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION( EigensolverNitscheDynamicScheme );

    typedef Scheme<TSparseSpace,TDenseSpace> BaseType;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    EigensolverNitscheDynamicScheme() : Scheme<TSparseSpace,TDenseSpace>() {}

    /// Destructor.
    ~EigensolverNitscheDynamicScheme() override {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    void CalculateSystemContributions(
        Element& rCurrentElement,
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        const ProcessInfo& CurrentProcessInfo
    ) override
    {
        KRATOS_TRY

        if (CurrentProcessInfo[BUILD_LEVEL] == 1)
        { 
            rCurrentElement.CalculateLocalSystem(LHS_Contribution,RHS_Contribution,CurrentProcessInfo);
        }

        rCurrentElement.EquationIdVector(EquationId,CurrentProcessInfo);

        KRATOS_CATCH("")
    }

    void CalculateSystemContributions(
        Condition& rCurrentCondition,
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemVectorType& RHS_Contribution,
        Condition::EquationIdVectorType& EquationId,
        const ProcessInfo& CurrentProcessInfo) override
    {
        KRATOS_TRY

        if (CurrentProcessInfo[BUILD_LEVEL] == 2) 
        {
            rCurrentCondition.CalculateLocalSystem(LHS_Contribution,RHS_Contribution,CurrentProcessInfo);
        }

        rCurrentCondition.EquationIdVector(EquationId,CurrentProcessInfo);

        KRATOS_CATCH("")
    }

    ///@}

}; /* Class Scheme */

///@}

///@name Type Definitions
///@{

///@}

}  /* namespace Kratos.*/

#endif /* KRATOS_EIGENSOLVER_NITSCHE_DYNAMIC_SCHEME  defined */

