//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//   License:        BSD License
//   Kratos default license: kratos/license.txt
//
//   Project Name:        $StructuralMechanicsApplication $
//   Last modified by:    $Author: michael.andre@tum.de   $
//   Date:                $Date:         September 2016   $
//   Revision:            $Revision:                0.0   $

#if !defined(KRATOS_EIGENSOLVER_DYNAMIC_SCHEME )
#define  KRATOS_EIGENSOLVER_DYNAMIC_SCHEME


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
#include "structural_mechanics_application_variables.h"

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
class EigensolverDynamicScheme : public Scheme<TSparseSpace,TDenseSpace>
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION( EigensolverDynamicScheme );

    typedef Scheme<TSparseSpace,TDenseSpace> BaseType;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    EigensolverDynamicScheme() : Scheme<TSparseSpace,TDenseSpace>() {}

    /// Destructor.
    ~EigensolverDynamicScheme() override {}

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
        { // mass matrix
            rCurrentElement.CalculateMassMatrix(LHS_Contribution,CurrentProcessInfo);
            std::size_t LocalSize = LHS_Contribution.size1();
            if (RHS_Contribution.size() != LocalSize)
                RHS_Contribution.resize(LocalSize,false);
            noalias(RHS_Contribution) = ZeroVector(LocalSize);
        }
        else if (CurrentProcessInfo[BUILD_LEVEL] == 2) // stiffness matrix
        {
            rCurrentElement.CalculateLocalSystem(LHS_Contribution,RHS_Contribution,CurrentProcessInfo);
        }
        else
        {
            KRATOS_ERROR <<"Invalid BUILD_LEVEL" << std::endl;
        }

        rCurrentElement.EquationIdVector(EquationId,CurrentProcessInfo);

        KRATOS_CATCH("")
    }

    void CalculateLHSContribution(
        Element& rCurrentElement,
        LocalSystemMatrixType& LHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        const ProcessInfo& CurrentProcessInfo) override
    {
        KRATOS_TRY

        LocalSystemVectorType RHS_Contribution;
        RHS_Contribution.resize(LHS_Contribution.size1(), false);
        CalculateSystemContributions(
                rCurrentElement,
                LHS_Contribution,
                RHS_Contribution,
                EquationId,
                CurrentProcessInfo);

        KRATOS_CATCH("")
    }

    void CalculateRHSContribution(
        Element& rCurrentElement,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        const ProcessInfo& CurrentProcessInfo) override
    {
        KRATOS_TRY

        rCurrentElement.CalculateRightHandSide(RHS_Contribution,CurrentProcessInfo);

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

        if (CurrentProcessInfo[BUILD_LEVEL] == 1)
        { // mass matrix
            rCurrentCondition.CalculateMassMatrix(LHS_Contribution,CurrentProcessInfo);
            std::size_t LocalSize = LHS_Contribution.size1();
            if (RHS_Contribution.size() != LocalSize)
            {
                RHS_Contribution.resize(LocalSize,false);
            }
            noalias(RHS_Contribution) = ZeroVector(LocalSize);
        }
        else if (CurrentProcessInfo[BUILD_LEVEL] == 2) // stiffness matrix
        {
            rCurrentCondition.CalculateLocalSystem(LHS_Contribution,RHS_Contribution,CurrentProcessInfo);
        }
        else
        {
            KRATOS_ERROR <<"Invalid BUILD_LEVEL" << std::endl;
        }

        rCurrentCondition.EquationIdVector(EquationId,CurrentProcessInfo);

        KRATOS_CATCH("")
    }

    void CalculateLHSContribution(
        Condition& rCurrentCondition,
        LocalSystemMatrixType& LHS_Contribution,
        Condition::EquationIdVectorType& EquationId,
        const ProcessInfo& CurrentProcessInfo) override
    {
        KRATOS_TRY

        LocalSystemVectorType RHS_Contribution;
        RHS_Contribution.resize(LHS_Contribution.size1(), false);
        CalculateSystemContributions(
                rCurrentCondition,
                LHS_Contribution,
                RHS_Contribution,
                EquationId,
                CurrentProcessInfo);

        KRATOS_CATCH("")
    }

    void CalculateRHSContribution(
        Condition& rCurrentCondition,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        const ProcessInfo& CurrentProcessInfo) override
    {
        KRATOS_TRY

        rCurrentCondition.CalculateRightHandSide(RHS_Contribution,CurrentProcessInfo);

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

#endif /* KRATOS_EIGENSOLVER_DYNAMIC_SCHEME  defined */

