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

#if !defined(KRATOS_MATRIX_BUILDER_SCHEME )
#define  KRATOS_MATRIX_BUILDER_SCHEME


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
// #include "structural_mechanics_application_variables.h"
#include "mor_application_variables.h"

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
class MatrixBuilderScheme : public Scheme<TSparseSpace,TDenseSpace>
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION( MatrixBuilderScheme );

    typedef Scheme<TSparseSpace,TDenseSpace> BaseType;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    MatrixBuilderScheme() : Scheme<TSparseSpace,TDenseSpace>() {}

    /// Destructor.
    ~MatrixBuilderScheme() override {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    void CalculateSystemContributions(
        Element::Pointer pCurrentElement,
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        ProcessInfo& CurrentProcessInfo
    ) override
    {
        KRATOS_TRY

        const int build_level = CurrentProcessInfo[BUILD_LEVEL];
        pCurrentElement->EquationIdVector(EquationId,CurrentProcessInfo);
        
        if ((0 < build_level) && (build_level < 100)) //stiffness matrix
        {
            pCurrentElement->CalculateLocalSystem(LHS_Contribution,RHS_Contribution,CurrentProcessInfo);
        }
        else if ((100 < build_level) && (build_level < 200)) //damping matrix
        {
            pCurrentElement->CalculateDampingMatrix(LHS_Contribution,CurrentProcessInfo);
            std::size_t LocalSize = LHS_Contribution.size1();
            if (RHS_Contribution.size() != LocalSize)
                RHS_Contribution.resize(LocalSize,false);
            noalias(RHS_Contribution) = ZeroVector(LocalSize);
        }
        else if ((200 < build_level) && (build_level < 300)) //mass matrix
        {
            pCurrentElement->CalculateMassMatrix(LHS_Contribution,CurrentProcessInfo);
            std::size_t LocalSize = LHS_Contribution.size1();
            if (RHS_Contribution.size() != LocalSize)
                RHS_Contribution.resize(LocalSize,false);
            noalias(RHS_Contribution) = ZeroVector(LocalSize);
        }
        else if ((300 < build_level) && (build_level < 400)) //various operations (output etc.)
        {
            const std::size_t local_size = EquationId.size();
            if ((LHS_Contribution.size1() != local_size) || (LHS_Contribution.size2() != local_size))
                LHS_Contribution.resize(local_size, local_size, false);
            noalias(LHS_Contribution) = ZeroMatrix(local_size, local_size);
            if (RHS_Contribution.size() != local_size)
                RHS_Contribution.resize(local_size,false);
            noalias(RHS_Contribution) = ZeroVector(local_size);
        }
        else
        {
            KRATOS_ERROR <<"Invalid BUILD_LEVEL" << std::endl;
        }


        KRATOS_CATCH("")
    }

    void Calculate_LHS_Contribution(
        Element::Pointer pCurrentElement,
        LocalSystemMatrixType& LHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        ProcessInfo& CurrentProcessInfo) override
    {
        KRATOS_TRY

        LocalSystemVectorType RHS_Contribution;
        RHS_Contribution.resize(LHS_Contribution.size1(), false);
        CalculateSystemContributions(
                pCurrentElement,
                LHS_Contribution,
                RHS_Contribution,
                EquationId,
                CurrentProcessInfo);

        KRATOS_CATCH("")
    }

    void Calculate_RHS_Contribution(
        Element::Pointer rCurrentElement,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        ProcessInfo& CurrentProcessInfo) override
    {
        KRATOS_TRY

        (rCurrentElement) -> CalculateRightHandSide(RHS_Contribution,CurrentProcessInfo);

        (rCurrentElement) -> EquationIdVector(EquationId,CurrentProcessInfo);

        KRATOS_CATCH("")
    }

    void Condition_CalculateSystemContributions(
        Condition::Pointer pCurrentCondition,
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemVectorType& RHS_Contribution,
        Condition::EquationIdVectorType& EquationId,
        ProcessInfo& CurrentProcessInfo) override
    {
        KRATOS_TRY

        const int build_level = CurrentProcessInfo[BUILD_LEVEL];
        pCurrentCondition->EquationIdVector(EquationId,CurrentProcessInfo);
        
        if ((0 < build_level) && (build_level < 100)) //stiffness matrix
        {
            pCurrentCondition->CalculateLocalSystem(LHS_Contribution,RHS_Contribution,CurrentProcessInfo);
        }
        else if ((100 < build_level) && (build_level < 200)) //damping matrix
        {
            pCurrentCondition->CalculateDampingMatrix(LHS_Contribution,CurrentProcessInfo);
            std::size_t LocalSize = LHS_Contribution.size1();
            if (RHS_Contribution.size() != LocalSize)
                RHS_Contribution.resize(LocalSize,false);
            noalias(RHS_Contribution) = ZeroVector(LocalSize);
        }
        else if ((200 < build_level) && (build_level < 300)) //mass matrix
        {
            pCurrentCondition->CalculateMassMatrix(LHS_Contribution,CurrentProcessInfo);
            std::size_t LocalSize = LHS_Contribution.size1();
            if (RHS_Contribution.size() != LocalSize)
                RHS_Contribution.resize(LocalSize,false);
            noalias(RHS_Contribution) = ZeroVector(LocalSize);
        }
        else if ((300 < build_level) && (build_level < 400)) //various operations for vectors only (output etc.)
        {
            const std::size_t local_size = EquationId.size();
            if ((LHS_Contribution.size1() != local_size) || (LHS_Contribution.size2() != local_size))
                LHS_Contribution.resize(local_size, local_size, false);
            noalias(LHS_Contribution) = ZeroMatrix(local_size, local_size);

            pCurrentCondition->CalculateRightHandSide(RHS_Contribution, CurrentProcessInfo);
        }
        else
        {
            KRATOS_ERROR <<"Invalid BUILD_LEVEL" << std::endl;
        }

        KRATOS_CATCH("")
    }

    void Condition_Calculate_LHS_Contribution(
            Condition::Pointer pCurrentCondition,
            LocalSystemMatrixType& LHS_Contribution,
            Condition::EquationIdVectorType& EquationId,
            ProcessInfo& CurrentProcessInfo) override
    {
        KRATOS_TRY

        LocalSystemVectorType RHS_Contribution;
        RHS_Contribution.resize(LHS_Contribution.size1(), false);
        Condition_CalculateSystemContributions(
                pCurrentCondition,
                LHS_Contribution,
                RHS_Contribution,
                EquationId,
                CurrentProcessInfo);

        KRATOS_CATCH("")
    }

    void Condition_Calculate_RHS_Contribution(
        Condition::Pointer rCurrentCondition,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        ProcessInfo& CurrentProcessInfo) override
    {
        KRATOS_TRY

        (rCurrentCondition) -> CalculateRightHandSide(RHS_Contribution,CurrentProcessInfo);

        (rCurrentCondition) -> EquationIdVector(EquationId,CurrentProcessInfo);

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

}; /* Class Scheme */

///@}

///@name Type Definitions
///@{

///@}

}  /* namespace Kratos.*/

#endif /* KRATOS_MATRIX_BUILDER_SCHEME  defined */

