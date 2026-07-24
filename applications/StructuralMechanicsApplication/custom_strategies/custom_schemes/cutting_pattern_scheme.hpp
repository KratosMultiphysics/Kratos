// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Ricky Aristio
//

#pragma once

// System includes

// External includes

// Project includes
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
#include "custom_elements/membrane_elements/membrane_cutting_pattern_element.hpp"

namespace Kratos
{
///@addtogroup StructuralMechanicsApplication
///@{

///@name Kratos Classes
///@{

/**
 * @class CuttingPatternScheme
 * @brief Redirects the standard Newton-Raphson element assembly to
 * MembraneCuttingPatternElement::OptimizationLeastSquare for cutting-pattern elements, instead of
 * CalculateLocalSystem
 */
template<class TSparseSpace, class TDenseSpace>
class CuttingPatternScheme
    : public ResidualBasedIncrementalUpdateStaticScheme<TSparseSpace, TDenseSpace>
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(CuttingPatternScheme);

    using BaseType = ResidualBasedIncrementalUpdateStaticScheme<TSparseSpace, TDenseSpace>;

    using LocalSystemMatrixType = typename BaseType::LocalSystemMatrixType;
    using LocalSystemVectorType = typename BaseType::LocalSystemVectorType;
    using EquationIdVectorType = typename BaseType::EquationIdVectorType;

    ///@}
    ///@name Life Cycle
    ///@{

    CuttingPatternScheme() : BaseType() {}

    ~CuttingPatternScheme() override = default;

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Switches cutting-pattern elements between the least-squares stress-flattening
     * system (OptimizationLeastSquare) and the standard elastic equilibrium system (Relaxation).
     */
    void SetUseRelaxation(bool UseRelaxation)
    {
        mUseRelaxation = UseRelaxation;
    }

    bool GetUseRelaxation() const
    {
        return mUseRelaxation;
    }

    void CalculateSystemContributions(
        Element& rCurrentElement,
        LocalSystemMatrixType& rLHS_Contribution,
        LocalSystemVectorType& rRHS_Contribution,
        EquationIdVectorType& rEquationId,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        if (auto* p_cutting_pattern_element = dynamic_cast<MembraneCuttingPatternElement*>(&rCurrentElement)) {
            if (mUseRelaxation) {
                p_cutting_pattern_element->Relaxation(rLHS_Contribution, rRHS_Contribution, rCurrentProcessInfo);
            } else {
                double response = 0.0;
                p_cutting_pattern_element->OptimizationLeastSquare(rLHS_Contribution, rRHS_Contribution, response, rCurrentProcessInfo);
            }
            rCurrentElement.EquationIdVector(rEquationId, rCurrentProcessInfo);
        } else {
            BaseType::CalculateSystemContributions(rCurrentElement, rLHS_Contribution, rRHS_Contribution, rEquationId, rCurrentProcessInfo);
        }
    }

    void CalculateLHSContribution(
        Element& rCurrentElement,
        LocalSystemMatrixType& rLHS_Contribution,
        EquationIdVectorType& rEquationId,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        if (auto* p_cutting_pattern_element = dynamic_cast<MembraneCuttingPatternElement*>(&rCurrentElement)) {
            LocalSystemVectorType dummy_rhs;
            if (mUseRelaxation) {
                p_cutting_pattern_element->Relaxation(rLHS_Contribution, dummy_rhs, rCurrentProcessInfo);
            } else {
                double response = 0.0;
                p_cutting_pattern_element->OptimizationLeastSquare(rLHS_Contribution, dummy_rhs, response, rCurrentProcessInfo);
            }
            rCurrentElement.EquationIdVector(rEquationId, rCurrentProcessInfo);
        } else {
            BaseType::CalculateLHSContribution(rCurrentElement, rLHS_Contribution, rEquationId, rCurrentProcessInfo);
        }
    }

    void CalculateRHSContribution(
        Element& rCurrentElement,
        LocalSystemVectorType& rRHS_Contribution,
        EquationIdVectorType& rEquationId,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        if (auto* p_cutting_pattern_element = dynamic_cast<MembraneCuttingPatternElement*>(&rCurrentElement)) {
            LocalSystemMatrixType dummy_lhs;
            if (mUseRelaxation) {
                p_cutting_pattern_element->Relaxation(dummy_lhs, rRHS_Contribution, rCurrentProcessInfo);
            } else {
                double response = 0.0;
                p_cutting_pattern_element->OptimizationLeastSquare(dummy_lhs, rRHS_Contribution, response, rCurrentProcessInfo);
            }
            rCurrentElement.EquationIdVector(rEquationId, rCurrentProcessInfo);
        } else {
            BaseType::CalculateRHSContribution(rCurrentElement, rRHS_Contribution, rEquationId, rCurrentProcessInfo);
        }
    }

    ///@}

private:
    ///@name Member Variables
    ///@{

    bool mUseRelaxation = false;

    ///@}

}; // Class CuttingPatternScheme

///@}

///@} addtogroup block

}  // namespace Kratos.
