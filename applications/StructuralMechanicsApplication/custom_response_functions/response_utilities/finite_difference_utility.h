// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Martin Fusseder, https://github.com/MFusseder
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/condition.h"
#include "structural_mechanics_application_variables.h"
#include "includes/adjoint_interface.hpp"

// --- STL Includes ---
#include <span> // std::span


namespace Kratos {


template <IAdjoint::ResidualTerm Term, class TEntity>
requires (std::is_same_v<TEntity,Element> || std::is_same_v<TEntity,Condition>)
class AdjointFiniteDifferenceUtility {
public:
    /// @brief Class uniquely defining what to perturb in an element.
    struct Perturbation {
        /// @brief Variable to perturb.
        /// @details @ref IAdjoint::DynamicVariable optionally offers to augment the variable
        ///          with an integer whose interpretation depends on the @ref mContext "specified context".
        ///          For example, in the context of an element and a @p YOUNG_MODULUS variable, it can refer
        ///          to the Young's modulus of a specific layer within a multi-layered composite material.
        IAdjoint::DynamicVariable mVariable;

        /// @brief Context specifying which container the @ref mVariable "variable" refers to.
        Globals::DataLocation mContext;

        /// @brief Magnitude of the perturbation to apply.
        double mMagnitude;
    }; // struct Perturbation

    using Value = std::conditional_t<
        Term == IAdjoint::ResidualTerm::Load,
        Vector,
        Matrix>;

    void FiniteDifferenceDerivative(
        const TEntity& rEntity,
        const Vector& rValue,
        std::span<const Perturbation> Perturbations,
        Matrix& rOutput,
        std::size_t iBuffer,
        const ProcessInfo& rProcessInfo) const
    requires (Term != IAdjoint::ResidualTerm::Load);

    void FiniteDifferenceDerivative(
        const TEntity& rEntity,
        std::span<const Perturbation> Perturbations,
        Matrix& rOutput,
        std::size_t iBuffer,
        const ProcessInfo& rProcessInfo) const
    requires (Term == IAdjoint::ResidualTerm::Load);

private:
    void FiniteDifferenceDerivative(
        const TEntity& rEntity,
        const Vector* pValue,
        std::span<const Perturbation> Perturbations,
        Matrix& rOutput,
        std::size_t iBuffer,
        const ProcessInfo& rProcessInfo) const;

    static void EvaluateTerm(
        TEntity& rEntity,
        const ProcessInfo& rProcessInfo,
        Value& rOutput);

    template <Globals::DataLocation TContext>
    void Perturb(
        TEntity& rEntity,
        const IAdjoint::DynamicVariable& rVariable,
        double PerturbationMagnitude,
        std::size_t iBuffer) const;

    void ComputeFiniteDifferences(
        Value& rOutput,
        const Value& rReferenceState,
        const Value& rPerturbedState,
        double PerturbationMagnitude,
        std::size_t OutputBegin) const;

    static void ResizeOutput(
        Matrix& rOutput,
        std::span<const std::size_t> SingleVariateDerivativeShape,
        std::size_t VariableCount);
}; // class AdjointFiniteDifferenceUtility


/** \brief FiniteDifferenceUtility
 *
 * This class calculates the derivatives of different element quantities (e.g. RHS, LHS, mass-matrix, ...)
 * with respect to a design variable (e.g. nodal-coordinate, property).
 */


class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) FiniteDifferenceUtility
{
public:

    typedef Variable<double> array_1d_component_type;
    typedef std::size_t IndexType;
    typedef std::size_t SizeType;

    static void CalculateRightHandSideDerivative(Element& rElement,
                                                const Vector& rRHS,
                                                const Variable<double>& rDesignVariable,
                                                const double& rPertubationSize,
                                                Matrix& rOutput,
                                                const ProcessInfo& rCurrentProcessInfo);

    template <typename TElementType>
    static void CalculateRightHandSideDerivative(TElementType& rElement,
                                                const Vector& rRHS,
                                                const array_1d_component_type& rDesignVariable,
                                                Node& rNode,
                                                const double& rPertubationSize,
                                                Vector& rOutput,
                                                const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY;

        if( rDesignVariable == SHAPE_SENSITIVITY_X || rDesignVariable == SHAPE_SENSITIVITY_Y || rDesignVariable == SHAPE_SENSITIVITY_Z )
        {
            const IndexType coord_dir =
                FiniteDifferenceUtility::GetCoordinateDirection(rDesignVariable);

            // define working variables
            Vector RHS_perturbed;

            if (rOutput.size() != rRHS.size())
                rOutput.resize(rRHS.size(), false);

            // perturb the design variable
            rNode.GetInitialPosition()[coord_dir] += rPertubationSize;
            rNode.Coordinates()[coord_dir] += rPertubationSize;

            // compute LHS after perturbation
            rElement.CalculateRightHandSide(RHS_perturbed, rCurrentProcessInfo);

            // compute derivative of RHS w.r.t. design variable with finite differences
            noalias(rOutput) = (RHS_perturbed - rRHS) / rPertubationSize;

            // unperturb the design variable
            rNode.GetInitialPosition()[coord_dir] -= rPertubationSize;
            rNode.Coordinates()[coord_dir] -= rPertubationSize;
        }
        else if( rDesignVariable == TEMPERATURE )
        {
            // define working variables
            Vector RHS_perturbed;

            rOutput.resize(rRHS.size(), false);

            // perturb the design variable
            rNode.FastGetSolutionStepValue(rDesignVariable) += rPertubationSize;

            // compute LHS after perturbation
            rElement.CalculateRightHandSide(RHS_perturbed, rCurrentProcessInfo);

            // compute derivative of RHS w.r.t. design variable with finite differences
            noalias(rOutput) = (RHS_perturbed - rRHS) / rPertubationSize;

            // unperturb the design variable
            rNode.FastGetSolutionStepValue(rDesignVariable) -= rPertubationSize;

        }
        else
        {
            KRATOS_WARNING("FiniteDifferenceUtility") << "Unsupported nodal design variable: " << rDesignVariable << ". " << std::endl
            << "Supported variables are SHAPE_SENSITIVITY_X, SHAPE_SENSITIVITY_Y, SHAPE_SENSITIVITY_Z, TEMPERATURE." << std::endl;
            if ( (rOutput.size() != 0) )
                rOutput.resize(0,false);
        }

        KRATOS_CATCH("");
    }

    static void CalculateLeftHandSideDerivative(Element& rElement,
                                                const Matrix& rLHS,
                                                const array_1d_component_type& rDesignVariable,
                                                Node& rNode,
                                                const double& rPertubationSize,
                                                Matrix& rOutput,
                                                const ProcessInfo& rCurrentProcessInfo);

    static void CalculateMassMatrixDerivative(Element& rElement,
                                                const Matrix& rMassMatrix,
                                                const array_1d_component_type& rDesignVariable,
                                                Node& rNode,
                                                const double& rPertubationSize,
                                                Matrix& rOutput,
                                                const ProcessInfo& rCurrentProcessInfo);

private:

    static std::size_t GetCoordinateDirection(const array_1d_component_type& rDesignVariable);

}; // class FiniteDifferenceUtility.



}  // namespace Kratos.


