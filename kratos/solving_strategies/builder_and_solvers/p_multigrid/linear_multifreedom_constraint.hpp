//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Máté Kelemen
//

#pragma once

// Project includes
#include "solving_strategies/builder_and_solvers/p_multigrid/multifreedom_constraint.hpp" // MultifreedomConstraint
#include "includes/smart_pointers.h" // KRATOS_CLASS_POINTER_DEFINITION
#include "includes/code_location.h" // KRATOS_CODE_LOCATION

// System includes
#include <limits> // std::numeric_limits


namespace Kratos {


/// @brief Class representing (part of) a linear constraint equation.
/// @details The relation matrix is stored in the GEOMETRIC_STIFFNESS_MATRIX
///          while the constraint gap in the INTERNAL_FORCES_VECTOR variable
///          of the instance's @ref DataValueContainer.
class KRATOS_API(KRATOS_CORE) LinearMultifreedomConstraint final
    : public MultifreedomConstraint
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(LinearMultifreedomConstraint);

    using MultifreedomConstraint::IndexType;

    using MultifreedomConstraint::DofPointerVectorType;

    using MultifreedomConstraint::MatrixType;

    using MultifreedomConstraint::VectorType;

    LinearMultifreedomConstraint() noexcept
        : LinearMultifreedomConstraint(std::numeric_limits<IndexType>::max())
    {}

    /// @copydoc MultifreedomConstraint::MultifreedomConstraint(const IndexType)
    LinearMultifreedomConstraint(const IndexType Id) noexcept
        : LinearMultifreedomConstraint(Id,
                                       DofPointerVectorType {},
                                       std::vector<std::size_t> {},
                                       MatrixType {},
                                       VectorType {})
    {}

    /// @brief Construct a constraint instance will all necessary information.
    /// @param Id Identifier of the constraint instance. Unrelated to the identifier of the constrain equation.
    /// @param rDofs DoFs participating in the constrain equation(s).
    /// @param rConstraintLabels Identifiers of the constraint equations.
    /// @param rRelationMatrix Matrix storing the coefficients of the participating DoFs. Each row of the matrix
    ///                        defines a constraint equation. The coefficients must be in the same order as the
    ///                        DoFs in @p rDofs. The number of rows must match the size of @p rConstraintLabels.
    /// @param rConstraintGaps Constants of the constraint equations. The constraint gap vector's size must match
    ///                        the number of rows in @p rRelationMatrix.
    LinearMultifreedomConstraint(const IndexType Id,
                                 DofPointerVectorType&& rDofs,
                                 const std::vector<std::size_t>& rConstraintLabels,
                                 const MatrixType& rRelationMatrix,
                                 const VectorType& rConstraintGaps);

    void CalculateLocalSystem(MatrixType& rRelationMatrix,
                              VectorType& rConstraintGaps,
                              const ProcessInfo& rProcessInfo) const override;

    /// @name Unsupported Interface
    /// @{

    MasterSlaveConstraint::Pointer
    Create(IndexType,
           DofPointerVectorType&,
           DofPointerVectorType&,
           const MatrixType&,
           const VectorType&) const override
    {KRATOS_ERROR << KRATOS_CODE_LOCATION.CleanFunctionName() << " is not supported.";}

    MasterSlaveConstraint::Pointer
    Create(IndexType,
           NodeType&,
           const VariableType&,
           NodeType&,
           const VariableType&,
           const double,
           const double) const override
    {KRATOS_ERROR << KRATOS_CODE_LOCATION.CleanFunctionName() << " is not supported.";}

    /// @}
}; // class LinearMultifreedomConstraint


} // namespace Kratos