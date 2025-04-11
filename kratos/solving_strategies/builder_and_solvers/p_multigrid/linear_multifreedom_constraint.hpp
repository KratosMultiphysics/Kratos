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

    LinearMultifreedomConstraint(const IndexType Id) noexcept
        : LinearMultifreedomConstraint(Id,
                                       DofPointerVectorType {},
                                       std::vector<std::size_t> {},
                                       MatrixType {},
                                       VectorType {})
    {}

    LinearMultifreedomConstraint(const IndexType Id,
                                 DofPointerVectorType&& rDofs,
                                 const std::vector<std::size_t>& rConstraintLabels,
                                 const MatrixType& rRelationMatrix,
                                 const VectorType& rConstraintGaps);

    void CalculateLocalSystem(MatrixType& rRelationMatrix,
                              VectorType& rConstraintGaps,
                              const ProcessInfo& rProcessInfo) const override;

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
}; // class LinearMultifreedomConstraint


} // namespace Kratos