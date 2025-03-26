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
#include "includes/master_slave_constraint.h" // MasterSlaveConstraint
#include "includes/code_location.h" // KRATOS_CODE_LOCATION

// System includes
#include <vector> // std::vector
#include <limits> // std::numeric_limits


namespace Kratos {


class KRATOS_API(KRATOS_CORE) MultifreedomConstraint
    : public MasterSlaveConstraint
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(MultifreedomConstraint);

    using MasterSlaveConstraint::IndexType;

    using MasterSlaveConstraint::VariableType;

    using MasterSlaveConstraint::DofPointerVectorType;

    using MasterSlaveConstraint::MatrixType;

    using MasterSlaveConstraint::VectorType;

    using DofArray = std::vector<Dof<double>*>;

    MultifreedomConstraint() noexcept
        : MultifreedomConstraint(std::numeric_limits<IndexType>::max())
    {}

    explicit MultifreedomConstraint(const IndexType Id) noexcept
        : MultifreedomConstraint(Id, DofPointerVectorType {}, std::vector<std::size_t> {})
    {}

    MultifreedomConstraint(const IndexType Id,
                           DofArray&& rDofs,
                           const std::vector<std::size_t>& rConstraintLabels);

    void GetDofList(DofPointerVectorType& rSlaveDofs,
                    DofPointerVectorType& rMasterDofs,
                    const ProcessInfo&) const override;

    void SetDofList(const DofPointerVectorType& rSlaveDofs,
                    const DofPointerVectorType& rMastertDofs,
                    const ProcessInfo&) override;

    void EquationIdVector(EquationIdVectorType&,
                          EquationIdVectorType&,
                          const ProcessInfo& rProcessInfo) const override;

    void ResetSlaveDofs(const ProcessInfo&) override {}

    void Apply(const ProcessInfo&) override {}

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

protected:
    const DofArray& GetDofs() const noexcept
    {
        return mDofs;
    }

private:
    friend class Serializer;

    void save(Serializer& rSerializer) const override;

    void load(Serializer& rDeserializer) override;

    DofArray mDofs;
}; // class MultifreedomConstraint


} // namespace Kratos
