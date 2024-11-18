//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Mate Kelemen
//

// Project includes
#include "includes/master_slave_constraint.h" // MasterSlaveConstraint
#include "includes/code_location.h" // KRATOS_CODE_LOCATION

// STL includes
#include <optional> // std::optional
#include <array> // std::array


namespace Kratos {


class LinkConstraint final : public MasterSlaveConstraint
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(LinkConstraint);

    LinkConstraint() noexcept = default;

    LinkConstraint(const IndexType Id,
                   Node& rFirst,
                   Node& rSecond,
                   const unsigned Dimensions);

    MasterSlaveConstraint::Pointer Clone(IndexType NewId) const override;

    void InitializeSolutionStep(const ProcessInfo& rProcessInfo) override;

    void InitializeNonLinearIteration(const ProcessInfo& rProcessInfo) override;

    void GetDofList(DofPointerVectorType& rSlaves,
                    DofPointerVectorType& rMasters,
                    const ProcessInfo& rProcessInfo) const override;

    void EquationIdVector(EquationIdVectorType& rSlaves,
                          EquationIdVectorType& rMasters,
                          const ProcessInfo& rProcessInfo) const override;

    const DofPointerVectorType& GetSlaveDofsVector() const override;

    const DofPointerVectorType& GetMasterDofsVector() const override;

    void CalculateLocalSystem(MatrixType& rRelationMatrix,
                              VectorType& rConstraintGaps,
                              const ProcessInfo& rProcessInfo) const override;

    int Check(const ProcessInfo& rProcessInfo) const override;

    /// @name Unsupported Virtual Functions
    /// @{

    MasterSlaveConstraint::Pointer Create(IndexType,
                                          DofPointerVectorType&,
                                          DofPointerVectorType&,
                                          const MatrixType&,
                                          const VectorType&) const override
    {KRATOS_ERROR << KRATOS_CODE_LOCATION.GetFunctionName() << " is not supported by LinkConstraint.";}

    MasterSlaveConstraint::Pointer Create(IndexType,
                                          NodeType&,
                                          const VariableType&,
                                          NodeType&,
                                          const VariableType&,
                                          const double,
                                          const double) const override
    {KRATOS_ERROR << KRATOS_CODE_LOCATION.GetFunctionName() << " is not supported by LinkConstraint.";}

    void SetDofList(const DofPointerVectorType&,
                    const DofPointerVectorType&,
                    const ProcessInfo&) override
    {KRATOS_ERROR << KRATOS_CODE_LOCATION.GetFunctionName() << " is not supported by LinkConstraint.";}

    void SetSlaveDofsVector(const DofPointerVectorType&) override
    {KRATOS_ERROR << KRATOS_CODE_LOCATION.GetFunctionName() << " is not supported by LinkConstraint.";}

    void SetMasterDofsVector(const DofPointerVectorType&) override
    {KRATOS_ERROR << KRATOS_CODE_LOCATION.GetFunctionName() << " is not supported by LinkConstraint.";}

    void ResetSlaveDofs(const ProcessInfo& rProcessInfo) override
    {KRATOS_ERROR << KRATOS_CODE_LOCATION.GetFunctionName() << " is not supported by LinkConstraint.";}

    void SetLocalSystem(const MatrixType&,
                        const VectorType&,
                        const ProcessInfo&) override
    {KRATOS_ERROR << KRATOS_CODE_LOCATION.GetFunctionName() << " is not supported by LinkConstraint.";}

    void Apply(const ProcessInfo& rProcessInfo) override
    {KRATOS_ERROR << KRATOS_CODE_LOCATION.GetFunctionName() << " is not supported by LinkConstraint.";}

    /// @}

private:
    unsigned mDimensions;

    /// @details MasterSlaveConstraint::GetSlaveDofsVector and MasterSlaveConstraint::GetMasterDofsVector
    ///          require arrays of mutable Dof pointers, which are only obtainable from mutable nodes,
    ///          so the nodes pointers stored here cannot be immutable. Risky business.
    std::optional<std::array<Node*,2>> mNodePair;

    /// @details Unfortunately, the MasterSlaveConstraint interface demands that
    ///          GetSlaveDofsVector and GetMasterDofsVector returns the array of
    ///          Dofs by reference, even though those must be computed dynamically
    ///          by LinkConstraint. As a result, these vectors must be stored as
    ///          member variables and updated in InitializeNonlinearIteration, instead
    ///          of being computed on the fly.
    std::optional<std::array<DofPointerVectorType,2>> mDofVectors; // {slave, masters}
}; // class LinkConstraint


} // namespace Kratos
