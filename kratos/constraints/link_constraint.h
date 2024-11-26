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
#include <array> // std::array
#include <memory> // std::unique_ptr


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

    LinkConstraint(LinkConstraint&& rRhs) noexcept;

    LinkConstraint(const LinkConstraint& rRhs);

    LinkConstraint& operator=(LinkConstraint&& rRhs) noexcept;

    LinkConstraint& operator=(const LinkConstraint& rRhs);

    ~LinkConstraint();

    MasterSlaveConstraint::Pointer Clone(IndexType NewId) const override;

    void Initialize(const ProcessInfo& rProcessInfo) override;

    void InitializeSolutionStep(const ProcessInfo& rProcessInfo) override;

    void InitializeNonLinearIteration(const ProcessInfo& rProcessInfo) override;

    void FinalizeNonLinearIteration(const ProcessInfo& rProcessInfo) override;

    void ResetSlaveDofs(const ProcessInfo& rProcessInfo) override;

    void Apply(const ProcessInfo& rProcessInfo) override;

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
    {KRATOS_ERROR << KRATOS_CODE_LOCATION.CleanFunctionName() << " is not supported by LinkConstraint.";}

    MasterSlaveConstraint::Pointer Create(IndexType,
                                          NodeType&,
                                          const VariableType&,
                                          NodeType&,
                                          const VariableType&,
                                          const double,
                                          const double) const override
    {KRATOS_ERROR << KRATOS_CODE_LOCATION.CleanFunctionName() << " is not supported by LinkConstraint.";}

    void SetDofList(const DofPointerVectorType&,
                    const DofPointerVectorType&,
                    const ProcessInfo&) override
    {KRATOS_ERROR << KRATOS_CODE_LOCATION.CleanFunctionName() << " is not supported by LinkConstraint.";}

    void SetSlaveDofsVector(const DofPointerVectorType&) override
    {KRATOS_ERROR << KRATOS_CODE_LOCATION.CleanFunctionName() << " is not supported by LinkConstraint.";}

    void SetMasterDofsVector(const DofPointerVectorType&) override
    {KRATOS_ERROR << KRATOS_CODE_LOCATION.CleanFunctionName() << " is not supported by LinkConstraint.";}

    void SetLocalSystem(const MatrixType&,
                        const VectorType&,
                        const ProcessInfo&) override
    {KRATOS_ERROR << KRATOS_CODE_LOCATION.CleanFunctionName() << " is not supported by LinkConstraint.";}

    /// @}

private:
    struct Impl;
    std::unique_ptr<Impl> mpImpl;
}; // class LinkConstraint


} // namespace Kratos
