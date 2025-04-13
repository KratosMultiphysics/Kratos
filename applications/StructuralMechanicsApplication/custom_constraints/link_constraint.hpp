//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Máté Kelemen
//

// Project includes
#include "solving_strategies/builder_and_solvers/p_multigrid/multifreedom_constraint.hpp" // MultifreedomConstraint
#include "includes/code_location.h" // KRATOS_CODE_LOCATION

// STL includes
#include <memory> // std::unique_ptr


namespace Kratos {


/** @brief A constraint enforcing the distance between two @ref Node "nodes" to remain constant.
 *  @details Let @f$ n^0_i @f$ and @f$ n^1_i @f$ denote the coordinates of the initial positions of @ref Node "node"
 *           @f$ n^0 @f$ and @f$ n^1 @f$ that get displaced by @f$ u^0_i @f$ and @f$ u^1_i @f$ respectively.
 *           @p LinkConstraint then represents the following constraint equation:
 *           @f[
 *              \sum_i \left( (n^0_i + u^0_i) - (n^1_i + u^1_i) \right)^2
 *              -
 *              \sum_i \left( n^0_i - n^1_i \right)^2
 *              = 0
 *           @f]
 *
 *           The constraint equation above is linearized in each nonlinear iteration @f$ k @f$
 *           to the following form:
 *           @f[
 *              \sum_i \left( (2 n^0_i + u^0_{i,k-1}) -  (2 n^1_i + u^1_{i,k-1}) \right)
 *              \left( u^0_{i,k} - u^1_{i,k} \right)
 *              = 0
 *           @f]
 *
 *  @note @p LinkConstraint is nonlinear and as such, should not be imposed via master-slave elimination (the slave
 *        DoF's coefficient may vanish). Consider using the @ref AugmentedLagrangeConstraintAssembler "penalty method",
 *        the method of Lagrange multipliers, or @ref AugmentedLagrangeConstraintAssembler "augmented lagrange" multipliers.
 */
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) LinkConstraint final
    : public MultifreedomConstraint
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(LinkConstraint);

    LinkConstraint() noexcept = default;

    LinkConstraint(const IndexType Id,
                   Node& rFirst,
                   Node& rSecond,
                   const std::size_t Dimensions,
                   bool IsMeshMoved);

    LinkConstraint(LinkConstraint&& rRhs) noexcept;

    LinkConstraint(const LinkConstraint& rRhs);

    LinkConstraint& operator=(LinkConstraint&& rRhs) noexcept;

    LinkConstraint& operator=(const LinkConstraint& rRhs);

    /// @internal
    ~LinkConstraint();

    MasterSlaveConstraint::Pointer Clone(IndexType NewId) const override;

    void Initialize(const ProcessInfo& rProcessInfo) override;

    void InitializeSolutionStep(const ProcessInfo& rProcessInfo) override;

    void InitializeNonLinearIteration(const ProcessInfo& rProcessInfo) override;

    void FinalizeNonLinearIteration(const ProcessInfo& rProcessInfo) override;

    void CalculateLocalSystem(MatrixType& rRelationMatrix,
                              VectorType& rConstraintGaps,
                              const ProcessInfo& rProcessInfo) const override;

    int Check(const ProcessInfo& rProcessInfo) const override;

    /// @name Unsupported Interface
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
