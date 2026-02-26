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
#include "includes/variables.h" // GEOMETRIC_STIFFNESS_MATRIX

// System includes
#include <vector> // std::vector
#include <limits> // std::numeric_limits


namespace Kratos {


/** @brief Class representing (part of) a linearizable constraint equation on @ref Dof "DoFs".
 *  @details MultifreedomConstraint offers a slightly different interface than @ref MasterSlaveConstraint,
 *           which it inherits from. The key difference is that @p MultifreedomConstraint does not impose
 *           a partitioning of its DoFs into slaves and masters.
 *
 *           Another important difference is how constraint equations are identified. When a constraint
 *           equations is broken up into several @p MasterSlaveConstraints, instances belonging to the
 *           same constraint equation are the ones that have identical slave DoF IDs.
 *           Since @p MultifreedomConstraint provides no slave DoFs, each instance stores an array of
 *           constraint equation identifiers. This vector is stored in the @ref CONSTRAINT_LABELS variable
 *           of the @ref DataValueContainer belonging to the constraint instance.
 *           @code
 *           MultifreedomConstraint& r_constraint = ...;
 *           auto& r_constraint_ids = r_constraint[CONSTRAINT_LABELS];
 *           @endcode
 *           Instances sharing constraint labels belong to the same constraint equation.
 *
 *  @note Since this class does not partition its DoFs into slaves and masters, any member function returning
 *        data related to slaves/masters will leave slaves empty and only populate masters.
 *        - @ref MultifreedomConstraint::GetDofList
 *        - @ref MultifreedomConstraint::SetDofList
 *        - @ref MultifreedomConstraint::EquationIdVector
 */
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

    /** @brief Construct with topological data.
     *  @param Id Identifier of the constraint instance. Note that, in general, this ID
     *            is not related to the ID of the constraint this instance is part of.
     *  @param rDofs Pointers to the @ref Dof "DoFs" this constraint is defined on.
     *  @param rConstraintLabels Identifiers of the constraint equations this instance is part of.
     *                           Instances sharing constraint labels are assumed to be part of the
     *                           same constraint equations, and are summed up during assembly.
     */
    MultifreedomConstraint(const IndexType Id,
                           DofArray&& rDofs,
                           const std::vector<std::size_t>& rConstraintLabels);

    /// @copydoc MasterSlaveConstraint::GetDofList
    /// @note Slave DoFs are left empty and only master DoFs are populated.
    void GetDofList(DofPointerVectorType& rSlaveDofs,
                    DofPointerVectorType& rMasterDofs,
                    const ProcessInfo&) const override;

    /// @copydoc MasterSlaveConstraint::SetDofList
    /// @note Slave DoFs are appended to the list of master DoFs.
    void SetDofList(const DofPointerVectorType& rSlaveDofs,
                    const DofPointerVectorType& rMastertDofs,
                    const ProcessInfo&) override;

    /// @copydoc MasterSlaveConstraint::EquationIdVector
    /// @note Only master DoFs are populated.
    void EquationIdVector(EquationIdVectorType& rSlaveDoFs,
                          EquationIdVectorType& rMasterDoFs,
                          const ProcessInfo& rProcessInfo) const override;

    /// @brief copydoc MasterSlaveConstraint::GetSlaveDofsVector
    const DofPointerVectorType& GetSlaveDofsVector() const override
    {return mDummy;}

    const DofPointerVectorType& GetMasterDofsVector() const override
    {return mDofs;}

    void SetMasterDofsVector(const DofPointerVectorType& rDofs) override
    {mDofs = rDofs;}


    /// This function is NoOp.
    void ResetSlaveDofs(const ProcessInfo&) override {}

    /// @copydoc MasterSlaveConstraint::Apply
    void Apply(const ProcessInfo&) override {}

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

protected:
    const DofArray& GetDofs() const noexcept
    {
        return mDofs;
    }

    const MatrixType& GetHessian() const
    {
        return this->GetData().GetValue(GEOMETRIC_STIFFNESS_MATRIX);
    }

    MatrixType& GetHessian()
    {
        return this->Data().GetValue(GEOMETRIC_STIFFNESS_MATRIX);
    }

private:
    friend class Serializer;

    void save(Serializer& rSerializer) const override;

    void load(Serializer& rDeserializer) override;

    DofArray mDofs, mDummy;
}; // class MultifreedomConstraint


} // namespace Kratos
