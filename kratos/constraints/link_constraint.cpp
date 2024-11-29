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
#include "constraints/link_constraint.h" // LinkConstraint
#include "includes/define.h"
#include "includes/process_info.h"
#include "includes/variables.h" // DISPLACEMENT_X, DISPLACEMENT_Y, DISPLACEMENT_Z
#include "utilities/atomic_utilities.h" // AtomicAdd

// STL includes
#include <algorithm> // std::max_element, std::equal, std::transform


namespace Kratos {


struct LinkConstraint::Impl
{
    using ValueType = double;

    static void ComputeRelationMatrix(LinkConstraint::MatrixType& rRelationMatrix,
                                      LinkConstraint::VectorType& rConstraintGaps,
                                      LinkConstraint::DofPointerVectorType& rSlaves,
                                      LinkConstraint::DofPointerVectorType& rMasters,
                                      const std::array<Node*, 2> rNodePair,
                                      const std::array<ValueType,6>& rLastPositions,
                                      const std::array<ValueType,6>& rLastDisplacements,
                                      const unsigned Dimensions)
    {
        // Make sure we don't try to allocate an underflown unsigned integer.
        KRATOS_ERROR_IF_NOT(Dimensions);

        // Static stuff.
        const std::array<const Variable<ValueType>*,3> displacement_components {
            &DISPLACEMENT_X,
            &DISPLACEMENT_Y,
            &DISPLACEMENT_Z
        };

        // Collect the set of related Dofs.
        // This constraint involves exactly 2 nodes, only applies to DISPLACEMENT,
        // and supports 3 dimensions at most. This means that the number of Dofs
        // is maximum 6. Dofs are inserted in the following order:
        // {u^0_x, ..., u^0_w, u^1_x, ..., u^1_w}
        // Where w is the axis with the highest index for the provided dimensions
        // (e.g.: z-axis in 3 dimensions). Unused components are left uninitialized.
        using DofType = Dof<ValueType>;
        std::array<DofType*, 6> dofs;

        // Dof coefficients in the constraint equation, in the
        // same order as "dofs".
        std::array<ValueType,6> constraint_equation {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        ValueType constraint_gap = 0;

        for (unsigned i_component=0u; i_component<Dimensions; ++i_component) {
            std::array<ValueType,2> displacement_component_pair;

            for (unsigned i_node=0u; i_node<2; ++i_node) {
                // For some ungodly reason, the MasterSlaveConstraint interface demands mutable
                // Dofs, but node just refuses to provide mutable access to individual ones.
                // What it does provide though, is mutable access to ALL its Dofs at once. Insane.
                const unsigned i_dof_in_node = rNodePair[i_node]->GetDofPosition(*displacement_components[i_component]);
                const unsigned iDof = i_node ? i_component + Dimensions : i_component;
                dofs[iDof] = rNodePair[i_node]->GetDofs()[i_dof_in_node].get();
            } // for i_node in range(2)

            constraint_equation[i_component] = rLastPositions[i_component] + rLastDisplacements[i_component]
                                               - rLastPositions[i_component + Dimensions] - rLastDisplacements[i_component + Dimensions];
            constraint_equation[i_component + Dimensions] = -constraint_equation[i_component];
            constraint_gap -= std::pow(rLastDisplacements[i_component] - rLastDisplacements[i_component + Dimensions], 2);
        } // for i_component in range(mDimensions)

        // Pick a slave.
        // The only hard criterion is that the slave's coefficient must
        // not vanish, but that's not a sufficiently well-defined rule
        // to unambiguously pick a slave Dof. So let's pick the Dof with
        // the largest absolute coefficient.
        const auto it_slave = std::max_element(constraint_equation.begin(),
                                               constraint_equation.end(),
                                               [](const ValueType left, const ValueType right) -> bool {
                                                   return std::abs(left) < std::abs(right);
                                               });
        const unsigned i_slave = std::distance(constraint_equation.begin(), it_slave);
        const ValueType slave_coefficient = constraint_equation[i_slave];

        constexpr ValueType tolerance = 1e-16;
        KRATOS_ERROR_IF(std::abs(slave_coefficient) < tolerance)
            << "coincident nodes in LinkConstraint";

        // Set output sizes.
        rRelationMatrix.resize(1, 2 * Dimensions - 1);
        rConstraintGaps.resize(1);
        rSlaves.resize(1);
        rMasters.resize(2 * Dimensions - 1);

        // Fill the relation matrix and Dof vectors.
        unsigned i_master = 0u;
        for (unsigned i_dof=0u; i_dof<2 * Dimensions; ++i_dof) {
            if (i_dof == i_slave) {
                rSlaves.front() = dofs[i_slave];
            } /*if i_master == i_slave*/ else {
                rMasters[i_master] = dofs[i_dof];
                rRelationMatrix(0, i_master) = -constraint_equation[i_dof] / slave_coefficient;
                ++i_master;
            }
        } // for i_dof in range(2 * Dimensions)

        rConstraintGaps[0] = constraint_gap;

        std::cout << "constraint equation: ";
        for (const auto c : constraint_equation) std::cout << c << ' ';
        std::cout << "\n";
    }


    void FetchNodePositions(std::array<ValueType,6>& rPositions,
                            std::array<ValueType,6>& rDisplacements) const {
        KRATOS_TRY

        const std::array<const Variable<Impl::ValueType>*,3> displacement_components {
            &DISPLACEMENT_X,
            &DISPLACEMENT_Y,
            &DISPLACEMENT_Z
        };

        for (unsigned i_component=0u; i_component<mDimensions; ++i_component) {
            for (unsigned i_node=0u; i_node<2; ++i_node) {
                const unsigned i_dof = i_component + (i_node ? mDimensions : 0u);
                rPositions[i_dof] = mNodePair[i_node]->Coordinates()[i_component];
                rDisplacements[i_dof] = mNodePair[i_node]->GetSolutionStepValue(*displacement_components[i_component]);
            } // for i_node in range(2)
        } // for i_component in range(mDimensions)

        KRATOS_CATCH("")
    }


    unsigned mDimensions;

    /// @details MasterSlaveConstraint::GetSlaveDofsVector and MasterSlaveConstraint::GetMasterDofsVector
    ///          require arrays of mutable Dof pointers, which are only obtainable from mutable nodes,
    ///          so the nodes pointers stored here cannot be immutable. Risky business.
    std::array<Node*,2> mNodePair;

    std::array<ValueType,6> mLastPositions;

    std::array<ValueType,6> mLastDisplacements;

    /// @details Unfortunately, the MasterSlaveConstraint interface demands that
    ///          GetSlaveDofsVector and GetMasterDofsVector returns the array of
    ///          Dofs by reference, even though those must be computed dynamically
    ///          by LinkConstraint. As a result, these vectors must be stored as
    ///          member variables and updated in InitializeNonlinearIteration, instead
    ///          of being computed on the fly.
    std::array<DofPointerVectorType,2> mDofVectors; // {slave, masters}

    LinkConstraint::MatrixType mRelationMatrix;

    LinkConstraint::VectorType mConstraintGaps;
}; // struct LinkConstraint::Impl


LinkConstraint::LinkConstraint(const IndexType Id,
                               Node& rFirst,
                               Node& rSecond,
                               const unsigned Dimensions)
    : MasterSlaveConstraint(Id),
      mpImpl(new Impl{Dimensions,
                      {&rFirst, &rSecond},
                      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                      {},
                      {},
                      {}})
{
}


LinkConstraint::LinkConstraint(LinkConstraint&& rRhs) noexcept
    : mpImpl(std::move(rRhs.mpImpl))
{
}


LinkConstraint::LinkConstraint(const LinkConstraint& rRhs)
    : mpImpl(new Impl(*rRhs.mpImpl))
{
}


LinkConstraint& LinkConstraint::operator=(LinkConstraint&& rRhs) noexcept
{
    mpImpl = std::move(rRhs.mpImpl);
    return *this;
}


LinkConstraint& LinkConstraint::operator=(const LinkConstraint& rRhs)
{
    mpImpl.reset(new Impl(*rRhs.mpImpl));
    return *this;
}


LinkConstraint::~LinkConstraint() = default;


MasterSlaveConstraint::Pointer LinkConstraint::Clone(IndexType NewId) const
{
    KRATOS_TRY
    auto p_other = LinkConstraint::Pointer(
        new LinkConstraint(NewId,
                           *mpImpl->mNodePair.front(),
                           *mpImpl->mNodePair.back(),
                           mpImpl->mDimensions)
    );
    return p_other;
    KRATOS_CATCH("")
}


void LinkConstraint::Initialize([[maybe_unused]] const ProcessInfo&)
{
    mpImpl->FetchNodePositions(mpImpl->mLastPositions,
                               mpImpl->mLastDisplacements);
    Impl::ComputeRelationMatrix(mpImpl->mRelationMatrix,
                                mpImpl->mConstraintGaps,
                                mpImpl->mDofVectors.front(),
                                mpImpl->mDofVectors.back(),
                                mpImpl->mNodePair,
                                mpImpl->mLastPositions,
                                mpImpl->mLastDisplacements,
                                mpImpl->mDimensions);
}


void LinkConstraint::InitializeSolutionStep(const ProcessInfo& rProcessInfo)
{
    KRATOS_TRY
    // Store displacements before nonlinear iterations begin.
    // At this point, displacements are at their last converged
    // state, when the constraint is supposed to be still satisfied.
    mpImpl->FetchNodePositions(mpImpl->mLastPositions,
                               mpImpl->mLastDisplacements);
    this->InitializeNonLinearIteration(rProcessInfo);
    KRATOS_CATCH("")
}


void LinkConstraint::InitializeNonLinearIteration([[maybe_unused]] const ProcessInfo&)
{
    KRATOS_TRY
    [[maybe_unused]] std::array<Impl::ValueType,6> dummy;
    mpImpl->FetchNodePositions(dummy, mpImpl->mLastDisplacements);
    Impl::ComputeRelationMatrix(mpImpl->mRelationMatrix,
                                mpImpl->mConstraintGaps,
                                mpImpl->mDofVectors.front(),
                                mpImpl->mDofVectors.back(),
                                mpImpl->mNodePair,
                                mpImpl->mLastPositions,
                                mpImpl->mLastDisplacements,
                                mpImpl->mDimensions);
    KRATOS_CATCH("")
}


void LinkConstraint::FinalizeNonLinearIteration([[maybe_unused]] const ProcessInfo&)
{
    //const std::array<const Variable<Impl::ValueType>*,3> displacement_components {
    //    &DISPLACEMENT_X,
    //    &DISPLACEMENT_Y,
    //    &DISPLACEMENT_Z
    //};
    //Impl::ValueType distance = 0;
    //for (unsigned i_component=0u; i_component<mpImpl->mDimensions; ++i_component) {
    //    distance += std::pow(mpImpl->mNodePair[0]->Coordinates()[i_component] + mpImpl->mNodePair[0]->GetSolutionStepValue(*displacement_components[i_component])
    //                         - mpImpl->mNodePair[1]->Coordinates()[i_component] - mpImpl->mNodePair[1]->GetSolutionStepValue(*displacement_components[i_component]), 2);
    //}
    //std::cout << "LinkConstraint: distance is " << distance << "\n";
}


void LinkConstraint::ResetSlaveDofs(const ProcessInfo& rProcessInfo)
{
    for (Dof<Impl::ValueType>* p_dof : mpImpl->mDofVectors.front()) {
        p_dof->GetSolutionStepValue() = 0;
    }
}


void LinkConstraint::Apply(const ProcessInfo& rProcessInfo)
{
    KRATOS_TRY
    Vector master_dofs_values(mpImpl->mDofVectors.back().size());

    for (IndexType i = 0; i < mpImpl->mDofVectors.back().size(); ++i) {
        master_dofs_values[i] = mpImpl->mDofVectors.back()[i]->GetSolutionStepValue();
    }

    // Apply the constraint to the slave dofs
    for (IndexType i = 0; i < mpImpl->mRelationMatrix.size1(); ++i) {
        Impl::ValueType tmp = mpImpl->mConstraintGaps[0];
        for(IndexType j = 0; j < mpImpl->mRelationMatrix.size2(); ++j) {
            tmp += mpImpl->mRelationMatrix(i,j) * master_dofs_values[j];
        }

        AtomicAdd(mpImpl->mDofVectors.front()[i]->GetSolutionStepValue(), tmp);
    }
    KRATOS_CATCH("")
}


void LinkConstraint::GetDofList(DofPointerVectorType& rSlaves,
                                DofPointerVectorType& rMasters,
                                [[maybe_unused]] const ProcessInfo& rProcessInfo) const
{
    rSlaves = this->GetSlaveDofsVector();
    rMasters = this->GetMasterDofsVector();
}


void LinkConstraint::EquationIdVector(EquationIdVectorType& rSlaves,
                                      EquationIdVectorType& rMasters,
                                      [[maybe_unused]] const ProcessInfo& rProcessInfo) const
{
    KRATOS_TRY
    rSlaves.resize(mpImpl->mDofVectors.front().size());
    rMasters.resize(mpImpl->mDofVectors.back().size());

    const auto dof_to_equation_id = [](const Dof<double>* pDof) -> Dof<double>::EquationIdType {
        return pDof->EquationId();
    };

    std::transform(mpImpl->mDofVectors.front().begin(),
                   mpImpl->mDofVectors.front().end(),
                   rSlaves.begin(),
                   dof_to_equation_id);

    std::transform(mpImpl->mDofVectors.back().begin(),
                   mpImpl->mDofVectors.back().end(),
                   rMasters.begin(),
                   dof_to_equation_id);
    KRATOS_CATCH("")
}


const LinkConstraint::DofPointerVectorType& LinkConstraint::GetSlaveDofsVector() const
{
    return mpImpl->mDofVectors.front();
}


const LinkConstraint::DofPointerVectorType& LinkConstraint::GetMasterDofsVector() const
{
    return mpImpl->mDofVectors.back();
}


void LinkConstraint::CalculateLocalSystem(MatrixType& rRelationMatrix,
                                          VectorType& rConstraintGaps,
                                          [[maybe_unused]] const ProcessInfo& rProcessInfo) const
{
    rRelationMatrix = mpImpl->mRelationMatrix;
    rConstraintGaps = mpImpl->mConstraintGaps;
}


int LinkConstraint::Check(const ProcessInfo& rProcessInfo) const
{
    // Check whatever the base class checks.
    MasterSlaveConstraint::Check(rProcessInfo);

    // Restrict dimensions to 1, 2, or 3.
    KRATOS_ERROR_IF_NOT(0 < mpImpl->mDimensions || mpImpl->mDimensions <= 3) << "unsupported dimensions (" << mpImpl->mDimensions << ")";

    // Make sure that the nodes have the necessary Dofs.
    const std::array<const Variable<Impl::ValueType>*,3> displacement_components {
        &DISPLACEMENT_X,
        &DISPLACEMENT_Y,
        &DISPLACEMENT_Z
    };

    for (unsigned i_component=0u; i_component<mpImpl->mDimensions; ++i_component) {
        for (unsigned i_node=0u; i_node<2; ++i_node) {
            KRATOS_ERROR_IF_NOT(mpImpl->mNodePair[i_node]->HasDofFor(*displacement_components[i_component]))
                << "node " << mpImpl->mNodePair[i_node]->Id()
                << " has no Dof for " << displacement_components[i_component]->Name();
        } // for i_node in range(2)
    } // for i_component in range(mDimensions)

    // Check for overlapping nodes.
    constexpr Impl::ValueType tolerance = 1e-16;
    KRATOS_ERROR_IF(norm_2(mpImpl->mNodePair.front()->Coordinates() - mpImpl->mNodePair.back()->Coordinates()) < tolerance)
        << "coincident nodes in LinkConstraint";

    return 0;
}

} // namespace Kratos
