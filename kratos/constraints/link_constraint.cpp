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
#include "includes/variables.h" // DISPLACEMENT_X, DISPLACEMENT_Y, DISPLACEMENT_Z

// STL includes
#include <algorithm> // std::max_element, std::equal, std::transform


namespace Kratos {


namespace detail {
void ComputeRelationMatrix(LinkConstraint::MatrixType& rRelationMatrix,
                           LinkConstraint::DofPointerVectorType& rSlaves,
                           LinkConstraint::DofPointerVectorType& rMasters,
                           const std::array<Node*, 2> rNodePair,
                           const unsigned Dimensions)
{
    std::cout << "got to " << KRATOS_CODE_LOCATION.GetLineNumber() << std::endl;
    // Set output sizes.
    rRelationMatrix.resize(1, 2 * Dimensions - 1);
    rSlaves.resize(1);
    rMasters.resize(2 * Dimensions - 1);

    // Static stuff.
    using ValueType = double;
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
    std::array<ValueType,6> constraint_equation;

    for (unsigned iComponent=0u; iComponent<Dimensions; ++iComponent) {
        std::array<ValueType,2> initial_coordinate_pair, displacement_component_pair;

        for (unsigned iNode=0u; iNode<2; ++iNode) {
            const unsigned iDof = iNode ? iComponent + Dimensions : iComponent;

            // For some ungodly reason, the MasterSlaveConstraint interface demands mutable
            // Dofs, but node just refuses to provide mutable access to individual Dofs.
            // What it does provide though, is mutable access to ALL its Dofs at once. Insane.
            const unsigned i_dof_in_node = rNodePair[iNode]->GetDofPosition(*displacement_components[iComponent]);
            dofs[iDof] = rNodePair[iNode]->GetDofs()[i_dof_in_node].get();

            initial_coordinate_pair[iNode] = rNodePair[iNode]->Coordinates()[iComponent];
            displacement_component_pair[iNode] = rNodePair[iNode]->GetSolutionStepValue(*displacement_components[iComponent]);
        } // for iNode in range(2)

        constraint_equation[iComponent] = displacement_component_pair[0]
                                          + 2 * (initial_coordinate_pair[0] - initial_coordinate_pair[1])
                                          - displacement_component_pair[1];
        constraint_equation[iComponent + Dimensions] = -constraint_equation[iComponent];
    } // for iComponent in range(mDimensions)

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
    std::cout << "got to " << KRATOS_CODE_LOCATION.GetLineNumber() << std::endl;
}


// Some wonderful strategies require a master/slave partitioning before
// LinkConstraint ever has a chance of producing an up-to-date one, so
// it must compute one that might get outdated once solution takes place ...
void ProducePartitioning(LinkConstraint::DofPointerVectorType& rSlaves,
                         LinkConstraint::DofPointerVectorType& rMasters,
                         const std::array<Node*,2>& rNodePair,
                         const unsigned Dimensions)
{
    [[maybe_unused]] LinkConstraint::MatrixType matrix;
    ComputeRelationMatrix(matrix, rSlaves, rMasters, rNodePair, Dimensions);
}
} // namespace detail


LinkConstraint::LinkConstraint(const IndexType Id,
                               Node& rFirst,
                               Node& rSecond,
                               const unsigned Dimensions)
    : MasterSlaveConstraint(Id),
      mDimensions(Dimensions),
      mNodePair({&rFirst, &rSecond}),
      mDofVectors()
{
    KRATOS_ERROR_IF_NOT(Dimensions);
    std::cout << "got to " << KRATOS_CODE_LOCATION.GetLineNumber() << std::endl;
    mDofVectors.emplace();
    detail::ProducePartitioning(mDofVectors.value()[0], mDofVectors.value()[1], mNodePair.value(), mDimensions);
}


MasterSlaveConstraint::Pointer LinkConstraint::Clone(IndexType NewId) const
{
    std::cout << "got to " << KRATOS_CODE_LOCATION.GetLineNumber() << std::endl;
    KRATOS_TRY

    const auto& r_node_pair = mNodePair.value();
    auto p_other = LinkConstraint::Pointer(
        new LinkConstraint(NewId,
                           *r_node_pair.front(),
                           *r_node_pair.back(),
                           mDimensions)
    );
    p_other->mDofVectors = mDofVectors;
    std::cout << "got to " << KRATOS_CODE_LOCATION.GetLineNumber() << std::endl;
    return p_other;

    KRATOS_CATCH("")
}


void LinkConstraint::InitializeSolutionStep(const ProcessInfo& rProcessInfo)
{
    KRATOS_TRY
    std::cout << "got to " << KRATOS_CODE_LOCATION.GetLineNumber() << std::endl;
    this->InitializeNonLinearIteration(rProcessInfo);
    std::cout << "got to " << KRATOS_CODE_LOCATION.GetLineNumber() << std::endl;
    KRATOS_CATCH("")
}


void LinkConstraint::InitializeNonLinearIteration(const ProcessInfo& rProcessInfo)
{
    std::cout << "got to " << KRATOS_CODE_LOCATION.GetLineNumber() << std::endl;
    KRATOS_TRY
    // Since MasterSlaveConstraint::GetSlaveDofsVector and MasterSlaveConstraint::GetMasterDofsVector
    // require a return by const reference, these vectors must be stored in this class instead
    // of computed on the fly. This also means that this calculation must happen outside these functions
    // because they're const ... so here we are.
    [[maybe_unused]] MatrixType relation_matrix;
    mDofVectors.emplace();
    detail::ComputeRelationMatrix(relation_matrix,
                                  mDofVectors.value()[0],
                                  mDofVectors.value()[1],
                                  mNodePair.value(),
                                  mDimensions);
    KRATOS_CATCH("")
    std::cout << "got to " << KRATOS_CODE_LOCATION.GetLineNumber() << std::endl;
}


void LinkConstraint::GetDofList(DofPointerVectorType& rSlaves,
                                DofPointerVectorType& rMasters,
                                const ProcessInfo& rProcessInfo) const
{
    std::cout << "got to " << KRATOS_CODE_LOCATION.GetLineNumber() << std::endl;
    KRATOS_TRY
    if (mDofVectors.has_value()) {
        rSlaves.resize(mDofVectors.value()[0].size());
        rMasters.resize(mDofVectors.value()[1].size());
    } else {
        detail::ProducePartitioning(rSlaves, rMasters, mNodePair.value(), mDimensions);
    }

    std::copy(mDofVectors.value()[0].begin(), mDofVectors.value()[0].end(), rSlaves.begin());
    std::copy(mDofVectors.value()[1].begin(), mDofVectors.value()[1].end(), rMasters.begin());
    std::cout << "got to " << KRATOS_CODE_LOCATION.GetLineNumber() << std::endl;
    KRATOS_CATCH("")
}


void LinkConstraint::EquationIdVector(EquationIdVectorType& rSlaves,
                                      EquationIdVectorType& rMasters,
                                      const ProcessInfo& rProcessInfo) const
{
    std::cout << "got to " << KRATOS_CODE_LOCATION.GetLineNumber() << std::endl;
    KRATOS_TRY
    rSlaves.resize(mDofVectors.value()[0].size());
    rMasters.resize(mDofVectors.value()[1].size());

    const auto dof_to_equation_id = [](const Dof<double>* pDof) -> Dof<double>::EquationIdType {
        return pDof->EquationId();
    };

    std::transform(mDofVectors.value()[0].begin(),
                   mDofVectors.value()[0].end(),
                   rSlaves.begin(),
                   dof_to_equation_id);

    std::transform(mDofVectors.value()[1].begin(),
                   mDofVectors.value()[1].end(),
                   rMasters.begin(),
                   dof_to_equation_id);
    std::cout << "got to " << KRATOS_CODE_LOCATION.GetLineNumber() << std::endl;
    KRATOS_CATCH("")
}


const LinkConstraint::DofPointerVectorType& LinkConstraint::GetSlaveDofsVector() const
{
    return mDofVectors.value()[0];
}


const LinkConstraint::DofPointerVectorType& LinkConstraint::GetMasterDofsVector() const
{
    std::cout << "got to " << KRATOS_CODE_LOCATION.GetLineNumber() << std::endl;
    return mDofVectors.value()[1];
}


void LinkConstraint::CalculateLocalSystem(MatrixType& rRelationMatrix,
                                          VectorType& rConstraintGaps,
                                          const ProcessInfo& rProcessInfo) const
{
    std::cout << "got to " << KRATOS_CODE_LOCATION.GetLineNumber() << std::endl;
    DofPointerVectorType slaves, masters;
    detail::ComputeRelationMatrix(rRelationMatrix,
                                  slaves,
                                  masters,
                                  mNodePair.value(),
                                  mDimensions);

    // The Dof partitioning into masters/slave should not have changed since
    // it was last computed (either in InitializeSolutionStep or InitializeNonLinearItertation).
    KRATOS_DEBUG_ERROR_IF_NOT(std::equal(slaves.begin(), slaves.end(), mDofVectors.value()[0].begin()));
    KRATOS_DEBUG_ERROR_IF_NOT(std::equal(masters.begin(), masters.end(), mDofVectors.value()[1].begin()));
}


int LinkConstraint::Check(const ProcessInfo& rProcessInfo) const
{
    std::cout << "got to " << KRATOS_CODE_LOCATION.GetLineNumber() << std::endl;
    // Check whatever the base class checks.
    MasterSlaveConstraint::Check(rProcessInfo);

    // Restrict dimensions to 1, 2, or 3.
    KRATOS_ERROR_IF_NOT(0 < mDimensions || mDimensions <= 3) << "unsupported dimensions (" << mDimensions << ")";

    // Make sure this instance was not default constructed.
    KRATOS_ERROR_IF_NOT(mNodePair.has_value());
    const auto& r_node_pair = mNodePair.value();

    // Make sure that the nodes have the necessary Dofs.
    using ValueType = double;
    const std::array<const Variable<ValueType>*,3> displacement_components {
        &DISPLACEMENT_X,
        &DISPLACEMENT_Y,
        &DISPLACEMENT_Z
    };

    for (unsigned i_component=0u; i_component<mDimensions; ++i_component) {
        for (unsigned i_node=0u; i_node<2; ++i_node) {
            KRATOS_ERROR_IF_NOT(r_node_pair[i_node]->HasDofFor(*displacement_components[i_component]))
                << "node " << r_node_pair[i_node]->Id()
                << " has no Dof for " << displacement_components[i_component]->Name();
        } // for i_node in range(2)
    } // for i_component in range(mDimensions)

    // Check for overlapping nodes.
    constexpr ValueType tolerance = 1e-16;
    KRATOS_ERROR_IF(norm_2(mNodePair.value()[0]->Coordinates() - mNodePair.value()[1]->Coordinates()) < tolerance)
        << "coincident nodes in LinkConstraint";

    std::cout << "got to " << KRATOS_CODE_LOCATION.GetLineNumber() << std::endl;
    return 0;
}

} // namespace Kratos
