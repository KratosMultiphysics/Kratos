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

// Project includes
#include "custom_constraints/link_constraint.hpp" // LinkConstraint
#include "includes/define.h"
#include "includes/process_info.h"
#include "includes/variables.h" // DISPLACEMENT_X, DISPLACEMENT_Y, DISPLACEMENT_Z
#include "includes/serializer.h" // Serializer
#include "includes/variables.h" // CONSTITUTIVE_MATRIX, INTERNAL_FORCES_VECTOR

// STL includes
#include <algorithm> // std::max_element, std::equal, std::transform


namespace Kratos {


struct LinkConstraint::Impl
{
    using ValueType = double;

    static void ComputeRelationMatrix(LinkConstraint::MatrixType& rRelationMatrix,
                                      LinkConstraint::MatrixType& rHessian,
                                      LinkConstraint::VectorType& rConstraintGaps,
                                      const std::array<ValueType,6>& rLastPositions,
                                      const std::array<ValueType,6>& rLastDisplacements,
                                      const unsigned Dimensions)
    {
        // Make sure we don't try to allocate an underflown unsigned integer.
        KRATOS_ERROR_IF_NOT(Dimensions);

        // Set output sizes.
        rRelationMatrix.resize(1, 2 * Dimensions);
        rHessian.resize(2 * Dimensions, 2 * Dimensions);
        rConstraintGaps.resize(1);

        // Initialize outputs.
        rConstraintGaps[0] = 0.0;
        std::fill(rRelationMatrix.data().begin(),
                  rRelationMatrix.data().end(),
                  static_cast<ValueType>(0));
        std::fill(rHessian.data().begin(),
                  rHessian.data().end(),
                  static_cast<ValueType>(0));

        // Compute constraint violation.
        double initial_norm = 0.0, current_norm = 0.0;
        for (std::size_t i_component=0ul; i_component<Dimensions; ++i_component) {
            const auto initial_diff = rLastPositions[i_component] - rLastPositions[i_component + Dimensions];
            const auto current_diff = initial_diff + rLastDisplacements[i_component] - rLastDisplacements[i_component + Dimensions];
            initial_norm += initial_diff * initial_diff;
            current_norm += current_diff * current_diff;
        } // for i_component in range(Dimensions)
        const auto constraint_loss = current_norm - initial_norm;
        const auto hessian_scale = 2.0 * constraint_loss;

        // Dof order:
        // {u^0_x, ..., u^0_w, u^1_x, ..., u^1_w}
        // Where w is the axis with the highest index for the provided dimensions
        // (e.g.: z-axis in 3 dimensions).

        for (std::size_t i_component=0u; i_component<Dimensions; ++i_component) {
            const auto value = rLastPositions[i_component] - rLastPositions[i_component + Dimensions]
                             + rLastDisplacements[i_component] - rLastDisplacements[i_component + Dimensions];
            rRelationMatrix(0, i_component) = 2 * value;
            rRelationMatrix(0, i_component + Dimensions) = -2 * value;

            rHessian(i_component, i_component) = hessian_scale;
            rHessian(i_component, i_component + Dimensions) = -hessian_scale;
            rHessian(i_component + Dimensions, i_component) = -hessian_scale;
            rHessian(i_component + Dimensions, i_component + Dimensions) = hessian_scale;
        } // for i_component in range(mDimensions)
    }


    static LinkConstraint::DofPointerVectorType CollectDofs(Node& rFirst,
                                                            Node& rSecond,
                                                            const std::size_t Dimensions)
    {
        DofPointerVectorType dofs;

        KRATOS_TRY
        using ValueType = double;

        const std::array<const Variable<ValueType>*,3> displacement_components {
            &DISPLACEMENT_X,
            &DISPLACEMENT_Y,
            &DISPLACEMENT_Z
        };

        std::array<Node*, 2> nodes {&rFirst, &rSecond};
        dofs.resize(2 * Dimensions);

        // Collect dofs from nodes in the following order:
        // {u^0_x, ..., u^0_w, u^1_x, ..., u^1_w}
        // Where w is the axis with the highest index for the provided dimensions
        // (e.g.: z-axis in 3 dimensions).
        for (std::size_t i_dimension=0ul; i_dimension<Dimensions; ++i_dimension) {
            const auto component_id = displacement_components[i_dimension]->Key();

            for (std::size_t i_node=0ul; i_node<2; ++i_node) {
                Node& r_node = *nodes[i_node];
                const auto it_dof = std::find_if(r_node.GetDofs().begin(),
                                                 r_node.GetDofs().end(),
                                                 [component_id](auto& rp_dof){
                                                         return rp_dof->GetVariable().Key() == component_id;
                                                 });
                KRATOS_ERROR_IF(it_dof == r_node.GetDofs().end())
                    << "cannot find DoF " << displacement_components[i_dimension]->Name()
                    << " in node " << r_node.Id();
                dofs[i_node * Dimensions + i_dimension] = it_dof->get();
            } // for i_node in range(2)
        } // for i_dimension in range(Dimensions)
        KRATOS_CATCH("")

        return dofs;
    }


    void FetchNodePositions(std::array<ValueType,6>& rPositions,
                            std::array<ValueType,6>& rDisplacements) const {
        KRATOS_TRY

        const std::array<const Variable<Impl::ValueType>*,3> displacement_components {
            &DISPLACEMENT_X,
            &DISPLACEMENT_Y,
            &DISPLACEMENT_Z
        };

        for (std::size_t i_component=0u; i_component<mDimensions; ++i_component) {
            for (std::size_t i_node=0u; i_node<2; ++i_node) {
                const std::size_t i_dof = i_component + (i_node ? mDimensions : 0u);
                rPositions[i_dof] = mNodePair[i_node]->Coordinates()[i_component];
                rDisplacements[i_dof] = mNodePair[i_node]->GetSolutionStepValue(*displacement_components[i_component]);
            } // for i_node in range(2)
        } // for i_component in range(mDimensions)

        KRATOS_CATCH("")
    }


    double GetSquareDistance() const {
        double output = 0.0;
        std::array<Impl::ValueType,6> positions, displacements;
        this->FetchNodePositions(positions, displacements);
        std::transform(positions.begin(), positions.end(), displacements.begin(), positions.begin(), std::plus<double>());
        for (std::size_t i_component=0ul; i_component<mDimensions; ++i_component) {
            const auto diff = positions[i_component] - positions[i_component + mDimensions];
            output += diff * diff;
        }
        return output;
    }


    std::size_t mDimensions;

    bool mIsMeshMoved;

    /// @details MasterSlaveConstraint::GetSlaveDofsVector and MasterSlaveConstraint::GetMasterDofsVector
    ///          require arrays of mutable Dof pointers, which are only obtainable from mutable nodes,
    ///          so the nodes pointers stored here cannot be immutable. Risky business.
    std::array<Node*,2> mNodePair;

    std::array<ValueType,6> mLastPositions;

    std::array<ValueType,6> mLastDisplacements;
}; // struct LinkConstraint::Impl


LinkConstraint::LinkConstraint(const IndexType Id,
                               Node& rFirst,
                               Node& rSecond,
                               const std::size_t Dimensions,
                               bool IsMeshMoved)
    : MultifreedomConstraint(/*Id               : */Id,
                             /*rDofs            : */Impl::CollectDofs(rFirst, rSecond, Dimensions),
                             /*rConstraintLabels: */std::vector<std::size_t> {Id}),
      mpImpl(new Impl{/*mDimensions         : */Dimensions,
                      /*mIsMeshMoved        : */IsMeshMoved,
                      /*mNodePair           : */{&rFirst, &rSecond},
                      /*mLastPositions      : */{0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                      /*mLastDisplacements  : */{0.0, 0.0, 0.0, 0.0, 0.0, 0.0}})
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
                           mpImpl->mDimensions,
                           mpImpl->mIsMeshMoved)
    );
    return p_other;
    KRATOS_CATCH("")
}


void LinkConstraint::Initialize([[maybe_unused]] const ProcessInfo&)
{
}


void LinkConstraint::InitializeSolutionStep(const ProcessInfo&)
{
    KRATOS_TRY
    // Store displacements before nonlinear iterations begin.
    // At this point, displacements are at their last converged
    // state, when the constraint is supposed to be still satisfied.
    mpImpl->FetchNodePositions(mpImpl->mLastPositions,
                               mpImpl->mLastDisplacements);

    if (!mpImpl->mIsMeshMoved)
        for (std::size_t i_component=0ul; i_component<mpImpl->mLastPositions.size(); ++i_component)
            mpImpl->mLastPositions[i_component] += mpImpl->mLastDisplacements[i_component];

    Impl::ComputeRelationMatrix(this->Data().GetValue(CONSTITUTIVE_MATRIX),
                                this->GetHessian(),
                                this->Data().GetValue(INTERNAL_FORCES_VECTOR),
                                mpImpl->mLastPositions,
                                mpImpl->mLastDisplacements,
                                mpImpl->mDimensions);
    KRATOS_CATCH("")
}


void LinkConstraint::InitializeNonLinearIteration(const ProcessInfo&)
{
    KRATOS_TRY
    [[maybe_unused]] std::array<Impl::ValueType,6> dummy;
    mpImpl->FetchNodePositions(dummy,  mpImpl->mLastDisplacements);

    Impl::ComputeRelationMatrix(this->Data().GetValue(CONSTITUTIVE_MATRIX),
                                this->GetHessian(),
                                this->Data().GetValue(INTERNAL_FORCES_VECTOR),
                                mpImpl->mLastPositions,
                                mpImpl->mLastDisplacements,
                                mpImpl->mDimensions);
    KRATOS_CATCH("")
}


void LinkConstraint::FinalizeNonLinearIteration(const ProcessInfo& rProcessInfo)
{
}


void LinkConstraint::FinalizeSolutionStep(const ProcessInfo&)
{
    this->Data().GetValue(CONSTITUTIVE_MATRIX).clear();
    this->Data().GetValue(INTERNAL_FORCES_VECTOR).clear();
}


void LinkConstraint::Finalize(const ProcessInfo&)
{
}


void LinkConstraint::CalculateLocalSystem(MatrixType& rRelationMatrix,
                                          VectorType& rConstraintGaps,
                                          [[maybe_unused]] const ProcessInfo& rProcessInfo) const
{
    rRelationMatrix = this->GetData().GetValue(CONSTITUTIVE_MATRIX);
    rConstraintGaps = this->GetData().GetValue(INTERNAL_FORCES_VECTOR);
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


void LinkConstraint::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, MultifreedomConstraint);
    rSerializer.save("Dimensions", mpImpl->mDimensions);
    rSerializer.save("Nodes", mpImpl->mNodePair);
    // The rest of the private members should be computed/cached after restarting.
}


void LinkConstraint::load(Serializer& rDeserializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rDeserializer, MultifreedomConstraint);
    rDeserializer.load("Dimensions", mpImpl->mDimensions);
    rDeserializer.load("Nodes", mpImpl->mNodePair);
}


std::ostream& operator<<(std::ostream& rStream, const LinkConstraint& rInstance)
{
    return rStream << "LinkConstraint between nodes "
                   << rInstance.mpImpl->mNodePair.front()->Id()
                   << " and "
                   << rInstance.mpImpl->mNodePair.back()->Id();
}


} // namespace Kratos
