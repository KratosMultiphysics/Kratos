//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

// System includes

// External includes

// Project includes
#include "geometries/geometry.h"
#include "includes/define.h"
#include "includes/kratos_flags.h"
#include "includes/node.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"
#include "utilities/coordinate_transformation_utilities.h"

// Application includes
#include "custom_utilities/fluid_calculation_utilities.h"
#include "fluid_dynamics_application_variables.h"

// Include base h
#include "fluid_adjoint_slip_utilities.h"

namespace Kratos
{

FluidAdjointSlipUtilities::FluidAdjointSlipUtilities(
    const IndexType Dimension,
    const IndexType BlockSize)
    : mDimension(Dimension),
      mBlockSize(BlockSize),
      mRotationTool(Dimension, mBlockSize)
{
    KRATOS_TRY

    if (Dimension == 2) {
        this->mAddNodalRotationDerivativesMethod = &FluidAdjointSlipUtilities::TemplatedAddNodalRotationDerivatives<2>;
    } else if (Dimension == 3) {
        this->mAddNodalRotationDerivativesMethod = &FluidAdjointSlipUtilities::TemplatedAddNodalRotationDerivatives<3>;
    } else {
        KRATOS_ERROR << "Unsupported dimensionality requested. Only 2D and 3D "
                        "supported. [ Dimension = "
                     << Dimension << " ].\n";
    }

    KRATOS_CATCH("");
}

void FluidAdjointSlipUtilities::CalculateRotatedSlipConditionAppliedSlipVariableDerivatives(
    Matrix& rOutput,
    const Matrix& rResidualDerivatives,
    const GeometryType& rGeometry) const
{
    if (rOutput.size1() != rResidualDerivatives.size1() ||
        rOutput.size2() != rResidualDerivatives.size2()) {
        rOutput.resize(rResidualDerivatives.size1(), rResidualDerivatives.size2(), false);
    }

    rOutput.clear();

    const IndexType number_of_nodes = rGeometry.PointsNumber();

    // add residual derivative contributions
    for (IndexType a = 0; a < number_of_nodes; ++a) {
        const auto& r_node = rGeometry[a];
        const IndexType block_index = a * mBlockSize;
        if (r_node.Is(SLIP)) {
            AddNodalRotationDerivatives(rOutput, rResidualDerivatives, block_index, r_node);
            if (!r_node.IsFixed(ADJOINT_FLUID_VECTOR_1_X)) {
                // in here we have to check whether SLIP and not INLET
                // At the inlet ADJOINT_FLUID_VECTOR_1_X is fixed.
                // There can be a situation where a node is
                // shared among both inlet and slip (corner node). In this case
                // AddNodalApplySlipConditionDerivatives will set all the values in the
                // column corresponding to ADJOINT_FLUID_VECTOR_1_X of the common node to zero except for the
                // nodal normal derivative contributions. If the normal is perpendicular to the X direction, then
                // this will set the diagonal in the column corresponding to this node's ADJOINT_FLUID_VECTOR_1_X dof
                // to zero. But there can be non-zero off-diagonal due to normal not being zero. Eventhough
                // ADJOINT_FLUID_VECTOR_1_X will be fixed in the block builder and solver
                // it will not recognize this column as empty, so diagonal won't be fixed by B&S to 1.0.
                // Since ADJOINT_FLUID_VECTOR_1_X of this node is fixed, B&S will set all the off-diagonals
                // to zero. This will create a zero column in the system matrix which makes it singular.
                // So this check is performed here to avoid this.
                // In here inlets are alsways assumed to be dirichlet.
                AddNodalApplySlipConditionDerivatives(rOutput, block_index, r_node);
            }
        } else {
            AddNodalResidualDerivatives(rOutput, rResidualDerivatives, block_index);
        }
    }
}

void FluidAdjointSlipUtilities::CalculateRotatedSlipConditionAppliedNonSlipVariableDerivatives(
    Matrix& rOutput,
    const Matrix& rResidualDerivatives,
    const GeometryType& rGeometry) const
{
    if (rOutput.size1() != rResidualDerivatives.size1() ||
        rOutput.size2() != rResidualDerivatives.size2()) {
        rOutput.resize(rResidualDerivatives.size1(), rResidualDerivatives.size2(), false);
    }

    rOutput.clear();

    const IndexType number_of_nodes = rGeometry.PointsNumber();

    // add residual derivative contributions
    for (IndexType a = 0; a < number_of_nodes; ++a) {
        const auto& r_node = rGeometry[a];
        const IndexType block_index = a * mBlockSize;
        if (r_node.Is(SLIP)) {
            AddNodalRotationDerivatives(rOutput, rResidualDerivatives, block_index, r_node);
            // since slip condition is only based on first derivative
            // variable, make the column zero for all derivatives
            if (!r_node.IsFixed(ADJOINT_FLUID_VECTOR_1_X)) {
                // in here we have to check whether SLIP and not INLET.
                // At the inlet ADJOINT_FLUID_VECTOR_1_X is fixed.
                // There can be a situation where a node is
                // shared among both inlet and slip (corner node). In this case
                // AddNodalApplySlipConditionDerivatives will set all the values in the
                // column corresponding to ADJOINT_FLUID_VECTOR_1_X of the common node to zero except for the
                // nodal normal derivative contributions. If the normal is perpendicular to the X direction, then
                // this will set the diagonal in the column corresponding to this node's ADJOINT_FLUID_VECTOR_1_X dof
                // to zero. But there can be non-zero off-diagonal due to normal not being zero. Eventhough
                // ADJOINT_FLUID_VECTOR_1_X will be fixed in the block builder and solver
                // it will not recognize this column as empty, so diagonal won't be fixed by B&S to 1.0.
                // Since ADJOINT_FLUID_VECTOR_1_X of this node is fixed, B&S will set all the off-diagonals
                // to zero. This will create a zero column in the system matrix which makes it singular.
                // So this check is performed here to avoid this.
                // In here inlets are alsways assumed to be dirichlet.
                ClearNodalResidualDerivatives(rOutput, block_index);
            }
        } else {
            AddNodalResidualDerivatives(rOutput, rResidualDerivatives, block_index);
        }
    }
}

void FluidAdjointSlipUtilities::AddNodalApplySlipConditionDerivatives(
    Matrix& rOutput,
    const IndexType NodeStartIndex,
    const NodeType& rNode) const
{
    KRATOS_TRY

    // Apply slip condition in primal scheme makes first momentum dof
    // fixed, making the velocity in the normal direction as rhs.

    // first clear the residual derivative
    ClearNodalResidualDerivatives(rOutput, NodeStartIndex);

    auto normal = rNode.FastGetSolutionStepValue(NORMAL);
    normal /= norm_2(normal);

    for (IndexType i = 0; i < mDimension; ++i) {
        rOutput(NodeStartIndex + i, NodeStartIndex) -= normal[i];
    }

    KRATOS_CATCH("");
}

void FluidAdjointSlipUtilities::AddNodalResidualDerivatives(
    Matrix& rOutput,
    const Matrix& rResidualDerivatives,
    const IndexType NodeStartIndex) const
{
    KRATOS_TRY

    // add non-rotated residual derivative contributions
    for (IndexType c = 0; c < rResidualDerivatives.size1(); ++c) {
        for (IndexType i = 0; i < mBlockSize; ++i) {
            rOutput(c, NodeStartIndex + i) += rResidualDerivatives(c, NodeStartIndex + i);
        }
    }

    KRATOS_CATCH("");
}

void FluidAdjointSlipUtilities::ClearNodalResidualDerivatives(
    Matrix& rOutput,
    const IndexType ResidualIndex) const
{
    for (IndexType c = 0; c < rOutput.size1(); ++c) {
        rOutput(c, ResidualIndex) = 0.0;
    }
}

} // namespace Kratos