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

// Include base h
#include "fluid_adjoint_utilities.h"

namespace Kratos
{
template <unsigned int TDim, unsigned int TBlockSize>
const CoordinateTransformationUtils<Matrix, Vector, double> FluidAdjointUtilities<TDim, TBlockSize>::mRotationTool(
    TDim, TBlockSize);

template <unsigned int TDim, unsigned int TBlockSize>
void FluidAdjointUtilities<TDim, TBlockSize>::CalculateRotatedSlipConditionAppliedSlipVariableDerivatives(
    Matrix& rOutput,
    const Matrix& rResidualDerivatives,
    const GeometryType& rGeometry)
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
        const IndexType block_index = a * TBlockSize;
        if (r_node.Is(SLIP)) {
            AddNodalRotationDerivatives(rOutput, rResidualDerivatives, block_index, r_node);
            AddNodalApplySlipConditionDerivatives(rOutput, block_index, r_node);
        } else {
            AddNodalResidualDerivatives(rOutput, rResidualDerivatives, block_index);
        }
    }
}

template <unsigned int TDim, unsigned int TBlockSize>
void FluidAdjointUtilities<TDim, TBlockSize>::CalculateRotatedSlipConditionAppliedNonSlipVariableDerivatives(
    Matrix& rOutput,
    const Matrix& rResidualDerivatives,
    const GeometryType& rGeometry)
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
        const IndexType block_index = a * TBlockSize;
        if (r_node.Is(SLIP)) {
            AddNodalRotationDerivatives(rOutput, rResidualDerivatives, block_index, r_node);
            // since slip condition is only based on first derivative
            // variable, make the column zero for all derivatives
            ClearNodalResidualDerivatives(rOutput, block_index);
        } else {
            AddNodalResidualDerivatives(rOutput, rResidualDerivatives, block_index);
        }
    }
}

template <unsigned int TDim, unsigned int TBlockSize>
void FluidAdjointUtilities<TDim, TBlockSize>::AddNodalRotationDerivatives(
    Matrix& rOutput,
    const Matrix& rResidualDerivatives,
    const IndexType NodeStartIndex,
    const NodeType& rNode)
{
    KRATOS_TRY

    BoundedVector<double, TDim> residual_derivative, aux_vector;
    BoundedMatrix<double, TDim, TDim> rotation_matrix;
    mRotationTool.LocalRotationOperatorPure(rotation_matrix, rNode);

    // add rotated residual derivative contributions
    for (IndexType c = 0; c < rResidualDerivatives.size1(); ++c) {
        // get the residual derivative relevant for node
        FluidCalculationUtilities::ReadSubVector<TDim>(
            residual_derivative, row(rResidualDerivatives, c), NodeStartIndex);

        // rotate residual derivative
        noalias(aux_vector) = prod(rotation_matrix, residual_derivative);

        // add rotated residual derivative to local matrix
        FluidCalculationUtilities::AddSubVector<TDim>(rOutput, aux_vector, c, NodeStartIndex);

        // add rest of the equation derivatives
        for (IndexType a = TDim; a < TBlockSize; ++a) {
            rOutput(c, NodeStartIndex + a) += rResidualDerivatives(c, NodeStartIndex + a);
        }
    }

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TBlockSize>
void FluidAdjointUtilities<TDim, TBlockSize>::AddNodalApplySlipConditionDerivatives(
    Matrix& rOutput,
    const IndexType NodeStartIndex,
    const NodeType& rNode)
{
    KRATOS_TRY

    // Apply slip condition in primal scheme makes first momentum dof
    // fixed, making the velocity in the normal direction as rhs.

    // first clear the residual derivative
    ClearNodalResidualDerivatives(rOutput, NodeStartIndex);

    auto normal = rNode.FastGetSolutionStepValue(NORMAL);
    normal /= norm_2(normal);

    for (IndexType i = 0; i < TDim; ++i) {
        rOutput(NodeStartIndex + i, NodeStartIndex) -= normal[i];
    }

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TBlockSize>
void FluidAdjointUtilities<TDim, TBlockSize>::AddNodalResidualDerivatives(
    Matrix& rOutput,
    const Matrix& rResidualDerivatives,
    const IndexType NodeStartIndex)
{
    KRATOS_TRY

    // add non-rotated residual derivative contributions
    for (IndexType c = 0; c < rResidualDerivatives.size1(); ++c) {
        for (IndexType i = 0; i < TBlockSize; ++i) {
            rOutput(c, NodeStartIndex + i) += rResidualDerivatives(c, NodeStartIndex + i);
        }
    }

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TBlockSize>
void FluidAdjointUtilities<TDim, TBlockSize>::ClearNodalResidualDerivatives(
    Matrix& rOutput,
    const IndexType ResidualIndex)
{
    for (IndexType c = 0; c < rOutput.size1(); ++c) {
        rOutput(c, ResidualIndex) = 0.0;
    }
}

// template instantiations

template class FluidAdjointUtilities<2, 3>;
template class FluidAdjointUtilities<2, 4>;
template class FluidAdjointUtilities<2, 5>;

template class FluidAdjointUtilities<3, 4>;
template class FluidAdjointUtilities<3, 5>;
template class FluidAdjointUtilities<3, 6>;

} // namespace Kratos
