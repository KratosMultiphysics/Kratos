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

#if !defined(KRATOS_FLUID_ADJOINT_SLIP_UTILITIES_H)
#define KRATOS_FLUID_ADJOINT_SLIP_UTILITIES_H

// System includes
#include <unordered_map>
#include <vector>

// External includes

// Project includes
#include "includes/model_part.h"
#include "includes/ublas_interface.h"
#include "utilities/coordinate_transformation_utilities.h"

// Application includes
#include "custom_utilities/fluid_calculation_utilities.h"


namespace Kratos
{
///@addtogroup FluidDynamicsApplication
///@{

///@name Kratos classes
///@{

class KRATOS_API(FLUID_DYNAMICS_APPLICATION) FluidAdjointSlipUtilities
{
public:
    ///@name Type Definitions
    ///@{

    using IndexType = std::size_t;

    using NodeType = ModelPart::NodeType;

    using GeometryType = ModelPart::GeometryType;

    using ConditionType = ModelPart::ConditionType;

    using ElementType = ModelPart::ElementType;

    using CoordinateTransformationUtilities = CoordinateTransformationUtils<Matrix, Vector, double>;

    ///@}
    ///@name Life Cycle
    ///@{

    explicit FluidAdjointSlipUtilities(
        const IndexType Dimension,
        const IndexType BlockSize);

    explicit FluidAdjointSlipUtilities(
        const IndexType Dimension,
        const IndexType BlockSize,
        const Variable<double>& rSlipDofVariable);

    ///@}
    ///@name Operations
    ///@{

    void Initialize(ModelPart& rModelPart);

    /**
     * @brief Calculates rotated slip applied state derivatives
     *
     * This method calculates rotated slip applied state derivatives
     * when non-rotated residual derivatives are given. Nodal rotations
     * are determined by SLIP flag. If SLIP flag is true, then it is assumed
     * to be rotated as well as slip condition applied.
     *
     * This method assumes first dofs to be VELOCITY, and derivatives starting with
     * w.r.t. VELOCITY
     *
     * @param rOutput                   Rotated and slip applied state derivatives
     * @param rResidualDerivatives      Non-rotated and non-slip applied state derivatives
     * @param rGeometry                 Geometry of which rResidualDerivatives are computed on
     */
    void CalculateRotatedSlipConditionAppliedSlipVariableDerivatives(
        Matrix& rOutput,
        const Matrix& rResidualDerivatives,
        const GeometryType& rGeometry) const;

    /**
     * @brief Calculates rotated slip applied shape derivatives
     *
     * This method calculates rotated slip applied shape derivatives
     * when non-rotated residuals and its derivatives are given. Nodal rotations
     * are determined by SLIP flag. If SLIP flag is true, then it is assumed
     * to be rotated as well as slip condition applied.
     *
     * This method assumes first dofs to be VELOCITY
     *
     * @param rOutput                   Rotated and slip applied shape derivatives
     * @param rGPSensitivityVector      Global pointer vector order of the derivatives given in rOutput
     * @param rResiduals                Non-rotated and non-slip applied residuals
     * @param rResidualDerivatives      Non-rotated and non-slip applied shape derivatives
     * @param TEntityType               Entity of which rResidualDerivatives are computed on
     */
    template<class TEntityType>
    void CalculateRotatedSlipConditionAppliedShapeVariableDerivatives(
        Matrix& rOutput,
        GlobalPointersVector<NodeType>& rGPSensitivityVector,
        const Vector& rResiduals,
        const Matrix& rResidualDerivatives,
        const TEntityType& rEntity,
        const ProcessInfo& rProcessInfo) const;

    /**
     * @brief Calculates rotated slip applied shape derivatives
     *
     * This method calculates rotated slip applied shape derivatives
     * when non-rotated residuals and its derivatives are given. Nodal rotations
     * are determined by SLIP flag. If SLIP flag is true, then it is assumed
     * to be rotated as well as slip condition applied.
     *
     * This method assumes first dofs to be VELOCITY
     *
     * @param rOutput                   Rotated and slip applied shape derivatives
     * @param rNodeIds                  Node ids vector order of the derivatives given in rOutput
     * @param rResiduals                Non-rotated and non-slip applied residuals
     * @param rResidualDerivatives      Non-rotated and non-slip applied shape derivatives
     * @param TEntityType               Entity of which rResidualDerivatives are computed on
     */
    template<class TEntityType>
    void CalculateRotatedSlipConditionAppliedShapeVariableDerivatives(
        Matrix& rOutput,
        std::vector<IndexType>& rNodeIds,
        const Vector& rResiduals,
        const Matrix& rResidualDerivatives,
        const TEntityType& rEntity,
        const ProcessInfo& rProcessInfo) const;

    /**
     * @brief Calculates rotated slip applied state derivatives
     *
     * This method calculates rotated slip applied state derivatives
     * when non-rotated residual derivatives are given. Nodal rotations
     * are determined by SLIP flag. If SLIP flag is true, then it is assumed
     * to be rotated as well as slip condition applied.
     *
     * This method assumes first dofs to be VELOCITY, and derivatives not starting with
     * w.r.t. VELOCITY because SLIP is based on the VELOCITY variable.
     *
     * @param rOutput                   Rotated and slip applied state derivatives
     * @param rResidualDerivatives      Non-rotated and non-slip applied state derivatives
     * @param rGeometry                 Geometry of which rResidualDerivatives are computed on
     */
    void CalculateRotatedSlipConditionAppliedNonSlipNonShapeVariableDerivatives(
        Matrix& rOutput,
        const Matrix& rResidualDerivatives,
        const GeometryType& rGeometry) const;

    /**
     * @brief Calculates rotated residual derivatives w.r.t shape variables
     *
     * This method calculates rotated residuals shape derivatives w.r.t. shape variables.
     *
     * @param rOutput                   Nodal rotated shape derivatives
     * @param rResiduals                Non-rotated residuals
     * @param rResidualDerivatives      Non-rotated shape derivatives
     * @param NodeStartIndex            Block starting column for that specific node
     * @param rDerivativesMap           <Node id, derivative node position> map for shape derivatives
     * @param rNode                     Node which has rotated residual derivative contributions
     */
    void AddNodalRotationShapeVariableDerivatives(
        Matrix& rOutput,
        const Vector& rResiduals,
        const Matrix& rResidualDerivatives,
        const IndexType NodeStartIndex,
        const std::unordered_map<IndexType, IndexType>& rDerivativesMap,
        const NodeType& rNode) const
    {
        (this->*(this->mAddNodalRotationShapeVariableDerivativesMethod))(rOutput, rResiduals, rResidualDerivatives, NodeStartIndex, rDerivativesMap, rNode);
    }

    /**
     * @brief Rotates residual derivatives
     *
     * This method is used to rotate part of the non-rotated residual derivatives
     * which corresponds to a given node.
     *
     * @param rOutput                   Nodal rotated state derivatives
     * @param rResidualDerivatives      Non-rotated state derivatives
     * @param NodeStartIndex            Block starting column for that specific node
     * @param rNode                     Node which has rotated residual derivative contributions
     */
    void AddNodalRotationNonShapeVariableDerivatives(
        Matrix& rOutput,
        const Matrix& rResidualDerivatives,
        const IndexType NodeStartIndex,
        const NodeType& rNode) const
    {
        (this->*(this->mAddNodalRotationNonShapeVariableDerivativesMethod))(rOutput, rResidualDerivatives, NodeStartIndex, rNode);
    }

    /**
     * @brief Slip condition residual derivatives
     *
     * This method is used to apply slip condition state derivatives to part of the non-rotated residual derivatives
     * which corresponds to a given node.
     *
     * @param rOutput                   Nodal slip condition applied state derivatives
     * @param NodeStartIndex            Block starting column for that specific node
     * @param rNode                     Node which has rotated residual derivative contributions
     */
    void AddNodalApplySlipConditionSlipVariableDerivatives(
        Matrix& rOutput,
        const IndexType NodeStartIndex,
        const NodeType& rNode) const;

    /**
     * @brief Slip condition residual derivatives w.r.t. shape
     *
     * This method is used to apply slip condition state derivatives to part of the non-rotated residual derivatives
     * which corresponds to a given node w.r.t shape.
     *
     * @param rOutput                   Nodal slip condition applied state derivatives
     * @param NodeStartIndex            Block starting column for that specific node
     * @param rDerivativesMap           <Node id, derivative node position> map for shape derivatives
     * @param rNode                     Node which has rotated residual derivative contributions
     */
    void AddNodalApplySlipConditionShapeVariableDerivatives(
        Matrix& rOutput,
        const IndexType NodeStartIndex,
        const std::unordered_map<IndexType, IndexType>& rDerivativesMap,
        const NodeType& rNode) const
    {
        (this->*(this->mAddNodalApplySlipConditionShapeVariableDerivativesMethod))(rOutput, NodeStartIndex, rDerivativesMap, rNode);
    }

    /**
     * @brief Adds the nodal state derivatives as it is
     *
     * This method is used to apply slip condition state derivatives to part of the non-rotated residual derivatives
     * which corresponds to a given node.
     *
     * @param rOutput                   Output state derivatives
     * @param rResidualDerivatives      Non-rotated state derivatives
     * @param NodeStartIndex            Block starting column for that specific node
     */
    void AddNodalResidualDerivatives(
        Matrix& rOutput,
        const Matrix& rResidualDerivatives,
        const IndexType NodeStartIndex) const;

    /**
     * @brief Clears a residual derivative
     *
     * @param rOutput               Output matrix
     * @param ResidualIndex         Residual index which needs to be set to zero.
     */
    void ClearNodalResidualDerivatives(
        Matrix& rOutput,
        const IndexType ResidualIndex) const;

    ///@}
private:
    ///@name Private Members
    ///@{

    const IndexType mDimension;
    const IndexType mBlockSize;

    const CoordinateTransformationUtilities mRotationTool;

    void (FluidAdjointSlipUtilities::*mAddNodalRotationNonShapeVariableDerivativesMethod)(
        Matrix&,
        const Matrix&,
        const IndexType,
        const NodeType&) const;

    void (FluidAdjointSlipUtilities::*mAddNodalRotationShapeVariableDerivativesMethod)(
        Matrix&,
        const Vector&,
        const Matrix&,
        const IndexType,
        const std::unordered_map<IndexType, IndexType>&,
        const NodeType&) const;

    void (FluidAdjointSlipUtilities::*mAddNodalApplySlipConditionShapeVariableDerivativesMethod)(
        Matrix&,
        const IndexType,
        const std::unordered_map<IndexType, IndexType>&,
        const NodeType&) const;

    const Variable<double>& mrSlipDofVariable;

    std::unordered_map<int, std::vector<int>> mNodalNeighboursMap;

    std::unordered_map<int, GlobalPointer<ModelPart::NodeType>> mGlobalPointerNodalMap;

    ///@}
    ///@name Private Operations
    ///@{

    template<unsigned int TDim>
    void TemplatedAddNodalRotationNonShapeVariableDerivatives(
        Matrix& rOutput,
        const Matrix& rResidualDerivatives,
        const IndexType NodeStartIndex,
        const NodeType& rNode) const
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
            for (IndexType a = TDim; a < mBlockSize; ++a) {
                rOutput(c, NodeStartIndex + a) += rResidualDerivatives(c, NodeStartIndex + a);
            }
        }

        KRATOS_CATCH("");
    }

    template<unsigned int TDim>
    void TemplatedAddNodalRotationShapeVariableDerivatives(
        Matrix& rOutput,
        const Vector& rResiduals,
        const Matrix& rResidualDerivatives,
        const IndexType NodeStartIndex,
        const std::unordered_map<IndexType, IndexType>& rDerivativesMap,
        const NodeType& rNode) const
    {
        KRATOS_TRY

        TemplatedAddNodalRotationNonShapeVariableDerivatives<TDim>(rOutput, rResidualDerivatives, NodeStartIndex, rNode);

        // // get the residual relevant for rNode
        BoundedVector<double, TDim> residual, aux_vector;
        FluidCalculationUtilities::ReadSubVector<TDim>(residual, rResiduals, NodeStartIndex);

        // first add rotation matrix derivative contributions w.r.t. rNode
        BoundedMatrix<double, TDim, TDim> rotation_matrix_derivative;
        const int current_node_index = rDerivativesMap.find(rNode.Id())->second * TDim;
        for (IndexType k = 0; k < TDim; ++k) {
            mRotationTool.CalculateRotationOperatorPureShapeSensitivities(
                rotation_matrix_derivative, 0, k, rNode);

            noalias(aux_vector) = prod(rotation_matrix_derivative, residual);
            FluidCalculationUtilities::AddSubVector<TDim>(
                rOutput, aux_vector, current_node_index + k, NodeStartIndex);
        }

        // add rotation matrix derivative contributions w.r.t. rNode neighbors
        const auto& r_neighbour_ids = mNodalNeighboursMap.find(rNode.Id())->second;
        for (IndexType b = 0; b < r_neighbour_ids.size(); ++b) {
            const int derivative_node_index =
                rDerivativesMap.find(r_neighbour_ids[b])->second * TDim;
            for (IndexType k = 0; k < TDim; ++k) {
                mRotationTool.CalculateRotationOperatorPureShapeSensitivities(
                    rotation_matrix_derivative, b + 1, k, rNode);

                noalias(aux_vector) = prod(rotation_matrix_derivative, residual);
                FluidCalculationUtilities::AddSubVector<TDim>(
                    rOutput, aux_vector, derivative_node_index + k, NodeStartIndex);
            }
        }

        KRATOS_CATCH("");
    }

    template<unsigned int TDim>
    void TemplatedAddNodalApplySlipConditionShapeVariableDerivatives(
        Matrix& rOutput,
        const IndexType NodeStartIndex,
        const std::unordered_map<IndexType, IndexType>& rDerivativesMap,
        const NodeType& rNode) const
    {
        KRATOS_TRY

        // Apply slip condition in primal scheme makes first momentum dof
        // fixed, making the velocity in the normal direction as rhs.

        // first clear the residual derivative
        for (IndexType c = 0; c < rOutput.size1(); ++c) {
            rOutput(c, NodeStartIndex) = 0.0;
        }

        array_1d<double, TDim> effective_velocity, normal;
        for (IndexType i = 0; i < TDim; ++i) {
            effective_velocity[i] = rNode.FastGetSolutionStepValue(MESH_VELOCITY)[i];
            effective_velocity[i] -= rNode.FastGetSolutionStepValue(VELOCITY)[i];
            normal[i] = rNode.FastGetSolutionStepValue(NORMAL)[i];
        }
        const double normal_magnitude = norm_2(normal);
        const Matrix& normal_shape_derivatives = rNode.GetValue(NORMAL_SHAPE_DERIVATIVE);

        // add unit normal derivative w.r.t. current node
        const int current_node_index = rDerivativesMap.find(rNode.Id())->second * TDim;
        for (IndexType k = 0; k < TDim; ++k) {
            const auto& nodal_normal_derivative = row(normal_shape_derivatives, k);

            double value = 0.0;
            value += inner_prod(nodal_normal_derivative, effective_velocity) / normal_magnitude;
            value -= inner_prod(normal, effective_velocity) *
                     inner_prod(normal, nodal_normal_derivative) /
                     std::pow(normal_magnitude, 3);

            rOutput(current_node_index + k, NodeStartIndex) = value;
        }

        // add unit normal derivative w.r.t. neighbour nodes
        const auto& r_neighbour_ids = mNodalNeighboursMap.find(rNode.Id())->second;
        for (IndexType b = 0; b < r_neighbour_ids.size(); ++b) {
            const int derivative_node_index =
                rDerivativesMap.find(r_neighbour_ids[b])->second * TDim;
            for (IndexType k = 0; k < TDim; ++k) {
                const auto& nodal_normal_derivative =
                    row(normal_shape_derivatives, (b + 1) * TDim + k);

                double value = 0.0;
                value += inner_prod(nodal_normal_derivative, effective_velocity) / normal_magnitude;
                value -= inner_prod(normal, effective_velocity) *
                         inner_prod(normal, nodal_normal_derivative) /
                         std::pow(normal_magnitude, 3);

                rOutput(derivative_node_index + k, NodeStartIndex) = value;
            }
        }

        KRATOS_CATCH("");
    }

    ///@}
};

///@}

///@}

} // namespace Kratos

#endif // KRATOS_FLUID_ADJOINT_SLIP_UTILITIES_H
