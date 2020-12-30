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

#if !defined(KRATOS_SIMPLE_STEADY_SENSITIVITY_BUILDER_SCHEME_H_INCLUDED)
#define KRATOS_SIMPLE_STEADY_SENSITIVITY_BUILDER_SCHEME_H_INCLUDED

// System includes
#include <unordered_map>

// External includes

// Project includes
#include "containers/global_pointers_unordered_map.h"
#include "containers/global_pointers_vector.h"
#include "includes/define.h"
#include "includes/global_pointer_variables.h"
#include "includes/model_part.h"
#include "processes/find_global_nodal_neighbours_for_entities_process.h"
#include "response_functions/adjoint_response_function.h"
#include "solving_strategies/schemes/sensitivity_builder_scheme.h"
#include "utilities/coordinate_transformation_utilities.h"
#include "utilities/normal_calculation_utils.h"
#include "utilities/openmp_utils.h"
#include "utilities/sensitivity_builder.h"

// Application includes
#include "custom_utilities/fluid_calculation_utilities.h"

namespace Kratos
{
///@name Kratos Classes
///@{

template<IndexType TDim>
class SimpleSteadySensitivityBuilderScheme : public SensitivityBuilderScheme
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(SimpleSteadySensitivityBuilderScheme);

    using BaseType = SensitivityBuilderScheme;

    using NodeType = BaseType::NodeType;

    using ConditionType = BaseType::ConditionType;

    using ElementType = BaseType::ElementType;

    using IndexType = std::size_t;

    static constexpr IndexType TBlockSize = TDim + 1;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    SimpleSteadySensitivityBuilderScheme()
        : SensitivityBuilderScheme()
    {
        // Allocate auxiliary memory.
        // This needs to be done in the constructor because, this scheme
        // is used to calculate sensitivities w.r.t. element quantities
        // and in there we don't usually pass the model part. Hence, we
        // can not call SimpleSteadySensitivityBuilderScheme::Initialize
        // method.
        const int number_of_threads = OpenMPUtils::GetNumThreads();
        mAuxVectors.resize(number_of_threads);
        mAuxMatrices.resize(number_of_threads);
        mRotatedSensitivityMatrices.resize(number_of_threads);
    }

    /// Destructor.
    ~SimpleSteadySensitivityBuilderScheme() = default;

    ///@}
    ///@name Operations
    ///@{

    void InitializeSolutionStep(
        ModelPart& rModelPart,
        ModelPart& rSensitivityModelPart,
        AdjointResponseFunction& rResponseFunction) override
    {
        KRATOS_TRY

        if (!mIsNodalNormalShapeDerivativesComputed) {
            mIsNodalNormalShapeDerivativesComputed = true;

            // Find nodal neighbors for conditions
            FindNodalNeighboursForEntitiesProcess<ModelPart::ConditionsContainerType> neighbor_process(
                rModelPart.GetCommunicator().GetDataCommunicator(), rModelPart,
                NEIGHBOUR_CONDITION_NODES);
            neighbor_process.Execute();
            mNodalNeighboursMap = neighbor_process.GetNeighbourIds(rModelPart.Nodes());

            // Compute condition NORMAL derivatives and store them in conditions
            NormalCalculationUtils().CalculateNormalShapeDerivativesOnSimplex(
                rModelPart.Conditions(), TDim);

            // Assign condition normal derivatives to corresponding nodes
            SensitivityBuilder::AssignEntityDerivativesToNodes<ModelPart::ConditionsContainerType>(
                rModelPart, TDim, NORMAL_SHAPE_DERIVATIVE,
                mNodalNeighboursMap, 1.0 / TDim, SLIP);
        }

        BaseType::InitializeSolutionStep(rModelPart, rSensitivityModelPart, rResponseFunction);

        KRATOS_CATCH("");
    }

    /**
     * @brief Calculates sensitivity from a given element
     *
     * This method needs to be implemented in the derrived class if sensitivity calculation is
     * based on nodal quantities, and elements has a contribution towards that.
     *
     * rSensitivity vector should have sensitivities in the order maching rGPSensitivity vector.
     * rGPSensitivityVector holds GlobalPointers indicating where these sensitivities need to
     * be stored accordingly.
     *
     * @param rCurrentElement           Element being passed
     * @param rResponseFunction         Response function used
     * @param rSensitivity              Sensitivity vector
     * @param rGPSensitivityVector      Global pointer vector
     * @param rVariable                 Design variable
     * @param rCurrentProcessInfo       Current process info
     */
    void CalculateSensitivity(
        ElementType& rCurrentElement,
        AdjointResponseFunction& rResponseFunction,
        Vector& rSensitivity,
        GlobalPointersVector<NodeType>& rGPSensitivityVector,
        const Variable<double>& rVariable,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        CalculateLocalSensitivityAndGlobalPointersVector(
            rCurrentElement, rResponseFunction, rSensitivity,
            rGPSensitivityVector, rVariable, rCurrentProcessInfo);
    }

    /**
     * @brief Calculates sensitivity from a given element
     *
     * This method needs to be implemented in the derrived class if sensitivity calculation is
     * based on element quantities, and elements has a contribution towards that.
     *
     * rSensitivity vector should have sensitivities in the order maching rGPSensitivity vector.
     * rGPSensitivityVector holds GlobalPointers indicating where these sensitivities need to
     * be stored accordingly.
     *
     * @param rCurrentElement           Element being passed
     * @param rResponseFunction         Response function used
     * @param rSensitivity              Sensitivity vector
     * @param rGPSensitivityVector      Global pointer vector
     * @param rVariable                 Design variable
     * @param rCurrentProcessInfo       Current process info
     */
    void CalculateSensitivity(
        ElementType& rCurrentElement,
        AdjointResponseFunction& rResponseFunction,
        Vector& rSensitivity,
        GlobalPointersVector<ElementType>& rGPSensitivityVector,
        const Variable<double>& rVariable,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        CalculateLocalSensitivityAndGlobalPointersVector(
            rCurrentElement, rResponseFunction, rSensitivity,
            rGPSensitivityVector, rVariable, rCurrentProcessInfo);
    }

    /**
     * @brief Calculates sensitivity from a given condition
     *
     * This method needs to be implemented in the derrived class if sensitivity calculation is
     * based on nodal quantities, and conditions has a contribution towards that.
     *
     * rSensitivity vector should have sensitivities in the order maching rGPSensitivity vector.
     * rGPSensitivityVector holds GlobalPointers indicating where these sensitivities need to
     * be stored accordingly.
     *
     * @param rCurrentCondition         Condition being passed
     * @param rResponseFunction         Response function used
     * @param rSensitivity              Sensitivity vector
     * @param rGPSensitivityVector      Global pointer vector
     * @param rVariable                 Design variable
     * @param rCurrentProcessInfo       Current process info
     */
    void CalculateSensitivity(
        ConditionType& rCurrentCondition,
        AdjointResponseFunction& rResponseFunction,
        Vector& rSensitivity,
        GlobalPointersVector<NodeType>& rGPSensitivityVector,
        const Variable<double>& rVariable,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        CalculateLocalSensitivityAndGlobalPointersVector(
            rCurrentCondition, rResponseFunction, rSensitivity,
            rGPSensitivityVector, rVariable, rCurrentProcessInfo);
    }

    /**
     * @brief Calculates sensitivity from a given condition
     *
     * This method needs to be implemented in the derrived class if sensitivity calculation is
     * based on condition quantities, and conditions has a contribution towards that.
     *
     * rSensitivity vector should have sensitivities in the order maching rGPSensitivity vector.
     * rGPSensitivityVector holds GlobalPointers indicating where these sensitivities need to
     * be stored accordingly.
     *
     * @param rCurrentCondition         Condition being passed
     * @param rResponseFunction         Response function used
     * @param rSensitivity              Sensitivity vector
     * @param rGPSensitivityVector      Global pointer vector
     * @param rVariable                 Design variable
     * @param rCurrentProcessInfo       Current process info
     */
    void CalculateSensitivity(
        ConditionType& rCurrentCondition,
        AdjointResponseFunction& rResponseFunction,
        Vector& rSensitivity,
        GlobalPointersVector<ConditionType>& rGPSensitivityVector,
        const Variable<double>& rVariable,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        CalculateLocalSensitivityAndGlobalPointersVector(
            rCurrentCondition, rResponseFunction, rSensitivity,
            rGPSensitivityVector, rVariable, rCurrentProcessInfo);
    }

    /**
     * @brief Calculates sensitivity from a given element
     *
     * This method needs to be implemented in the derrived class if sensitivity calculation is
     * based on nodal quantities, and elements has a contribution towards that.
     *
     * rSensitivity vector should have sensitivities in the order maching rGPSensitivity vector.
     * rGPSensitivityVector holds GlobalPointers indicating where these sensitivities need to
     * be stored accordingly.
     *
     * @param rCurrentElement           Element being passed
     * @param rResponseFunction         Response function used
     * @param rSensitivity              Sensitivity vector
     * @param rGPSensitivityVector      Global pointer vector
     * @param rVariable                 Design variable
     * @param rCurrentProcessInfo       Current process info
     */
    void CalculateSensitivity(
        ElementType& rCurrentElement,
        AdjointResponseFunction& rResponseFunction,
        Vector& rSensitivity,
        GlobalPointersVector<NodeType>& rGPSensitivityVector,
        const Variable<array_1d<double, 3>>& rVariable,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        CalculateLocalSensitivityAndGlobalPointersVector(
            rCurrentElement, rResponseFunction, rSensitivity,
            rGPSensitivityVector, rVariable, rCurrentProcessInfo);
    }

    /**
     * @brief Calculates sensitivity from a given element
     *
     * This method needs to be implemented in the derrived class if sensitivity calculation is
     * based on element quantities, and elements has a contribution towards that.
     *
     * rSensitivity vector should have sensitivities in the order maching rGPSensitivity vector.
     * rGPSensitivityVector holds GlobalPointers indicating where these sensitivities need to
     * be stored accordingly.
     *
     * @param rCurrentElement           Element being passed
     * @param rResponseFunction         Response function used
     * @param rSensitivity              Sensitivity vector
     * @param rGPSensitivityVector      Global pointer vector
     * @param rVariable                 Design variable
     * @param rCurrentProcessInfo       Current process info
     */
    void CalculateSensitivity(
        ElementType& rCurrentElement,
        AdjointResponseFunction& rResponseFunction,
        Vector& rSensitivity,
        GlobalPointersVector<ElementType>& rGPSensitivityVector,
        const Variable<array_1d<double, 3>>& rVariable,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        CalculateLocalSensitivityAndGlobalPointersVector(
            rCurrentElement, rResponseFunction, rSensitivity,
            rGPSensitivityVector, rVariable, rCurrentProcessInfo);
    }

    /**
     * @brief Calculates sensitivity from a given condition
     *
     * This method needs to be implemented in the derrived class if sensitivity calculation is
     * based on nodal quantities, and conditions has a contribution towards that.
     *
     * rSensitivity vector should have sensitivities in the order maching rGPSensitivity vector.
     * rGPSensitivityVector holds GlobalPointers indicating where these sensitivities need to
     * be stored accordingly.
     *
     * @param rCurrentCondition         Condition being passed
     * @param rResponseFunction         Response function used
     * @param rSensitivity              Sensitivity vector
     * @param rGPSensitivityVector      Global pointer vector
     * @param rVariable                 Design variable
     * @param rCurrentProcessInfo       Current process info
     */
    void CalculateSensitivity(
        ConditionType& rCurrentCondition,
        AdjointResponseFunction& rResponseFunction,
        Vector& rSensitivity,
        GlobalPointersVector<NodeType>& rGPSensitivityVector,
        const Variable<array_1d<double, 3>>& rVariable,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        CalculateLocalSensitivityAndGlobalPointersVector(
            rCurrentCondition, rResponseFunction, rSensitivity,
            rGPSensitivityVector, rVariable, rCurrentProcessInfo);
    }

    /**
     * @brief Calculates sensitivity from a given condition
     *
     * This method needs to be implemented in the derrived class if sensitivity calculation is
     * based on condition quantities, and conditions has a contribution towards that.
     *
     * rSensitivity vector should have sensitivities in the order maching rGPSensitivity vector.
     * rGPSensitivityVector holds GlobalPointers indicating where these sensitivities need to
     * be stored accordingly.
     *
     * @param rCurrentCondition         Condition being passed
     * @param rResponseFunction         Response function used
     * @param rSensitivity              Sensitivity vector
     * @param rGPSensitivityVector      Global pointer vector
     * @param rVariable                 Design variable
     * @param rCurrentProcessInfo       Current process info
     */
    void CalculateSensitivity(
        ConditionType& rCurrentCondition,
        AdjointResponseFunction& rResponseFunction,
        Vector& rSensitivity,
        GlobalPointersVector<ConditionType>& rGPSensitivityVector,
        const Variable<array_1d<double, 3>>& rVariable,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        CalculateLocalSensitivityAndGlobalPointersVector(
            rCurrentCondition, rResponseFunction, rSensitivity,
            rGPSensitivityVector, rVariable, rCurrentProcessInfo);
    }

    void CalculateResidualSensitivityMatrix(
        ElementType& rElement,
        Vector& rAdjointValues,
        Matrix& rOutput,
        GlobalPointersVector<NodeType>& rGPSensitivityVector,
        const Variable<array_1d<double, 3>>& rVariable,
        const ProcessInfo& rCurrentProcessInfo)
    {
        CalculateResidualSensitivityMatrix<ElementType, array_1d<double, 3>>(
            rElement, rAdjointValues, rOutput, rGPSensitivityVector, rVariable,
            rCurrentProcessInfo);
    }

    void CalculateResidualSensitivityMatrix(
        ConditionType& rElement,
        Vector& rAdjointValues,
        Matrix& rOutput,
        GlobalPointersVector<NodeType>& rGPSensitivityVector,
        const Variable<array_1d<double, 3>>& rVariable,
        const ProcessInfo& rCurrentProcessInfo)
    {
        CalculateResidualSensitivityMatrix<ElementType, array_1d<double, 3>>(
            rElement, rAdjointValues, rOutput, rGPSensitivityVector, rVariable,
            rCurrentProcessInfo);
    }

    void Clear() override
    {
        BaseType::Clear();
        mAuxVectors.clear();
        mAuxMatrices.clear();
        mRotatedSensitivityMatrices.clear();
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "SimpleSteadySensitivityBuilderScheme";
    }

    ///@}

private:
    ///@name Member Variables
    ///@{

    bool mIsNodalNormalShapeDerivativesComputed = false;
    std::vector<Matrix> mAuxMatrices;
    std::vector<Vector> mAuxVectors;
    std::vector<Matrix> mRotatedSensitivityMatrices;

    std::unordered_map<int, std::vector<int>> mNodalNeighboursMap;

    ///@}
    ///@name Private Operations
    ///@{

    template <typename TEntityType, typename TDerivativeEntityType, typename TDataType>
    void CalculateLocalSensitivityAndGlobalPointersVector(
        TEntityType& rEntity,
        AdjointResponseFunction& rResponseFunction,
        Vector& rSensitivityVector,
        GlobalPointersVector<TDerivativeEntityType>& rGPSensitivityVector,
        const Variable<TDataType>& rVariable,
        const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        const auto k = OpenMPUtils::ThisThread();

        rEntity.CalculateSensitivityMatrix(rVariable, mSensitivityMatrices[k], rProcessInfo);
        rEntity.GetValuesVector(mAdjointVectors[k]);

        KRATOS_ERROR_IF(mAdjointVectors[k].size() != mSensitivityMatrices[k].size2())
            << "mAdjointVectors.size(): " << mAdjointVectors[k].size()
            << " incompatible with mSensitivityMatrices[k].size1(): "
            << mSensitivityMatrices[k].size2() << ". Variable: " << rVariable << std::endl;

        rResponseFunction.CalculatePartialSensitivity(
            rEntity, rVariable, mSensitivityMatrices[k], mPartialSensitivity[k], rProcessInfo);

        KRATOS_ERROR_IF(mPartialSensitivity[k].size() != mSensitivityMatrices[k].size1())
            << "mPartialSensitivity.size(): " << mPartialSensitivity[k].size()
            << " incompatible with mSensitivityMatrices.size1(): "
            << mSensitivityMatrices[k].size1() << ". Variable: " << rVariable << std::endl;

        if (rSensitivityVector.size() != mSensitivityMatrices[k].size1()) {
            rSensitivityVector.resize(mSensitivityMatrices[k].size1(), false);
        }

        noalias(rSensitivityVector) = prod(mSensitivityMatrices[k], mAdjointVectors[k]) + mPartialSensitivity[k];

        if (rGPSensitivityVector.size() != 1) {
            rGPSensitivityVector.resize(1);
        }

        rGPSensitivityVector(0) = GlobalPointer<TDerivativeEntityType>(&rEntity, mRank);

        KRATOS_CATCH("");
    }

    template <typename TEntityType, typename TDataType>
    void CalculateLocalSensitivityAndGlobalPointersVector(
        TEntityType& rEntity,
        AdjointResponseFunction& rResponseFunction,
        Vector& rSensitivityVector,
        GlobalPointersVector<NodeType>& rGPSensitivityVector,
        const Variable<TDataType>& rVariable,
        const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        const auto k = OpenMPUtils::ThisThread();

        auto& adjoint_vector = mAdjointVectors[k];
        auto& rotated_sensitivity_matrix = mRotatedSensitivityMatrices[k];

        // get adjoint solution vector
        rEntity.GetValuesVector(adjoint_vector);
        const IndexType residuals_size = adjoint_vector.size();

        if (residuals_size != 0) {
            this->CalculateResidualSensitivityMatrix<TEntityType, TDataType>(
                rEntity, adjoint_vector, rotated_sensitivity_matrix,
                rGPSensitivityVector, rVariable, rProcessInfo);

            if (rSensitivityVector.size() != rotated_sensitivity_matrix.size1()) {
                rSensitivityVector.resize(rotated_sensitivity_matrix.size1(), false);
            }
            noalias(rSensitivityVector) = prod(rotated_sensitivity_matrix, adjoint_vector);

            // add objective derivative contributions
            auto& objective_partial_sensitivity = mPartialSensitivity[k];
            rResponseFunction.CalculatePartialSensitivity(
                rEntity, rVariable, rotated_sensitivity_matrix,
                objective_partial_sensitivity, rProcessInfo);

            KRATOS_DEBUG_ERROR_IF(objective_partial_sensitivity.size() >
                                  rSensitivityVector.size())
                << "rSensitivityVector does not have sufficient rows to add "
                   "objective partial sensitivity. [ rSensitivityVector.size() "
                   "= "
                << rSensitivityVector.size() << ", PartialSensitivityVectorSize = "
                << objective_partial_sensitivity.size() << " ].\n";

            // this assumes objective_partial_sensitivity also follows the same order of gps given by
            // gp_index_map. This has to be done in this way because AdjointResponseFunction does not have
            // an interface to pass a gp_vector.
            // may be we can extend AdjointResponseFunction interface?
            for (IndexType c = 0; c < objective_partial_sensitivity.size(); ++c) {
                rSensitivityVector[c] += objective_partial_sensitivity[c];
            }
        }

        KRATOS_CATCH("");
    }

    template <typename TEntityType, typename TDataType>
    void CalculateResidualSensitivityMatrix(
        TEntityType& rEntity,
        Vector& rAdjointValues,
        Matrix& rOutput,
        GlobalPointersVector<NodeType>& rGPSensitivityVector,
        const Variable<TDataType>& rVariable,
        const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        const auto k = OpenMPUtils::ThisThread();

        auto& sensitivity_matrix = mSensitivityMatrices[k];
        auto& aux_vector = mAuxVectors[k];
        auto& aux_matrix = mAuxMatrices[k];

        // get adjoint solution vector
        rEntity.GetValuesVector(rAdjointValues);
        const IndexType residuals_size = rAdjointValues.size();

        // calculate entity residual derivatives
        rEntity.CalculateSensitivityMatrix(rVariable, sensitivity_matrix, rProcessInfo);

        KRATOS_ERROR_IF(residuals_size != sensitivity_matrix.size2())
            << "mAdjointVectors.size(): " << residuals_size
            << " incompatible with sensitivity_matrix.size1(): "
            << sensitivity_matrix.size2() << ". Variable: " << rVariable << std::endl;

        // a map to store node id as key and the node order as the value
        std::unordered_map<IndexType, IndexType> gp_index_map;
        auto& r_geometry = rEntity.GetGeometry();
        const IndexType number_of_nodes = r_geometry.PointsNumber();
        const IndexType local_size = residuals_size / number_of_nodes;

        if (rGPSensitivityVector.size() != number_of_nodes) {
            rGPSensitivityVector.resize(number_of_nodes);
        }

        // add geometry gps
        for (IndexType i = 0; i < number_of_nodes; ++i) {
            const auto& r_node = r_geometry[i];
            auto& r_gp = this->mGlobalPointerNodalMap[r_node.Id()];
            rGPSensitivityVector(i) = r_gp;
            gp_index_map[r_node.Id()] = i;
        }

        // add relevant neighbour gps
        bool found_slip = false;
        IndexType local_index = number_of_nodes;
        for (IndexType i = 0; i < number_of_nodes; ++i) {
            const auto& r_node = r_geometry[i];
            if (r_node.Is(SLIP)) {
                found_slip = true;
                const auto& neighbour_gps = r_node.GetValue(NEIGHBOUR_CONDITION_NODES);
                const auto& neighbour_ids = mNodalNeighboursMap[r_node.Id()];
                for (IndexType i = 0; i < neighbour_ids.size(); ++i) {
                    const auto neighbour_id = neighbour_ids[i];
                    const auto p_itr = gp_index_map.find(neighbour_id);
                    if (p_itr == gp_index_map.end()) {
                        rGPSensitivityVector.push_back(neighbour_gps(i));
                        gp_index_map[neighbour_id] = local_index++;
                    }
                }
            }
        }

        const IndexType derivatives_size = local_index * TDim;
        if (rOutput.size1() != derivatives_size || rOutput.size2() != residuals_size) {
            rOutput.resize(derivatives_size, residuals_size, false);
        }
        rOutput.clear();

        if (found_slip) {
            // calculate entity residual.
            // following methods will throw an error in old adjoint elements/conditions since
            // they does not support SLIP condition based primal solutions
            rEntity.CalculateLocalSystem(aux_matrix, aux_vector, rProcessInfo);
            rEntity.CalculateLocalVelocityContribution(aux_matrix, aux_vector, rProcessInfo);
        }

        // add residual derivative contributions
        for (IndexType a = 0; a < number_of_nodes; ++a) {
            const auto& r_node = r_geometry[a];
            if (r_node.Is(SLIP)) {
                AddNodalRotatedResidualDerivativeToMatrix(
                    rOutput, sensitivity_matrix, aux_vector, a * local_size,
                    gp_index_map, r_node);
            } else {
                AddNodalResidualDerivativeToMatrix(rOutput, sensitivity_matrix, a * local_size);
            }
        }

        KRATOS_CATCH("");
    }

    void AddNodalRotatedResidualDerivativeToMatrix(
        Matrix& rOutput,
        const Matrix& rResidualDerivatives,
        const Vector& rResiduals,
        const IndexType NodeStartIndex,
        const std::unordered_map<IndexType, IndexType>& rDerivativesMap,
        const NodeType& rNode)
    {
        KRATOS_TRY

        using coordinate_transformation_utils = CoordinateTransformationUtils<Matrix, Vector, double>;

        // // get the residual relevant for rNode
        BoundedVector<double, TDim> residual, residual_derivative, aux_vector;
        FluidCalculationUtilities::ReadSubVector<TDim>(residual, rResiduals, NodeStartIndex);

        // get the rotation matrix relevant for rNode
        BoundedMatrix<double, TDim, TDim> rotation_matrix;
        coordinate_transformation_utils::LocalRotationOperatorPure(rotation_matrix, rNode);

        // add rotated residual derivative contributions
        for (IndexType c = 0; c < rResidualDerivatives.size1(); ++c) {
            // get the residual derivative relevant for node
            FluidCalculationUtilities::ReadSubVector<TDim>(
                residual_derivative, row(rResidualDerivatives, c), NodeStartIndex);

            // rotate residual derivative
            noalias(aux_vector) = prod(rotation_matrix, residual_derivative);

            // add rotated residual derivative to local matrix
            FluidCalculationUtilities::AddSubVector<TDim>(
                rOutput, aux_vector, c, NodeStartIndex);

            // add continuity equation derivatives
            rOutput(c, NodeStartIndex + TDim) +=
                rResidualDerivatives(c, NodeStartIndex + TDim);
        }

        // first add rotation matrix derivative contributions w.r.t. rNode
        BoundedMatrix<double, TDim, TDim> rotation_matrix_derivative;
        const int current_node_index = rDerivativesMap.find(rNode.Id())->second * TDim;
        for (IndexType k = 0; k < TDim; ++k) {
            coordinate_transformation_utils::CalculateRotationOperatorPureShapeSensitivities(
                rotation_matrix_derivative, 0, k, rNode);

            noalias(aux_vector) = prod(rotation_matrix_derivative, residual);
            FluidCalculationUtilities::AddSubVector<TDim>(
                rOutput, aux_vector, current_node_index + k, NodeStartIndex);
        }

        // add rotation matrix derivative contributions w.r.t. rNode neighbors
        const auto& r_neighbour_ids = mNodalNeighboursMap[rNode.Id()];
        for (IndexType b = 0; b < r_neighbour_ids.size(); ++b) {
            const int derivative_node_index =
                rDerivativesMap.find(r_neighbour_ids[b])->second * TDim;
            for (IndexType k = 0; k < TDim; ++k) {
                coordinate_transformation_utils::CalculateRotationOperatorPureShapeSensitivities(
                    rotation_matrix_derivative, b + 1, k, rNode);

                noalias(aux_vector) = prod(rotation_matrix_derivative, residual);
                FluidCalculationUtilities::AddSubVector<TDim>(
                    rOutput, aux_vector, derivative_node_index + k, NodeStartIndex);
            }
        }

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

    void AddNodalResidualDerivativeToMatrix(
        Matrix& rOutput,
        const Matrix& rResidualDerivatives,
        const IndexType NodeStartIndex)
    {
        KRATOS_TRY

        // add non-rotated residual derivative contributions
        for (IndexType c = 0; c < rResidualDerivatives.size1(); ++c) {
            for (IndexType i = 0; i < TBlockSize; ++i) {
                rOutput(c, NodeStartIndex + i) +=
                    rResidualDerivatives(c, NodeStartIndex + i);
            }
        }

        KRATOS_CATCH("");
    }

    ///@}

}; /* Class SimpleSteadySensitivityBuilderScheme */

///@}

} /* namespace Kratos.*/

#endif /* KRATOS_SIMPLE_STEADY_SENSITIVITY_BUILDER_SCHEME_H_INCLUDED defined */
