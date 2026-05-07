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
#include "utilities/sensitivity_utilities.h"

// Application includes
#include "custom_utilities/fluid_calculation_utilities.h"

namespace Kratos
{
///@name Kratos Classes
///@{

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

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    SimpleSteadySensitivityBuilderScheme(
        const IndexType Dimension,
        const IndexType BlockSize)
        : SensitivityBuilderScheme(),
          mDimension(Dimension),
          mBlockSize(BlockSize),
          mRotationalTool(Dimension, mBlockSize)
    {
        KRATOS_TRY

        // Allocate auxiliary memory.
        // This needs to be done in the constructor because, this scheme
        // is used to calculate sensitivities w.r.t. element quantities
        // and in there we don't usually pass the model part. Hence, we
        // can not call SimpleSteadySensitivityBuilderScheme::Initialize
        // method.
        const int number_of_threads = ParallelUtilities::GetNumThreads();
        mAuxVectors.resize(number_of_threads);
        mAuxMatrices.resize(number_of_threads);
        mRotatedSensitivityMatrices.resize(number_of_threads);
        mSensitivityMatrices.resize(number_of_threads);

        if (Dimension == 2) {
            this->mAddNodalRotationDerivativesMethod = &SimpleSteadySensitivityBuilderScheme::TemplatedAddNodalRotationDerivatives<2>;
            this->mAddNodalApplySlipConditionDerivativesMethod = &SimpleSteadySensitivityBuilderScheme::TemplatedAddNodalApplySlipConditionDerivatives<2>;
        } else if (Dimension == 3) {
            this->mAddNodalRotationDerivativesMethod = &SimpleSteadySensitivityBuilderScheme::TemplatedAddNodalRotationDerivatives<3>;
            this->mAddNodalApplySlipConditionDerivativesMethod = &SimpleSteadySensitivityBuilderScheme::TemplatedAddNodalApplySlipConditionDerivatives<3>;
        } else {
            KRATOS_ERROR << "Unsupported dimensionality requested. Only 2D and 3D "
                            "supported. [ Dimension = "
                        << Dimension << " ].\n";
        }

        KRATOS_INFO(this->Info()) << this->Info() << " created [ Dimensionality = " << mDimension << ", BlockSize = " << mBlockSize << " ].\n";

        KRATOS_CATCH("");
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
                rModelPart,
                NEIGHBOUR_CONDITION_NODES);
            neighbor_process.Execute();
            mNodalNeighboursMap = neighbor_process.GetNeighbourIds(rModelPart.Nodes());

            // Assign condition normal derivatives to corresponding nodes
            // NORMAL_SHAPE_DERIVATIVE on conditions should already be available before
            // executing the following command. (This is done in adjoint_vmsmonolithic_solver.py)
            SensitivityUtilities::AssignEntityDerivativesToNodes<ModelPart::ConditionsContainerType>(
                rModelPart, mDimension, NORMAL_SHAPE_DERIVATIVE,
                mNodalNeighboursMap, 1.0 / mDimension, SLIP);
        }

        BaseType::InitializeSolutionStep(rModelPart, rSensitivityModelPart, rResponseFunction);

        KRATOS_CATCH("");
    }

    /**
     * @brief Calculates sensitivity from a given element
     *
     * This method needs to be implemented in the derived class if sensitivity calculation is
     * based on nodal quantities, and elements has a contribution towards that.
     *
     * rSensitivity vector should have sensitivities in the order matching rGPSensitivity vector.
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
     * This method needs to be implemented in the derived class if sensitivity calculation is
     * based on element quantities, and elements has a contribution towards that.
     *
     * rSensitivity vector should have sensitivities in the order matching rGPSensitivity vector.
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
     * This method needs to be implemented in the derived class if sensitivity calculation is
     * based on nodal quantities, and conditions has a contribution towards that.
     *
     * rSensitivity vector should have sensitivities in the order matching rGPSensitivity vector.
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
     * This method needs to be implemented in the derived class if sensitivity calculation is
     * based on condition quantities, and conditions has a contribution towards that.
     *
     * rSensitivity vector should have sensitivities in the order matching rGPSensitivity vector.
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
     * This method needs to be implemented in the derived class if sensitivity calculation is
     * based on nodal quantities, and elements has a contribution towards that.
     *
     * rSensitivity vector should have sensitivities in the order matching rGPSensitivity vector.
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
     * This method needs to be implemented in the derived class if sensitivity calculation is
     * based on element quantities, and elements has a contribution towards that.
     *
     * rSensitivity vector should have sensitivities in the order matching rGPSensitivity vector.
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
     * This method needs to be implemented in the derived class if sensitivity calculation is
     * based on nodal quantities, and conditions has a contribution towards that.
     *
     * rSensitivity vector should have sensitivities in the order matching rGPSensitivity vector.
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
     * This method needs to be implemented in the derived class if sensitivity calculation is
     * based on condition quantities, and conditions has a contribution towards that.
     *
     * rSensitivity vector should have sensitivities in the order matching rGPSensitivity vector.
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
        const auto k = OpenMPUtils::ThisThread();
        CalculateResidualSensitivityMatrix<ElementType, array_1d<double, 3>>(
            rElement, rAdjointValues, mSensitivityMatrices[k], rOutput, rGPSensitivityVector, rVariable,
            rCurrentProcessInfo);
    }

    void CalculateResidualSensitivityMatrix(
        ConditionType& rCondition,
        Vector& rAdjointValues,
        Matrix& rOutput,
        GlobalPointersVector<NodeType>& rGPSensitivityVector,
        const Variable<array_1d<double, 3>>& rVariable,
        const ProcessInfo& rCurrentProcessInfo)
    {
        const auto k = OpenMPUtils::ThisThread();
        CalculateResidualSensitivityMatrix<ConditionType, array_1d<double, 3>>(
            rCondition, rAdjointValues, mSensitivityMatrices[k], rOutput, rGPSensitivityVector, rVariable,
            rCurrentProcessInfo);
    }

    void Clear() override
    {
        BaseType::Clear();
        mAuxVectors.clear();
        mAuxMatrices.clear();
        mRotatedSensitivityMatrices.clear();
        mSensitivityMatrices.clear();
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

protected:
    ///@name Protected Operations
    ///@{

    virtual void CalculateLHSAndRHS(
        ElementType& rElement,
        Matrix& rLHS,
        Vector& rRHS,
        const ProcessInfo& rProcessInfo)
    {
        // following calls uses the same method calls as in the primal scheme to be consistent
        rElement.CalculateLocalSystem(rLHS, rRHS, rProcessInfo);
        rElement.CalculateLocalVelocityContribution(rLHS, rRHS, rProcessInfo);
    }

    virtual void CalculateLHSAndRHS(
        ConditionType& rCondition,
        Matrix& rLHS,
        Vector& rRHS,
        const ProcessInfo& rProcessInfo)
    {
        // following calls uses the same method calls as in the primal scheme to be consistent
        rCondition.CalculateLocalSystem(rLHS, rRHS, rProcessInfo);
        rCondition.CalculateLocalVelocityContribution(rLHS, rRHS, rProcessInfo);
    }

    ///@}

private:
    ///@name Member Variables
    ///@{

    const IndexType mDimension;
    const IndexType mBlockSize;

    bool mIsNodalNormalShapeDerivativesComputed = false;
    std::vector<Matrix> mAuxMatrices;
    std::vector<Vector> mAuxVectors;
    std::vector<Matrix> mRotatedSensitivityMatrices;
    std::vector<Matrix> mSensitivityMatrices;

    std::unordered_map<int, std::vector<int>> mNodalNeighboursMap;

    const CoordinateTransformationUtils<Matrix, Vector, double> mRotationalTool;

    void (SimpleSteadySensitivityBuilderScheme::*mAddNodalRotationDerivativesMethod)(
        Matrix&,
        const Matrix&,
        const Vector&,
        const IndexType,
        const std::unordered_map<IndexType, IndexType>&,
        const NodeType&) const;

    void (SimpleSteadySensitivityBuilderScheme::*mAddNodalApplySlipConditionDerivativesMethod)(
        Matrix&,
        const IndexType,
        const std::unordered_map<IndexType, IndexType>&,
        const NodeType&) const;

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
        auto& sensitivity_matrix = mSensitivityMatrices[k];

        // get adjoint solution vector
        rEntity.GetValuesVector(adjoint_vector);
        const IndexType residuals_size = adjoint_vector.size();

        if (residuals_size != 0) {
            this->CalculateResidualSensitivityMatrix<TEntityType, TDataType>(
                rEntity, adjoint_vector, sensitivity_matrix, rotated_sensitivity_matrix,
                rGPSensitivityVector, rVariable, rProcessInfo);

            if (rSensitivityVector.size() != rotated_sensitivity_matrix.size1()) {
                rSensitivityVector.resize(rotated_sensitivity_matrix.size1(), false);
            }
            noalias(rSensitivityVector) = prod(rotated_sensitivity_matrix, adjoint_vector);

            // add objective derivative contributions
            auto& objective_partial_sensitivity = mPartialSensitivity[k];
            rResponseFunction.CalculatePartialSensitivity(
                rEntity, rVariable, sensitivity_matrix,
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
        Matrix& rEntityResidualDerivatives,
        Matrix& rEntityRotatedResidualDerivatives,
        GlobalPointersVector<NodeType>& rGPSensitivityVector,
        const Variable<TDataType>& rVariable,
        const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        using DofsVectorType = std::vector<Dof<double>::Pointer>;

        const auto k = OpenMPUtils::ThisThread();

        auto& aux_vector = mAuxVectors[k];
        auto& aux_matrix = mAuxMatrices[k];

        // get adjoint solution vector
        rEntity.GetValuesVector(rAdjointValues);
        const IndexType residuals_size = rAdjointValues.size();

        // calculate entity residual derivatives
        rEntity.CalculateSensitivityMatrix(rVariable, rEntityResidualDerivatives, rProcessInfo);

        KRATOS_ERROR_IF(residuals_size != rEntityResidualDerivatives.size2())
            << "mAdjointVectors.size(): " << residuals_size
            << " incompatible with rEntityResidualDerivatives.size1(): "
            << rEntityResidualDerivatives.size2() << ". Variable: " << rVariable << std::endl;

        // a map to store node id as key and the node order as the value
        // dofs are used in here, because if wall distance is calculated in the condition
        // then there will be derivative contributions coming to domain internal node as well
        // hence, list of dofs are used in here.
        DofsVectorType dofs;
        rEntity.GetDofList(dofs, rProcessInfo);

        // get derivative node ids
        std::vector<int> derivative_node_ids;
        for (IndexType i = 0; i < dofs.size(); ++i) {
            if (std::find(derivative_node_ids.begin(), derivative_node_ids.end(), dofs[i]->Id()) == derivative_node_ids.end()) {
                derivative_node_ids.push_back(dofs[i]->Id());
            }
        }

        std::unordered_map<IndexType, IndexType> gp_index_map;
        auto& r_geometry = rEntity.GetGeometry();
        const IndexType number_of_nodes = r_geometry.PointsNumber();
        const IndexType number_of_derivative_nodes = derivative_node_ids.size();

        if (rGPSensitivityVector.size() != number_of_derivative_nodes) {
            rGPSensitivityVector.resize(number_of_derivative_nodes);
        }

        // add geometry gps
        for (IndexType i = 0; i < number_of_derivative_nodes; ++i) {
            auto& r_gp = this->mGlobalPointerNodalMap[derivative_node_ids[i]];
            rGPSensitivityVector(i) = r_gp;
            gp_index_map[derivative_node_ids[i]] = i;
        }

        if (rVariable == SHAPE_SENSITIVITY) {
            KRATOS_DEBUG_ERROR_IF(number_of_derivative_nodes * mDimension != rEntityResidualDerivatives.size1())
                << "Entity sensitivity matrix size mismatch. [ rEntityResidualDerivatives.size = ( " << rEntityResidualDerivatives.size1()
                << ", " << rEntityResidualDerivatives.size2() << " ), required size = ( " << (number_of_derivative_nodes * mDimension)
                << ", " << residuals_size << " ) ].\n";

            // add relevant neighbour gps
            bool found_slip = false;
            IndexType local_index = number_of_derivative_nodes;
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

            const IndexType derivatives_size = local_index * mDimension;
            if (rEntityRotatedResidualDerivatives.size1() != derivatives_size || rEntityRotatedResidualDerivatives.size2() != residuals_size) {
                rEntityRotatedResidualDerivatives.resize(derivatives_size, residuals_size, false);
            }
            rEntityRotatedResidualDerivatives.clear();

            if (found_slip) {
                // calculate entity residual.
                // following methods will throw an error in old adjoint elements/conditions since
                // they does not support SLIP condition based primal solutions
                this->CalculateLHSAndRHS(rEntity, aux_matrix, aux_vector, rProcessInfo);
            }

            // add residual derivative contributions
            for (IndexType a = 0; a < number_of_nodes; ++a) {
                const auto& r_node = r_geometry[a];
                const IndexType block_index = a * mBlockSize;
                if (r_node.Is(SLIP)) {
                    AddNodalRotationDerivatives(rEntityRotatedResidualDerivatives, rEntityResidualDerivatives, aux_vector, block_index, gp_index_map, r_node);
                    AddNodalApplySlipConditionDerivatives(rEntityRotatedResidualDerivatives, block_index, gp_index_map, r_node);
                } else {
                    AddNodalResidualDerivatives(rEntityRotatedResidualDerivatives, rEntityResidualDerivatives, block_index);
                }
            }
        } else {
            const IndexType derivatives_size = number_of_nodes * mDimension;
            if (rEntityRotatedResidualDerivatives.size1() != derivatives_size || rEntityRotatedResidualDerivatives.size2() != residuals_size) {
                rEntityRotatedResidualDerivatives.resize(derivatives_size, residuals_size, false);
            }
            noalias(rEntityRotatedResidualDerivatives) = rEntityResidualDerivatives;
        }

        KRATOS_CATCH("");
    }

    void AddNodalRotationDerivatives(
        Matrix& rOutput,
        const Matrix& rResidualDerivatives,
        const Vector& rResiduals,
        const IndexType NodeStartIndex,
        const std::unordered_map<IndexType, IndexType>& rDerivativesMap,
        const NodeType& rNode) const
    {
        (this->*(this->mAddNodalRotationDerivativesMethod))(rOutput, rResidualDerivatives, rResiduals, NodeStartIndex, rDerivativesMap, rNode);
    }

    template<unsigned int TDim>
    void TemplatedAddNodalRotationDerivatives(
        Matrix& rOutput,
        const Matrix& rResidualDerivatives,
        const Vector& rResiduals,
        const IndexType NodeStartIndex,
        const std::unordered_map<IndexType, IndexType>& rDerivativesMap,
        const NodeType& rNode) const
    {
        KRATOS_TRY

        // // get the residual relevant for rNode
        BoundedVector<double, TDim> residual, residual_derivative, aux_vector;
        FluidCalculationUtilities::ReadSubVector<TDim>(residual, rResiduals, NodeStartIndex);

        // get the rotation matrix relevant for rNode
        BoundedMatrix<double, TDim, TDim> rotation_matrix;
        mRotationalTool.LocalRotationOperatorPure(rotation_matrix, rNode);

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

            // add rest of the equation derivatives
            for (IndexType a = TDim; a < mBlockSize; ++a) {
                rOutput(c, NodeStartIndex + a) +=
                    rResidualDerivatives(c, NodeStartIndex + a);
            }
        }

        // first add rotation matrix derivative contributions w.r.t. rNode
        BoundedMatrix<double, TDim, TDim> rotation_matrix_derivative;
        const int current_node_index = rDerivativesMap.find(rNode.Id())->second * TDim;
        for (IndexType k = 0; k < TDim; ++k) {
            mRotationalTool.CalculateRotationOperatorPureShapeSensitivities(
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
                mRotationalTool.CalculateRotationOperatorPureShapeSensitivities(
                    rotation_matrix_derivative, b + 1, k, rNode);

                noalias(aux_vector) = prod(rotation_matrix_derivative, residual);
                FluidCalculationUtilities::AddSubVector<TDim>(
                    rOutput, aux_vector, derivative_node_index + k, NodeStartIndex);
            }
        }

        KRATOS_CATCH("");
    }

    void AddNodalApplySlipConditionDerivatives(
        Matrix& rOutput,
        const IndexType NodeStartIndex,
        const std::unordered_map<IndexType, IndexType>& rDerivativesMap,
        const NodeType& rNode) const
    {
        (this->*(this->mAddNodalApplySlipConditionDerivativesMethod))(rOutput, NodeStartIndex, rDerivativesMap, rNode);
    }

    template<unsigned int TDim>
    void TemplatedAddNodalApplySlipConditionDerivatives(
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

    void AddNodalResidualDerivatives(
        Matrix& rOutput,
        const Matrix& rResidualDerivatives,
        const IndexType NodeStartIndex) const
    {
        KRATOS_TRY

        // add non-rotated residual derivative contributions
        for (IndexType c = 0; c < rResidualDerivatives.size1(); ++c) {
            for (IndexType i = 0; i < mBlockSize; ++i) {
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
