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

// External includes

// Project includes
#include "containers/global_pointers_vector.h"
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/schemes/sensitivity_builder_scheme.h"
#include "utilities/openmp_utils.h"

// Application includes
#include "custom_utilities/fluid_adjoint_slip_utilities.h"

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
          mAdjointSlipUtilities(Dimension, BlockSize),
          mDimension(Dimension),
          mBlockSize(BlockSize)
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

            mAdjointSlipUtilities.Initialize(rModelPart);
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
        const auto k = OpenMPUtils::ThisThread();
        CalculateResidualSensitivityMatrix<ElementType, array_1d<double, 3>>(
            rElement, rAdjointValues, mAuxVectors[k], mSensitivityMatrices[k], rOutput, mAuxMatrices[k], rGPSensitivityVector, rVariable,
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
            rCondition, rAdjointValues, mAuxVectors[k], mSensitivityMatrices[k], rOutput, mAuxMatrices[k], rGPSensitivityVector, rVariable,
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
    ///@name Protected Member Variables
    ///@{

    FluidAdjointSlipUtilities mAdjointSlipUtilities;

    ///@}
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
        auto& residuals = mAuxVectors[k];
        auto& aux_matrix = mAuxMatrices[k];

        this->CalculateResidualSensitivityMatrix<TEntityType, TDataType>(
                rEntity, adjoint_vector, residuals,sensitivity_matrix, rotated_sensitivity_matrix, aux_matrix,
                rGPSensitivityVector, rVariable, rProcessInfo);

        if (adjoint_vector.size() != 0) {

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
        Vector& rEntityResiduals,
        Matrix& rEntityResidualDerivatives,
        Matrix& rEntityRotatedResidualDerivatives,
        Matrix& rAuxMatrix,
        GlobalPointersVector<NodeType>& rGPSensitivityVector,
        const Variable<TDataType>& rVariable,
        const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        // get adjoint solution vector
        rEntity.GetValuesVector(rAdjointValues);
        const IndexType residuals_size = rAdjointValues.size();

        if (residuals_size != 0) {
            // calculate entity residual derivatives
            rEntity.CalculateSensitivityMatrix(rVariable, rEntityResidualDerivatives, rProcessInfo);

            if (rVariable == SHAPE_SENSITIVITY) {
                const auto& r_geometry = rEntity.GetGeometry();
                const IndexType number_of_nodes = r_geometry.PointsNumber();

                bool found_slip = false;
                for (IndexType i = 0; i < number_of_nodes; ++i) {
                    const auto& r_node = r_geometry[i];
                    if (r_node.Is(SLIP)) {
                        found_slip = true;
                        break;
                    }
                }

                if (found_slip) {
                    // calculate entity residual.
                    // following methods will throw an error in old adjoint elements/conditions since
                    // they does not support SLIP condition based primal solutions
                    this->CalculateLHSAndRHS(rEntity, rAuxMatrix, rEntityResiduals, rProcessInfo);
                } else {
                    rEntityResiduals = ZeroVector(residuals_size);
                }

                mAdjointSlipUtilities.CalculateRotatedSlipConditionAppliedShapeVariableDerivatives(rEntityRotatedResidualDerivatives, rGPSensitivityVector, rEntityResiduals, rEntityResidualDerivatives, rEntity, rProcessInfo);
            } else {
                if (rEntityRotatedResidualDerivatives.size1() != rEntityResidualDerivatives.size1() || rEntityRotatedResidualDerivatives.size2() != rEntityResidualDerivatives.size2()) {
                    rEntityRotatedResidualDerivatives.resize(rEntityResidualDerivatives.size1(), rEntityResidualDerivatives.size2(), false);
                }
                noalias(rEntityRotatedResidualDerivatives) = rEntityResidualDerivatives;
            }
        }

        KRATOS_CATCH("");
    }
    ///@}

}; /* Class SimpleSteadySensitivityBuilderScheme */

///@}

} /* namespace Kratos.*/

#endif /* KRATOS_SIMPLE_STEADY_SENSITIVITY_BUILDER_SCHEME_H_INCLUDED defined */
