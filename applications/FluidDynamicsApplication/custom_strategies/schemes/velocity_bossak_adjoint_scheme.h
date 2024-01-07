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

#if !defined(KRATOS_VELOCITY_BOSSAK_ADJOINT_SCHEME_H_INCLUDED)
#define KRATOS_VELOCITY_BOSSAK_ADJOINT_SCHEME_H_INCLUDED

// System includes
#include <vector>
#include <string>
#include <unordered_map>

// External includes

// Project includes
#include "includes/checks.h"
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "response_functions/adjoint_response_function.h"
#include "solving_strategies/schemes/residual_based_adjoint_bossak_scheme.h"
#include "solving_strategies/schemes/scheme.h"
#include "utilities/indirect_scalar.h"
#include "utilities/parallel_utilities.h"
#include "utilities/variable_utils.h"

// Application includes
#include "custom_processes/element_refinement_process.h"
#include "custom_utilities/fluid_adjoint_slip_utilities.h"

namespace Kratos
{
///@name Kratos Classes
///@{

template <class TSparseSpace, class TDenseSpace>
class VelocityBossakAdjointScheme : public ResidualBasedAdjointBossakScheme<TSparseSpace, TDenseSpace>
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(VelocityBossakAdjointScheme);

    using BaseType = ResidualBasedAdjointBossakScheme<TSparseSpace, TDenseSpace>;

    using LocalSystemVectorType = typename BaseType::LocalSystemVectorType;

    using LocalSystemMatrixType = typename BaseType::LocalSystemMatrixType;

    using DofsArrayType = typename BaseType::DofsArrayType;

    using BossakConstants = typename BaseType::BossakConstants;

    using TSystemVectorType = typename BaseType::TSystemVectorType;

    using TSystemMatrixType = typename BaseType::TSystemMatrixType;

    using EquationIdVectorType = std::vector<IndexType>;

    using AssembledVectorType = std::unordered_map<IndexType, double>;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    VelocityBossakAdjointScheme(
        Parameters Settings,
        AdjointResponseFunction::Pointer pResponseFunction,
        const IndexType Dimension,
        const IndexType BlockSize,
        ElementRefinementProcess::Pointer pElementRefinementProcess = nullptr,
        const IndexType EchoLevel = 0)
        : BaseType(Settings, pResponseFunction),
          mAdjointSlipUtilities(Dimension, BlockSize),
          mpElementRefinementProcess(pElementRefinementProcess),
          mEchoLevel(EchoLevel)
    {
        KRATOS_INFO(this->Info()) << this->Info() << " created [ Dimensionality = " << Dimension << ", BlockSize = " << BlockSize << " ].\n";
    }

    /// Destructor.
    ~VelocityBossakAdjointScheme() override = default;

    ///@}
    ///@name Operations
    ///@{

    void Initialize(ModelPart& rModelPart) override
    {
        KRATOS_TRY;

        BaseType::Initialize(rModelPart);

        // Allocate auxiliary memory.
        int num_threads = ParallelUtilities::GetNumThreads();
        mAuxiliaryMatrix.resize(num_threads);
        mRotatedMatrix.resize(num_threads);

        if (mpElementRefinementProcess) {
            mpElementRefinementProcess->ExecuteInitialize();

            // now set NUMBER_OF_NEIGHBOUR_ELEMENTS for the refined model parts
            IndexPartition<int>(num_threads).for_each([&](const int) {
                auto& r_thread_local_refined_model_part = mpElementRefinementProcess->GetThreadLocalModelPart();

                // initialize everything to zero
                for (auto& r_node : r_thread_local_refined_model_part.Nodes()) {
                    r_node.SetValue(NUMBER_OF_NEIGHBOUR_ELEMENTS, 0);
                }

                for (auto& r_element : r_thread_local_refined_model_part.Elements()) {
                    for (auto& r_node : r_element.GetGeometry()) {
                        r_node.GetValue(NUMBER_OF_NEIGHBOUR_ELEMENTS) += 1;
                    }
                }
            });

            // now create data holders for each coarse element nodal previous time step values
            for (const auto& r_element : rModelPart.Elements()) {
                mRefinedNodalData[r_element.Id()];
            }

            const IndexType number_of_dofs_per_node = GetNumberOfDofsPerNode(rModelPart);
            const Vector zero_vector_values(number_of_dofs_per_node, 0.0);
            block_for_each(rModelPart.Elements(), [&](Element& rElement){
                auto& r_thread_local_refined_model_part = mpElementRefinementProcess->GetThreadLocalModelPart();
                auto& r_nodal_data_map = mRefinedNodalData.find(rElement.Id())->second;

                for (auto& r_node : r_thread_local_refined_model_part.Nodes()) {
                    auto& r_nodal_data = r_nodal_data_map[r_node.Id()];
                    r_nodal_data.Mu = zero_vector_values;
                    r_nodal_data.D = zero_vector_values;
                    r_nodal_data.E = zero_vector_values;
                }

                rElement.SetValue(RESPONSE_FUNCTION_INTERPOLATION_ERROR_AUXILIARY_1, zero_vector_values);
                rElement.SetValue(RESPONSE_FUNCTION_INTERPOLATION_ERROR_AUXILIARY_2, zero_vector_values);
            });

        }

        KRATOS_CATCH("");
    }

    void InitializeSolutionStep(
        ModelPart& rModelPart,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b) override
    {
        KRATOS_TRY

        BaseType::InitializeSolutionStep(rModelPart, A, Dx, b);

        // calculate constants
        mC1 = 1.0 / (mBossak.Gamma * mBossak.DeltaTime);
        mC2 = 1.0 / (mBossak.Gamma * mBossak.Gamma * mBossak.DeltaTime);
        mC3 = (mBossak.Gamma - 1.0) / mBossak.Gamma;
        mC4 = 1.0 - mBossak.Alpha;
        mC5 = (1.0 - mBossak.Alpha) / mBossak.Alpha;

        if (mpElementRefinementProcess) {
            int num_threads = ParallelUtilities::GetNumThreads();

            IndexPartition<int>(num_threads).for_each([&](const int) {
                auto& r_thread_local_refined_model_part = mpElementRefinementProcess->GetThreadLocalModelPart();
                r_thread_local_refined_model_part.GetProcessInfo() = rModelPart.GetProcessInfo();
            });
        }

        KRATOS_CATCH("")
    }

    void FinalizeSolutionStep(
        ModelPart& rModelPart,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b) override
    {
        KRATOS_TRY

        BaseType::FinalizeSolutionStep(rModelPart, A, Dx, b);

        if (mpElementRefinementProcess) {
            CalculateResponseFunctionInterpolationError(rModelPart);
        }

        KRATOS_CATCH("")
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "VelocityBossakAdjointScheme";
    }

    ///@}

protected:
    ///@name Protected Operations
    ///@{

    void CalculateGradientContributions(
        Element& rElement,
        LocalSystemMatrixType& rLHS_Contribution,
        LocalSystemVectorType& rRHS_Contribution,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        CalculateEntityGradientContributions(
            rElement, rLHS_Contribution, rRHS_Contribution, rCurrentProcessInfo);
    }

    void CalculateGradientContributions(
        Condition& rCondition,
        LocalSystemMatrixType& rLHS_Contribution,
        LocalSystemVectorType& rRHS_Contribution,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        CalculateEntityGradientContributions(
            rCondition, rLHS_Contribution, rRHS_Contribution, rCurrentProcessInfo);
    }

    void CalculateFirstDerivativeContributions(
        Element& rElement,
        LocalSystemMatrixType& rLHS_Contribution,
        LocalSystemVectorType& rRHS_Contribution,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        CalculateEntityFirstDerivativeContributions(
            rElement, rLHS_Contribution, rRHS_Contribution, rCurrentProcessInfo);
    }

    void CalculateFirstDerivativeContributions(
        Condition& rCondition,
        LocalSystemMatrixType& rLHS_Contribution,
        LocalSystemVectorType& rRHS_Contribution,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        CalculateEntityFirstDerivativeContributions(
            rCondition, rLHS_Contribution, rRHS_Contribution, rCurrentProcessInfo);
    }

    void CalculateSecondDerivativeContributions(
        Element& rElement,
        LocalSystemMatrixType& rLHS_Contribution,
        LocalSystemVectorType& rRHS_Contribution,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        CalculateEntitySecondDerivativeContributions(
            rElement, rLHS_Contribution, rRHS_Contribution, rCurrentProcessInfo);
    }

    void CalculateSecondDerivativeContributions(
        Condition& rCondition,
        LocalSystemMatrixType& rLHS_Contribution,
        LocalSystemVectorType& rRHS_Contribution,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        CalculateEntitySecondDerivativeContributions(
            rCondition, rLHS_Contribution, rRHS_Contribution, rCurrentProcessInfo);
    }

    void CalculateTimeSchemeContributions(
        Element& rElement,
        LocalSystemVectorType& rAdjointTimeSchemeValues2,
        LocalSystemVectorType& rAdjointTimeSchemeValues3,
        AdjointResponseFunction& rAdjointResponseFunction,
        const BossakConstants& rBossakConstants,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        CalculateEntityTimeSchemeContributions(rElement, rAdjointTimeSchemeValues2,
                                               rAdjointTimeSchemeValues3,
                                               rCurrentProcessInfo);
    }

    void CalculateTimeSchemeContributions(
        Condition& rCondition,
        LocalSystemVectorType& rAdjointTimeSchemeValues2,
        LocalSystemVectorType& rAdjointTimeSchemeValues3,
        AdjointResponseFunction& rAdjointResponseFunction,
        const BossakConstants& rBossakConstants,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        CalculateEntityTimeSchemeContributions(rCondition, rAdjointTimeSchemeValues2,
                                               rAdjointTimeSchemeValues3,
                                               rCurrentProcessInfo);
    }

    void CalculateAuxiliaryVariableContributions(
        Element& rElement,
        LocalSystemVectorType& rAdjointAuxiliaryValues,
        AdjointResponseFunction& rAdjointResponseFunction,
        const BossakConstants& rBossakConstants,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        CalculateEntityAuxiliaryVariableContributions(
            rElement, rAdjointAuxiliaryValues, rCurrentProcessInfo);
    }

    void CalculateAuxiliaryVariableContributions(
        Condition& rCondition,
        LocalSystemVectorType& rAdjointAuxiliaryValues,
        AdjointResponseFunction& rAdjointResponseFunction,
        const BossakConstants& rBossakConstants,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        CalculateEntityAuxiliaryVariableContributions(
            rCondition, rAdjointAuxiliaryValues, rCurrentProcessInfo);
    }

    void CheckAndResizeThreadStorage(unsigned SystemSize) override
    {
        const int k = OpenMPUtils::ThisThread();

        BaseType::CheckAndResizeThreadStorage(SystemSize);

        if (mAuxiliaryMatrix[k].size1() != SystemSize ||
            mAuxiliaryMatrix[k].size2() != SystemSize) {
            mAuxiliaryMatrix[k].resize(SystemSize, SystemSize, false);
        }

        if (mRotatedMatrix[k].size1() != SystemSize || mRotatedMatrix[k].size2() != SystemSize) {
            mRotatedMatrix[k].resize(SystemSize, SystemSize, false);
        }
    }

    ///@}

private:
    ///@name Private Members
    ///@{

    using AssemblyDataType = std::unordered_map<IndexType, std::tuple<IndexType, double, double>>;

    using BaseType::mBossak;

    std::vector<Matrix> mAuxiliaryMatrix;
    std::vector<Matrix> mRotatedMatrix;

    const FluidAdjointSlipUtilities mAdjointSlipUtilities;

    ElementRefinementProcess::Pointer mpElementRefinementProcess;

    const IndexType mEchoLevel;

    struct RefinedNodeAdjointDataStorage
    {
        Vector Mu;
        Vector D;
        Vector E;
    };

    // corse_element_id, coarse element's refined node id, storage
    std::unordered_map<IndexType, std::unordered_map<IndexType, RefinedNodeAdjointDataStorage>> mRefinedNodalData;

    // constants
    // mC1 = 1.0 / (mBossak.Gamma * mBossak.DeltaTime)
    double mC1;

    // mC2 = 1.0 / (mBossak.Gamma * mBossak.Gamma * mBossak.DeltaTime)
    double mC2;

    // mC3 = (mGamma - 1) / mGamma
    double mC3;

    // mC4 = 1.0 - mBossak.Alpha
    double mC4;

    // mC5 = (1.0 - mBossak.Alpha) / mBossak.Alpha
    double mC5;

    ///@}
    ///@name Private Operations
    ///@{

    IndexType GetNumberOfDofsPerNode(ModelPart& rModelPart)
    {
        KRATOS_TRY

        IndexType number_of_dofs_per_node{0};
        if (rModelPart.NumberOfElements() > 0) {
            EquationIdVectorType equation_ids;
            auto& r_dummy_element = rModelPart.Elements().front();
            r_dummy_element.EquationIdVector(equation_ids, rModelPart.GetProcessInfo());
            number_of_dofs_per_node = equation_ids.size() / r_dummy_element.GetGeometry().size();
        }

        return rModelPart.GetCommunicator().GetDataCommunicator().MaxAll(number_of_dofs_per_node);

        KRATOS_CATCH("");
    }

    template<class TEntityType>
    void CalculateEntityGradientContributions(
        TEntityType& rCurrentEntity,
        LocalSystemMatrixType& rLHS_Contribution,
        LocalSystemVectorType& rRHS_Contribution,
        const ProcessInfo& rCurrentProcessInfo)
    {
        const int k = OpenMPUtils::ThisThread();

        auto& aux_matrix = this->mAuxiliaryMatrix[k];

        rCurrentEntity.CalculateLeftHandSide(aux_matrix, rCurrentProcessInfo);
        this->mpResponseFunction->CalculateGradient(
            rCurrentEntity, aux_matrix, rRHS_Contribution, rCurrentProcessInfo);

        mAdjointSlipUtilities.CalculateRotatedSlipConditionAppliedNonSlipNonShapeVariableDerivatives(
            rLHS_Contribution, aux_matrix, rCurrentEntity.GetGeometry());

        noalias(rRHS_Contribution) = -rRHS_Contribution;
    }

    template<class TEntityType>
    void CalculateEntityFirstDerivativeContributions(
        TEntityType& rCurrentEntity,
        LocalSystemMatrixType& rLHS_Contribution,
        LocalSystemVectorType& rRHS_Contribution,
        const ProcessInfo& rCurrentProcessInfo)
    {
        const int k = OpenMPUtils::ThisThread();
        auto& aux_matrix = this->mAuxiliaryMatrix[k];
        auto& rotated_matrix = this->mRotatedMatrix[k];
        auto& response_first_derivatives = this->mFirstDerivsResponseGradient[k];

        rCurrentEntity.CalculateFirstDerivativesLHS(aux_matrix, rCurrentProcessInfo);
        this->mpResponseFunction->CalculateFirstDerivativesGradient(
            rCurrentEntity, aux_matrix, response_first_derivatives, rCurrentProcessInfo);

        mAdjointSlipUtilities.CalculateRotatedSlipConditionAppliedSlipVariableDerivatives(
            rotated_matrix, aux_matrix, rCurrentEntity.GetGeometry());

        noalias(rLHS_Contribution) += this->mBossak.C6 * rotated_matrix;
        noalias(rRHS_Contribution) -= this->mBossak.C6 * response_first_derivatives;
    }

    template<class TEntityType>
    void CalculateEntitySecondDerivativeContributions(
        TEntityType& rCurrentEntity,
        LocalSystemMatrixType& rLHS_Contribution,
        LocalSystemVectorType& rRHS_Contribution,
        const ProcessInfo& rCurrentProcessInfo)
    {
        const int k = OpenMPUtils::ThisThread();
        auto& aux_matrix = this->mAuxiliaryMatrix[k];
        auto& rotated_matrix = this->mRotatedMatrix[k];
        auto& response_second_derivatives = this->mSecondDerivsResponseGradient[k];

        rCurrentEntity.CalculateSecondDerivativesLHS(aux_matrix, rCurrentProcessInfo);
        aux_matrix *= (1.0 - this->mBossak.Alpha);
        this->mpResponseFunction->CalculateSecondDerivativesGradient(
            rCurrentEntity, aux_matrix, response_second_derivatives, rCurrentProcessInfo);

        mAdjointSlipUtilities.CalculateRotatedSlipConditionAppliedNonSlipNonShapeVariableDerivatives(
            rotated_matrix, aux_matrix, rCurrentEntity.GetGeometry());

        noalias(rLHS_Contribution) += this->mBossak.C7 * rotated_matrix;
        noalias(rRHS_Contribution) -= this->mBossak.C7 * response_second_derivatives;
    }

    template<class TEntityType>
    void CalculateEntityTimeSchemeContributions(
        TEntityType& rCurrentEntity,
        LocalSystemVectorType& rAdjointTimeSchemeValues2,
        LocalSystemVectorType& rAdjointTimeSchemeValues3,
        const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY

        const int k = OpenMPUtils::ThisThread();
        auto& adjoint_values = this->mAdjointValuesVector[k];
        auto& aux_matrix = this->mAuxiliaryMatrix[k];
        auto& entity_first_derivatives = this->mFirstDerivsLHS[k];
        auto& entity_second_derivatives = this->mSecondDerivsLHS[k];
        auto& response_first_derivatives = this->mFirstDerivsResponseGradient[k];
        auto& response_second_derivatives = this->mSecondDerivsResponseGradient[k];

        const auto& r_const_entity_ref = rCurrentEntity;
        r_const_entity_ref.GetValuesVector(adjoint_values);
        this->CheckAndResizeThreadStorage(adjoint_values.size());

        /// starting to build residual for next time step calculations
        rCurrentEntity.CalculateFirstDerivativesLHS(aux_matrix, rProcessInfo);
        this->mpResponseFunction->CalculateFirstDerivativesGradient(
            rCurrentEntity, aux_matrix, response_first_derivatives, rProcessInfo);

        mAdjointSlipUtilities.CalculateRotatedSlipConditionAppliedSlipVariableDerivatives(
            entity_first_derivatives, aux_matrix, rCurrentEntity.GetGeometry());

        rCurrentEntity.CalculateSecondDerivativesLHS(aux_matrix, rProcessInfo);
        aux_matrix *= (1.0 - this->mBossak.Alpha);
        this->mpResponseFunction->CalculateSecondDerivativesGradient(
            rCurrentEntity, aux_matrix, response_second_derivatives, rProcessInfo);

        mAdjointSlipUtilities.CalculateRotatedSlipConditionAppliedNonSlipNonShapeVariableDerivatives(
            entity_second_derivatives, aux_matrix, rCurrentEntity.GetGeometry());

        if (rAdjointTimeSchemeValues2.size() != response_first_derivatives.size()) {
            rAdjointTimeSchemeValues2.resize(response_first_derivatives.size(), false);
        }

        noalias(rAdjointTimeSchemeValues2) =
            -response_first_derivatives - prod(entity_first_derivatives, adjoint_values);

        if (rAdjointTimeSchemeValues3.size() != response_second_derivatives.size()) {
            rAdjointTimeSchemeValues3.resize(response_second_derivatives.size(), false);
        }

        noalias(rAdjointTimeSchemeValues3) =
            -response_second_derivatives - prod(entity_second_derivatives, adjoint_values);

        KRATOS_CATCH("");
    }

    template <class TEntityType>
    void CalculateEntityAuxiliaryVariableContributions(
        TEntityType& rCurrentEntity,
        LocalSystemVectorType& rAdjointAuxiliaryValues,
        const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY

        const int k = OpenMPUtils::ThisThread();
        auto& adjoint_values = this->mAdjointValuesVector[k];
        auto& aux_matrix = this->mAuxiliaryMatrix[k];
        auto& entity_second_derivatives = this->mSecondDerivsLHS[k];
        auto& response_second_derivatives = this->mSecondDerivsResponseGradient[k];

        const auto& r_const_entity_ref = rCurrentEntity;
        r_const_entity_ref.GetValuesVector(adjoint_values);
        this->CheckAndResizeThreadStorage(adjoint_values.size());

        rCurrentEntity.CalculateSecondDerivativesLHS(aux_matrix, rProcessInfo);
        aux_matrix *= this->mBossak.Alpha;
        this->mpResponseFunction->CalculateSecondDerivativesGradient(
            rCurrentEntity, aux_matrix, response_second_derivatives, rProcessInfo);

        mAdjointSlipUtilities.CalculateRotatedSlipConditionAppliedNonSlipNonShapeVariableDerivatives(
            entity_second_derivatives, aux_matrix, rCurrentEntity.GetGeometry());

        if (rAdjointAuxiliaryValues.size() != entity_second_derivatives.size1()) {
            rAdjointAuxiliaryValues.resize(entity_second_derivatives.size1(), false);
        }
        noalias(rAdjointAuxiliaryValues) =
            prod(entity_second_derivatives, adjoint_values) + response_second_derivatives;

        KRATOS_CATCH("");
    }

    void CalculateResponseFunctionInterpolationError(ModelPart& rModelPart)
    {
        KRATOS_TRY

        const IndexType number_of_dofs_per_node = GetNumberOfDofsPerNode(rModelPart);
        const Vector& zero_vector_values = ZeroVector(number_of_dofs_per_node);

        // thread local storage
        struct TLS
        {
            // constructor
            TLS(ElementRefinementProcess::Pointer pElementRefinementProcess,
                IndexType NumberOfDofsPerNode)
                : mpElementRefinementProcess(pElementRefinementProcess),
                  mNumberOfDofsPerNode(NumberOfDofsPerNode)
            {
            }

            // copy constructor
            // this is called when creating TLS in each thread
            TLS(const TLS& rTLS)
                : mpElementRefinementProcess(rTLS.mpElementRefinementProcess),
                  mNumberOfDofsPerNode(rTLS.mNumberOfDofsPerNode)
            {
                pRefinedModelPart = &mpElementRefinementProcess->GetThreadLocalModelPart();

                for (auto& r_element : pRefinedModelPart->Elements()) {
                    r_element.EquationIdVector(EquationIds, pRefinedModelPart->GetProcessInfo());
                    for (IndexType i = 0; i < EquationIds.size(); ++i) {
                        const IndexType equation_id = EquationIds[i];
                        auto& r_data = AssemblyData[equation_id];
                        std::get<0>(r_data) = (i % mNumberOfDofsPerNode);
                    }
                }
            }

            // TLS data containers
            Matrix LHS;
            Matrix SecondDerivatives;
            Matrix RotatedSecondDerivatives;
            Vector Values;
            Vector AuxiliaryValues1;
            Vector AuxiliaryValues2;
            Vector AuxiliaryValues3;
            EquationIdVectorType EquationIds;

            // equation_id, <dof_index, aux_value_1, aux_value_2>
            AssemblyDataType AssemblyData;
            ModelPart* pRefinedModelPart;

            // data to be passed when TLS is copied
            const ElementRefinementProcess::Pointer mpElementRefinementProcess;
            const IndexType mNumberOfDofsPerNode;
        };

        block_for_each(rModelPart.Elements(), TLS(mpElementRefinementProcess, number_of_dofs_per_node), [&](Element& rCoarseElement, TLS& rTLS) {
            const auto& r_process_info = rTLS.pRefinedModelPart->GetProcessInfo();

            for (auto& p : rTLS.AssemblyData) {
                std::get<1>(p.second) = 0.0;
                std::get<2>(p.second) = 0.0;
            }

            auto p_refined_nodal_data_itr = mRefinedNodalData.find(rCoarseElement.Id());
            KRATOS_DEBUG_ERROR_IF(p_refined_nodal_data_itr == mRefinedNodalData.end()) << "Coarse element id " << rCoarseElement.Id() << " not found in the refined nodal data map.\n";
            auto& r_refined_nodal_data_map = p_refined_nodal_data_itr->second;

            mpElementRefinementProcess->InterpolateThreadLocalRefinedMeshFromCoarseElement(rCoarseElement);

            for (auto& r_refined_element : rTLS.pRefinedModelPart->Elements()) {
                r_refined_element.Initialize(r_process_info);
                r_refined_element.InitializeSolutionStep(r_process_info);
            }

            for (auto& r_condition : rTLS.pRefinedModelPart->Conditions()) {
                r_condition.Initialize(r_process_info);
                r_condition.InitializeSolutionStep(r_process_info);
            }

            for (auto& r_refined_element : rTLS.pRefinedModelPart->Elements()) {
                CalculateAndAssembleLocalEntityContributions(
                    r_refined_element, rTLS.EquationIds, rTLS.LHS, rTLS.SecondDerivatives, rTLS.RotatedSecondDerivatives,
                    rTLS.Values, rTLS.AuxiliaryValues1, rTLS.AuxiliaryValues2,
                    rTLS.AssemblyData, r_refined_nodal_data_map, r_process_info);
            }

            for (auto& r_refined_condition : rTLS.pRefinedModelPart->Conditions()) {
                if (r_refined_condition.Is(ACTIVE)) {
                    CalculateAndAssembleLocalEntityContributions(
                        r_refined_condition, rTLS.EquationIds, rTLS.LHS, rTLS.SecondDerivatives, rTLS.RotatedSecondDerivatives,
                        rTLS.Values, rTLS.AuxiliaryValues1, rTLS.AuxiliaryValues2,
                        rTLS.AssemblyData, r_refined_nodal_data_map, r_process_info);
                }
            }

            UpdateResponseFunctionInterpolationErrorTimeStepValues(
                *rTLS.pRefinedModelPart, rTLS.Values,
                rTLS.AuxiliaryValues1, rTLS.AuxiliaryValues2, rTLS.AuxiliaryValues3,
                rTLS.SecondDerivatives, rTLS.RotatedSecondDerivatives, rTLS.LHS, r_refined_nodal_data_map);

            // now add \mu^{n,h}_{mp} contributions from refined grid to b^{n,h}_{mp}
            for (auto& r_refined_element : rTLS.pRefinedModelPart->Elements()) {
                r_refined_element.EquationIdVector(rTLS.EquationIds, r_process_info);

                IndexType local_index = 0;
                for (auto& r_node : r_refined_element.GetGeometry()) {
                    const auto& r_values = r_refined_nodal_data_map.find(r_node.Id())->second.Mu;
                    for (IndexType i = 0; i < r_values.size(); ++i) {
                        std::get<2>(rTLS.AssemblyData[rTLS.EquationIds[local_index++]]) += r_values[i];
                    }
                }
            }

            auto& r_response_function_interpolation_auxiliary_1 = rCoarseElement.GetValue(RESPONSE_FUNCTION_INTERPOLATION_ERROR_AUXILIARY_1);
            auto& r_response_function_interpolation_auxiliary_2 = rCoarseElement.GetValue(RESPONSE_FUNCTION_INTERPOLATION_ERROR_AUXILIARY_2);

            noalias(r_response_function_interpolation_auxiliary_1) = zero_vector_values;
            noalias(r_response_function_interpolation_auxiliary_2) = zero_vector_values;

            for (auto& p : rTLS.AssemblyData) {
                const IndexType dof_id = std::get<0>(p.second);
                r_response_function_interpolation_auxiliary_1[dof_id] += std::abs(std::get<1>(p.second));
                r_response_function_interpolation_auxiliary_2[dof_id] += std::abs(std::get<2>(p.second));
            }
        });

        KRATOS_CATCH("");
    }

    // This method can be called in parallel regions since this is used inside open mp loops of original model part
    // and iterated over refined elements of the original model part elements
    template<class TEntityType>
    void CalculateAndAssembleLocalEntityContributions(
        TEntityType& rRefinedEntity,
        EquationIdVectorType& rEquationIdsVector,
        Matrix& rLHS,
        Matrix& rSecondDerivatives,
        Matrix& rRotatedSecondDerivatives,
        Vector& rValues,
        Vector& rAuxiliaryValues1,
        Vector& rAuxiliaryValues2,
        AssemblyDataType& rAssemblyData,
        const std::unordered_map<IndexType, RefinedNodeAdjointDataStorage>& rRefinedNodalMap,
        const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY

        rRefinedEntity.GetValuesVector(rValues);
        if (rAuxiliaryValues1.size() != rValues.size()) {
            rAuxiliaryValues1.resize(rValues.size());
        }

        if (rAuxiliaryValues2.size() != rValues.size()) {
            rAuxiliaryValues2.resize(rValues.size(), false);
        }

        this->CalculateLHSContribution(rRefinedEntity, rLHS, rEquationIdsVector, rProcessInfo);
        rLHS *= (1.0 / mBossak.C6);
        noalias(rAuxiliaryValues1) = prod(rLHS, rValues);

        rRefinedEntity.CalculateSecondDerivativesLHS(rSecondDerivatives, rProcessInfo);
        rSecondDerivatives *= (1.0 - mBossak.Alpha);

        mAdjointSlipUtilities.CalculateRotatedSlipConditionAppliedNonSlipNonShapeVariableDerivatives(
            rRotatedSecondDerivatives, rSecondDerivatives, rRefinedEntity.GetGeometry());

        noalias(rAuxiliaryValues2) = prod(rRotatedSecondDerivatives, rValues);

        this->AddAuxiliaryVariableContributionsForAdaptiveMeshRefinementFromPreviousTimeStep(
            rAuxiliaryValues1, rAuxiliaryValues2, rRefinedNodalMap, rRefinedEntity);

        for (IndexType i = 0; i < rValues.size(); ++i) {
            const auto equation_id = rEquationIdsVector[i];
            auto& r_data = rAssemblyData.find(equation_id)->second;
            std::get<1>(r_data) += rAuxiliaryValues1[i];
            std::get<2>(r_data) += rAuxiliaryValues2[i];
        }

        KRATOS_CATCH("");
    }

    void AddAuxiliaryVariableContributionsForAdaptiveMeshRefinementFromPreviousTimeStep(
        Vector& rAuxiliaryValues1,
        Vector& rAuxiliaryValues2,
        const std::unordered_map<IndexType, RefinedNodeAdjointDataStorage>& rRefinedNodalMap,
        const Element& rElement)
    {
        KRATOS_TRY

        IndexType local_index{0};
        for (const auto& r_node : rElement.GetGeometry()) {
            const auto& r_nodal_data = rRefinedNodalMap.find(r_node.Id())->second;
            const double weight = 1.0 / r_node.GetValue(NUMBER_OF_NEIGHBOUR_ELEMENTS);

            for (IndexType j = 0; j < r_nodal_data.E.size(); ++j){
                const double e = r_nodal_data.E[j];
                const double mu = r_nodal_data.Mu[j];
                rAuxiliaryValues1[local_index] += (mC1 * e + mC2 * mu) * weight;
                rAuxiliaryValues2[local_index] += (e - mC3 * mu) * weight;
                ++local_index;
            }
        }

        KRATOS_CATCH("");
    }

    void AddAuxiliaryVariableContributionsForAdaptiveMeshRefinementFromPreviousTimeStep(
        Vector& rAuxiliaryValues1,
        Vector& rAuxiliaryValues2,
        const std::unordered_map<IndexType, RefinedNodeAdjointDataStorage>& rRefinedNodalMap,
        const Condition& rCondition)
    {
        // do nothing here since there are no auxiliary contributions from conditions. Condition contributions
        // are also included in the element values
    }

    void UpdateResponseFunctionInterpolationErrorTimeStepValues(
        ModelPart& rThreadLocalModelPart,
        Vector& rValues,
        Vector& rAuxiliaryVector,
        Vector& rSecondDerivsResponseGradient,
        Vector& rSecondDerivsResponseGradientOld,
        Matrix& rSecondDerivsLHS,
        Matrix& rRotatedSecondDerivsLHS,
        Matrix& rAuxiliaryMatrix,
        std::unordered_map<IndexType, RefinedNodeAdjointDataStorage>& rRefinedNodalMap)
    {
        KRATOS_TRY

        const auto& r_process_info = rThreadLocalModelPart.GetProcessInfo();

        // first add all previous time step values for mu computation
        for (auto& r_map_pair : rRefinedNodalMap) {
            noalias(r_map_pair.second.Mu) = r_map_pair.second.Mu * mC3 -
                                            r_map_pair.second.D -
                                            r_map_pair.second.E;
        }

        // now clear maps
        for (auto& r_map_pair : rRefinedNodalMap) {
            r_map_pair.second.D.clear();
            r_map_pair.second.E.clear();
        }

        for (auto& r_refined_element : rThreadLocalModelPart.Elements()) {
            UpdateResponseFunctionInterpolationErrorTimeStepValuesFromEntity(
                r_process_info, rValues, rAuxiliaryVector, rSecondDerivsResponseGradient,
                rSecondDerivsResponseGradientOld, rSecondDerivsLHS, rRotatedSecondDerivsLHS,
                rAuxiliaryMatrix, rRefinedNodalMap, r_refined_element);
        }

        for (auto& r_refined_condition : rThreadLocalModelPart.Conditions()) {
            if (r_refined_condition.Is(ACTIVE)) {
                UpdateResponseFunctionInterpolationErrorTimeStepValuesFromEntity(
                    r_process_info, rValues, rAuxiliaryVector, rSecondDerivsResponseGradient,
                    rSecondDerivsResponseGradientOld, rSecondDerivsLHS, rRotatedSecondDerivsLHS,
                    rAuxiliaryMatrix, rRefinedNodalMap, r_refined_condition);
            }
        }

        KRATOS_CATCH("");
    }

    template<class EntityType>
    void UpdateResponseFunctionInterpolationErrorTimeStepValuesFromEntity(
        const ProcessInfo& rProcessInfo,
        Vector& rValues,
        Vector& rAuxiliaryVector,
        Vector& rSecondDerivsResponseGradient,
        Vector& rSecondDerivsResponseGradientOld,
        Matrix& rSecondDerivsLHS,
        Matrix& rRotatedSecondDerivsLHS,
        Matrix& rAuxiliaryMatrix,
        std::unordered_map<IndexType, RefinedNodeAdjointDataStorage>& rRefinedNodalMap,
        EntityType& rRefinedEntity)
    {
        KRATOS_TRY

        // get \lambda
        rRefinedEntity.GetValuesVector(rValues);
        if (rAuxiliaryVector.size() != rValues.size()) {
            rAuxiliaryVector.resize(rValues.size(), false);
        }

        // calculate \frac{\partial R^n}{\partial \dot{w}^n}
        rRefinedEntity.CalculateSecondDerivativesLHS(rSecondDerivsLHS, rProcessInfo);

        mAdjointSlipUtilities.CalculateRotatedSlipConditionAppliedNonSlipNonShapeVariableDerivatives(
            rRotatedSecondDerivsLHS, rSecondDerivsLHS, rRefinedEntity.GetGeometry());

        if (rAuxiliaryMatrix.size1() != rSecondDerivsLHS.size1() || rAuxiliaryMatrix.size2() != rSecondDerivsLHS.size2()) {
            rAuxiliaryMatrix.resize(rSecondDerivsLHS.size1(), rSecondDerivsLHS.size2(), false);
        }

        // calculate \frac{\partial J^n}{\partial \dot{w}^n}
        noalias(rAuxiliaryMatrix) = rSecondDerivsLHS * mC4;
        this->mpResponseFunction->CalculateSecondDerivativesGradient(rRefinedEntity, rAuxiliaryMatrix, rSecondDerivsResponseGradient, rProcessInfo);

        // calculate \frac{\partial J^n}{\partial \dot{w}^{n-1}}
        noalias(rAuxiliaryMatrix) = rSecondDerivsLHS * mBossak.Alpha;
        this->mpResponseFunction->CalculateSecondDerivativesGradient(rRefinedEntity, rAuxiliaryMatrix, rSecondDerivsResponseGradientOld, rProcessInfo);

        // calculate \lambda^n \frac{\partial R^n}{\partial \dot{w}^{n-1}}
        noalias(rAuxiliaryMatrix) = rRotatedSecondDerivsLHS * mBossak.Alpha;
        noalias(rAuxiliaryVector) = prod(rAuxiliaryMatrix, rValues);

        IndexType local_index{0};
        for (const auto& r_node : rRefinedEntity.GetGeometry()) {
            auto& r_nodal_data = rRefinedNodalMap.find(r_node.Id())->second;
            for (IndexType i = 0; i < r_nodal_data.D.size(); ++i) {
                const double temp = rAuxiliaryVector[local_index];

                // store \frac{\partial f^n}{\partial \dot{w}^{n-1}} in nodes
                r_nodal_data.D[i]  += rSecondDerivsResponseGradientOld[local_index];

                // store \lambda^n \frac{\partial R^n}{\partial \dot{w}^{n-1}} in nodes
                r_nodal_data.E[i]  += temp;

                // add -\frac{\partial f^n}{\partial \dot{w}^n} - \lambda^n\frac{\partial R^n}{\partial \dot{w}^n}
                r_nodal_data.Mu[i] -= (rSecondDerivsResponseGradient[local_index] + temp * mC5);
                ++local_index;
            }
        }

        KRATOS_CATCH("");
    }

    ///@}

}; /* Class VelocityBossakAdjointScheme */

///@}

} /* namespace Kratos.*/

#endif /* KRATOS_VELOCITY_BOSSAK_ADJOINT_SCHEME_H_INCLUDED defined */
