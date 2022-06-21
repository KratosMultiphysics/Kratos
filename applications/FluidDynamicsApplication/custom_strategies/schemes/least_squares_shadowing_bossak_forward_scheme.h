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

#if !defined(KRATOS_LSS_BOSSAK_FORWARD_SCHEME_H_INCLUDED)
#define KRATOS_LSS_BOSSAK_FORWARD_SCHEME_H_INCLUDED

// System includes
#include <string>

// External includes

// Project includes
#include "includes/define.h"
#include "solving_strategies/schemes/scheme.h"
#include "utilities/parallel_utilities.h"
#include "utilities/time_discretization.h"

// Application includes
#include "custom_utilities/fluid_adjoint_slip_utilities.h"
#include "custom_utilities/fluid_least_squares_shadowing_utilities.h"
#include "fluid_dynamics_application_variables.h"

namespace Kratos
{
///@name Kratos Classes
///@{

template <class TSparseSpace, class TDenseSpace>
class LeastSquaresShadowingBossakForwardScheme : public Scheme<TSparseSpace, TDenseSpace>
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(LeastSquaresShadowingBossakForwardScheme);

    using IndexType = std::size_t;

    using BaseType = Scheme<TSparseSpace, TDenseSpace>;

    using LocalSystemVectorType = typename BaseType::LocalSystemVectorType;

    using LocalSystemMatrixType = typename BaseType::LocalSystemMatrixType;

    using TSystemVectorType = typename BaseType::TSystemVectorType;

    using TSystemMatrixType = typename BaseType::TSystemMatrixType;

    using EquationIdVectorType = std::vector<IndexType>;

    using AssembledVectorType = std::unordered_map<IndexType, double>;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    explicit LeastSquaresShadowingBossakForwardScheme(
        const double BossakAlpha,
        const double DeltaTimeDialationAlpha,
        const IndexType Dimension,
        const IndexType BlockSize,
        const IndexType ShapeDerivativeNodeId,
        const IndexType ShapeDerivativeDirection,
        const FluidLeastSquaresShadowingUtilities& rFluidLeastSquaresShadowingUtilities,
        const Variable<double>& rDeltaTimeTotalShapeDerivativeVariable,
        const IndexType EchoLevel)
        : BaseType(),
          mDeltaTimeDialationAlpha(DeltaTimeDialationAlpha),
          mDimension(Dimension),
          mShapeDerivativeNodeId(ShapeDerivativeNodeId),
          mShapeDerivativeDirection(ShapeDerivativeDirection),
          mrDeltaTimeTotalShapeDerivativeVariable(rDeltaTimeTotalShapeDerivativeVariable),
          mBossak(BossakAlpha),
          mBossakConstants(mBossak),
          mFluidLeastSquaresShadowingUtilities(rFluidLeastSquaresShadowingUtilities),
          mEchoLevel(EchoLevel),
          mAdjointSlipUtilities(Dimension, BlockSize)
    {
        KRATOS_TRY

        // Allocate auxiliary memory.
        const int number_of_threads = ParallelUtilities::GetNumThreads();
        mTLS.resize(number_of_threads);

        KRATOS_ERROR_IF_NOT(mFluidLeastSquaresShadowingUtilities.GetPrimalVariablePointersList().size() == BlockSize)
                << "Provided block size does not match with the number of primal variables provided in least squares shadowing utilities. [ block size = "
                << BlockSize << ", number of primal variables in least squares shadowing utilities = " << mFluidLeastSquaresShadowingUtilities.GetPrimalVariablePointersList().size() << " ].\n";

        KRATOS_ERROR_IF(mShapeDerivativeDirection >= Dimension)
            << "Requested shape derivative direction index is larger than the problem dimension. [ shape derivative direction index = "
            << mShapeDerivativeDirection << ", problem dimension = " << Dimension << " ].\n";

        KRATOS_INFO_IF(this->Info(), mEchoLevel > 0) << "Created [ Dimensionality = " << Dimension << ", BlockSize = " << BlockSize << " ].\n";

        KRATOS_CATCH("");
    }

    /// Destructor.
    ~LeastSquaresShadowingBossakForwardScheme() override = default;

    ///@}
    ///@name Operations
    ///@{

    int Check(const ModelPart& rModelPart) const override
    {
        KRATOS_TRY

        IndexPartition<IndexType>(rModelPart.NumberOfElements()).for_each([&](const IndexType iElement) {
            const auto& r_element = *(rModelPart.ElementsBegin() + iElement);
            mFluidLeastSquaresShadowingUtilities.CheckVariables(r_element);
        });

        const IndexType buffer_size = rModelPart.GetBufferSize();
        KRATOS_ERROR_IF(buffer_size < 2) << "Buffer size needs to be greater than 1 in " << rModelPart.FullName() << " [ buffer size = " << buffer_size << " ].\n";

        const double delta_time = rModelPart.GetProcessInfo()[DELTA_TIME];
        KRATOS_ERROR_IF(delta_time < 0.0) << "LeastSquaresShadowingBossakForwardScheme should be run in forward time with positive delta time. [ delta_time = " << delta_time << " ].\n";

        return BaseType::Check(rModelPart);

        KRATOS_CATCH("");
    }

    void Initialize(ModelPart& rModelPart) override
    {
        KRATOS_TRY

        BaseType::Initialize(rModelPart);

        mAdjointSlipUtilities.Initialize(rModelPart);

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

        // Calculate the delta time total shape derivative
        struct TLS {
            Matrix mRotatedResidualPartialTimeDerivative;
            Matrix mResidualPartialTimeDerivative;
            Vector mAdjointSolution;
        };

        auto& r_process_info = rModelPart.GetProcessInfo();

        const double elemental_eta = block_for_each<SumReduction<double>>(rModelPart.Elements(), TLS(), [&](ModelPart::ElementType& rElement, TLS& rTLS) -> double {
            mFluidLeastSquaresShadowingUtilities.GetAdjointValues(rTLS.mAdjointSolution, rElement, 0);
            rElement.CalculateSensitivityMatrix(TIME_STEP_SENSITIVITY, rTLS.mResidualPartialTimeDerivative, r_process_info);

            mAdjointSlipUtilities.CalculateRotatedSlipConditionAppliedNonSlipNonShapeVariableDerivatives(
                rTLS.mRotatedResidualPartialTimeDerivative, rTLS.mResidualPartialTimeDerivative, rElement.GetGeometry());

            return inner_prod(rTLS.mAdjointSolution, row(rTLS.mRotatedResidualPartialTimeDerivative, 0));
        });

        const double delta_time = r_process_info[DELTA_TIME];
        const double coeff = delta_time * delta_time * mBossak.GetGamma();
        const auto& r_primal_variable_pointers_list = mFluidLeastSquaresShadowingUtilities.GetPrimalVariablePointersList();
        const auto& r_adjoint_first_derivative_variable_pointers_list = mFluidLeastSquaresShadowingUtilities.GetAdjointFirstDerivativeVariablePointersList();
        const IndexType number_of_variables = r_primal_variable_pointers_list.size();

        const double nodal_eta = block_for_each<SumReduction<double>>(rModelPart.GetCommunicator().LocalMesh().Nodes(), [&](ModelPart::NodeType& rNode) -> double {
            double value = 0.0;
            for (IndexType i = 0; i < number_of_variables; ++i) {
                value += (rNode.FastGetSolutionStepValue(*r_primal_variable_pointers_list[i], 0) - rNode.FastGetSolutionStepValue(*r_primal_variable_pointers_list[i], 1)) * rNode.FastGetSolutionStepValue(*r_adjoint_first_derivative_variable_pointers_list[i]);
            }
            return value;
        });

        const double eta = elemental_eta + nodal_eta / coeff;
        r_process_info[mrDeltaTimeTotalShapeDerivativeVariable] = rModelPart.GetCommunicator().GetDataCommunicator().SumAll(eta) / (-2.0 * mDeltaTimeDialationAlpha * mDeltaTimeDialationAlpha);

        // set bossak constants
        mBossakConstants.SetConstants(delta_time);

        KRATOS_INFO_IF(this->Info(), mEchoLevel > 0) << "Computed time dialation total shape derivative in " << rModelPart.FullName()
            << " and stored in " << mrDeltaTimeTotalShapeDerivativeVariable.Name() << ".\n";

        KRATOS_CATCH("")
    }

    void CalculateSystemContributions(
        Element& rCurrentElement,
        LocalSystemMatrixType& rLHS_Contribution,
        LocalSystemVectorType& rRHS_Contribution,
        Element::EquationIdVectorType& rEquationId,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        CalculateEntitySystemContributions(rCurrentElement, rLHS_Contribution, rRHS_Contribution,
                                           rEquationId, rCurrentProcessInfo);
    }

    void CalculateLHSContribution(
        Element& rCurrentElement,
        LocalSystemMatrixType& rLHS_Contribution,
        Element::EquationIdVectorType& rEquationId,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        const int thread_id = OpenMPUtils::ThisThread();
        TLS& r_tls = mTLS[thread_id];

        CalculateEntityLHSContribution(
            rLHS_Contribution, rEquationId, rCurrentElement, r_tls.mResidualFirstDerivatives, r_tls.mRotatedResidualFirstDerivatives,
            r_tls.mResidualSecondDerivatives, r_tls.mRotatedResidualSecondDerivatives, rCurrentProcessInfo);
    }

    void CalculateSystemContributions(
        Condition& rCurrentCondition,
        LocalSystemMatrixType& rLHS_Contribution,
        LocalSystemVectorType& rRHS_Contribution,
        Condition::EquationIdVectorType& rEquationId,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        CalculateEntitySystemContributions(rCurrentCondition, rLHS_Contribution, rRHS_Contribution,
                                           rEquationId, rCurrentProcessInfo);
    }

    void CalculateLHSContribution(
        Condition& rCurrentCondition,
        LocalSystemMatrixType& rLHS_Contribution,
        Condition::EquationIdVectorType& rEquationId,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        const int thread_id = OpenMPUtils::ThisThread();
        TLS& r_tls = mTLS[thread_id];

        CalculateEntityLHSContribution(
            rLHS_Contribution, rEquationId, rCurrentCondition, r_tls.mResidualFirstDerivatives, r_tls.mRotatedResidualFirstDerivatives,
            r_tls.mResidualSecondDerivatives, r_tls.mRotatedResidualSecondDerivatives, rCurrentProcessInfo);

    }

    void FinalizeSolutionStep(
        ModelPart& rModelPart,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b) override
    {
        KRATOS_TRY

        BaseType::FinalizeSolutionStep(rModelPart, A, Dx, b);

        // update the vdot
        auto& r_process_info = rModelPart.GetProcessInfo();
        const double delta_time = r_process_info[DELTA_TIME];
        const double coeff_1 = 1 / (delta_time * mBossak.GetGamma());
        const double coeff_2 = (mBossak.GetGamma() - 1) / mBossak.GetGamma();
        const double coeff_3 = r_process_info[mrDeltaTimeTotalShapeDerivativeVariable] / (delta_time * delta_time * mBossak.GetGamma());

        const auto& r_primal_variable_pointers_list = mFluidLeastSquaresShadowingUtilities.GetPrimalVariablePointersList();
        const auto& r_lss_variable_pointers_list = mFluidLeastSquaresShadowingUtilities.GetLSSVariablePointersList();
        const auto& r_lss_first_derivative_variable_pointers_list = mFluidLeastSquaresShadowingUtilities.GetLSSFirstDerivativeVariablePointersList();
        const IndexType number_of_variables = r_primal_variable_pointers_list.size();

        block_for_each(rModelPart.Nodes(), [&](ModelPart::NodeType& rNode) {
            for (IndexType i = 0; i < number_of_variables; ++i){
                double value = 0.0;

                value += coeff_1 * rNode.FastGetSolutionStepValue(*r_lss_variable_pointers_list[i]);
                value -= coeff_1 * rNode.FastGetSolutionStepValue(*r_lss_variable_pointers_list[i], 1);
                value += coeff_2 * rNode.FastGetSolutionStepValue(*r_lss_first_derivative_variable_pointers_list[i], 1);
                value -= coeff_3 * rNode.FastGetSolutionStepValue(*r_primal_variable_pointers_list[i]);
                value += coeff_3 * rNode.FastGetSolutionStepValue(*r_primal_variable_pointers_list[i], 1);

                rNode.FastGetSolutionStepValue(*r_lss_first_derivative_variable_pointers_list[i]) = value;
            }
        });

        KRATOS_INFO_IF(this->Info(), mEchoLevel > 0) << "Computed primal first derivative total shape derivatives in " << rModelPart.FullName() << ".\n";

        KRATOS_CATCH("")
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "LeastSquaresShadowingBossakForwardScheme";
    }

    ///@}

protected:
    ///@name Protected Member Variables
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
private:
    ///@name Private Classes
    ///@{

    struct BossakConstants
    {
        double mC1;
        double mC2;
        double mC3;
        double mC4;
        const TimeDiscretization::Bossak& mrBossak;

        explicit BossakConstants(const TimeDiscretization::Bossak& rBossak) : mrBossak(rBossak)
        {
        }

        void SetConstants(const double DeltaTime)
        {
            mC1 = 1 / (DeltaTime * mrBossak.GetGamma());
            mC2 = (mrBossak.GetGamma() - 1 + mrBossak.GetAlphaM()) / mrBossak.GetGamma();
            mC3 = (1 - mrBossak.GetAlphaM()) / (DeltaTime * mrBossak.GetGamma());
            mC4 = (1 - mrBossak.GetAlphaM()) / (DeltaTime * DeltaTime * mrBossak.GetGamma());
        }
    };

    struct TLS
    {
        Vector mLSSValues;

        // LHS TLS
        Matrix mResidualFirstDerivatives;
        Matrix mRotatedResidualFirstDerivatives;
        Matrix mResidualSecondDerivatives;
        Matrix mRotatedResidualSecondDerivatives;

        // RHS TLS
        Vector mCurrentPrimalSolution;
        Vector mPreviousPrimalSolution;
        Vector mPreviousLSSSolution;
        Vector mPreviousLSSSecondDerivativeSolution;
        Vector mResiduals;
        Matrix mAuxMatrix;
        Matrix mMassMatrix;
        Matrix mRotatedMassMatrix;
        Matrix mResidualTimeStepDerivatives;
        Matrix mRotatedResidualTimeStepDerivatives;
        Matrix mResidualShapeDerivatives;
        Matrix mRotatedResidualShapeDerivatives;
        std::vector<IndexType> mDerivativeNodeIds;
    };

    ///@}
    ///@name Private Member Variables
    ///@{

    const double mDeltaTimeDialationAlpha;

    const IndexType mDimension;

    const IndexType mShapeDerivativeNodeId;

    const IndexType mShapeDerivativeDirection;

    const Variable<double>& mrDeltaTimeTotalShapeDerivativeVariable;

    const TimeDiscretization::Bossak mBossak;

    BossakConstants mBossakConstants;

    const FluidLeastSquaresShadowingUtilities mFluidLeastSquaresShadowingUtilities;

    const IndexType mEchoLevel;

    FluidAdjointSlipUtilities mAdjointSlipUtilities;

    std::vector<TLS> mTLS;

    ///@}
    ///@name Private Operations
    ///@{

    template<class TEntityType>
    void CalculateEntitySystemContributions(
        TEntityType& rEntity,
        LocalSystemMatrixType& rLHS,
        LocalSystemVectorType& rRHS,
        typename TEntityType::EquationIdVectorType& rEquationId,
        const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY;

        const int thread_id = OpenMPUtils::ThisThread();
        TLS& r_tls = mTLS[thread_id];

        CalculateEntityLHSContribution<TEntityType>(
            rLHS, rEquationId, rEntity, r_tls.mResidualFirstDerivatives, r_tls.mRotatedResidualFirstDerivatives,
            r_tls.mResidualSecondDerivatives, r_tls.mRotatedResidualSecondDerivatives, rCurrentProcessInfo);

        CalculateEntityRHSContribution<TEntityType>(
            rRHS, rEntity, r_tls.mCurrentPrimalSolution, r_tls.mPreviousPrimalSolution, r_tls.mPreviousLSSSolution,
            r_tls.mPreviousLSSSecondDerivativeSolution, r_tls.mResiduals, r_tls.mAuxMatrix, r_tls.mMassMatrix,
            r_tls.mRotatedMassMatrix, r_tls.mResidualTimeStepDerivatives, r_tls.mRotatedResidualTimeStepDerivatives,
            r_tls.mResidualShapeDerivatives, r_tls.mRotatedResidualShapeDerivatives, r_tls.mDerivativeNodeIds, rCurrentProcessInfo);

        // Calculate system contributions in residual form.
        if (rLHS.size1() != 0) {
            rEntity.GetValuesVector(r_tls.mLSSValues);
            noalias(rRHS) -= prod(rLHS, r_tls.mLSSValues);
        }

        KRATOS_CATCH("");
    }

    template<class TEntityType>
    void CalculateEntityLHSContribution(
        Matrix& rLHS,
        typename TEntityType::EquationIdVectorType& rEquationId,
        TEntityType& rEntity,
        Matrix& rEntityResidualFirstDerivatives,
        Matrix& rEntityRotatedResidualFirstDerivatives,
        Matrix& rEntityResidualSecondDerivatives,
        Matrix& rEntityRotatedResidualSecondDerivatives,
        const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        rEntity.EquationIdVector(rEquationId, rCurrentProcessInfo);

        rEntity.CalculateFirstDerivativesLHS(rEntityResidualFirstDerivatives, rCurrentProcessInfo);
        mAdjointSlipUtilities.CalculateRotatedSlipConditionAppliedSlipVariableDerivatives(
            rEntityRotatedResidualFirstDerivatives, rEntityResidualFirstDerivatives, rEntity.GetGeometry());

        rEntity.CalculateSecondDerivativesLHS(rEntityResidualSecondDerivatives, rCurrentProcessInfo);
        rEntityResidualSecondDerivatives *= (1.0 - this->mBossak.GetAlphaM());
        mAdjointSlipUtilities.CalculateRotatedSlipConditionAppliedNonSlipNonShapeVariableDerivatives(
            rEntityRotatedResidualSecondDerivatives, rEntityResidualSecondDerivatives, rEntity.GetGeometry());

        if (rLHS.size1() != rEntityRotatedResidualFirstDerivatives.size1() || rLHS.size2() != rEntityRotatedResidualFirstDerivatives.size2()) {
            rLHS.resize(rEntityRotatedResidualFirstDerivatives.size1(), rEntityRotatedResidualFirstDerivatives.size2(), false);
        }

        noalias(rLHS) = trans(rEntityRotatedResidualFirstDerivatives + mBossakConstants.mC1 * rEntityRotatedResidualSecondDerivatives);

        KRATOS_CATCH("");
    }

    template<class TEntityType>
    void CalculateEntityRHSContribution(
        Vector& rRHS,
        TEntityType& rEntity,
        Vector& rCurrentPrimalSolution,
        Vector& rPreviousPrimalSolution,
        Vector& rPreviousLSSSolution,
        Vector& rPreviousLSSSecondDerivativeSolution,
        Vector& rResiduals,
        Matrix& rAuxMatrix,
        Matrix& rMassMatrix,
        Matrix& rRotatedMassMatrix,
        Matrix& rResidualTimeStepDerivatives,
        Matrix& rRotatedResidualTimeStepDerivatives,
        Matrix& rResidualShapeDerivatives,
        Matrix& rRotatedResidualShapeDerivatives,
        std::vector<IndexType>& rDerivativeNodeIds,
        const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        const double eta = rCurrentProcessInfo[mrDeltaTimeTotalShapeDerivativeVariable];

        // In LSS methods, adjoint elements needs to implement their Mass matrices.
        // They cannot be zero.
        rEntity.CalculateMassMatrix(rMassMatrix, rCurrentProcessInfo);
        mAdjointSlipUtilities.CalculateRotatedSlipConditionAppliedNonSlipNonShapeVariableDerivatives(
            rRotatedMassMatrix, rMassMatrix, rEntity.GetGeometry());

        // Calculates primal residual
        rEntity.CalculateLocalSystem(rAuxMatrix, rResiduals, rCurrentProcessInfo);
        rEntity.CalculateLocalVelocityContribution(rAuxMatrix, rResiduals, rCurrentProcessInfo);

        rEntity.CalculateSensitivityMatrix(SHAPE_SENSITIVITY, rResidualShapeDerivatives, rCurrentProcessInfo);
        mAdjointSlipUtilities.CalculateRotatedSlipConditionAppliedShapeVariableDerivatives(
            rRotatedResidualShapeDerivatives, rDerivativeNodeIds, rResiduals, rResidualShapeDerivatives, rEntity, rCurrentProcessInfo);

        rEntity.CalculateSensitivityMatrix(TIME_STEP_SENSITIVITY, rResidualTimeStepDerivatives, rCurrentProcessInfo);
        mAdjointSlipUtilities.CalculateRotatedSlipConditionAppliedNonSlipNonShapeVariableDerivatives(
            rRotatedResidualTimeStepDerivatives, rResidualTimeStepDerivatives, rEntity.GetGeometry());

        rEntity.GetValuesVector(rPreviousLSSSolution, 1);
        rEntity.GetSecondDerivativesVector(rPreviousLSSSecondDerivativeSolution, 1);
        mFluidLeastSquaresShadowingUtilities.GetPrimalValues(rCurrentPrimalSolution, rEntity, 0);
        mFluidLeastSquaresShadowingUtilities.GetPrimalValues(rPreviousPrimalSolution, rEntity, 1);

        if (rRHS.size() != rMassMatrix.size1()) {
            rRHS.resize(rMassMatrix.size1());
        }

        noalias(rRHS) = prod(rRotatedMassMatrix, rPreviousLSSSecondDerivativeSolution) * mBossakConstants.mC2;
        noalias(rRHS) = rRHS - rPreviousLSSSolution * mBossakConstants.mC3;
        noalias(rRHS) = rRHS - prod(rRotatedMassMatrix, rCurrentPrimalSolution - rPreviousPrimalSolution) * mBossakConstants.mC4 * eta;
        noalias(rRHS) = rRHS - row(rResidualTimeStepDerivatives, 0) * eta;

        const auto& p_itr = find(rDerivativeNodeIds.begin(), rDerivativeNodeIds.end(), mShapeDerivativeNodeId);
        if (p_itr != rDerivativeNodeIds.end()) {
            noalias(rRHS) = rRHS - row(rRotatedResidualShapeDerivatives, std::distance(rDerivativeNodeIds.begin(),  p_itr) * mDimension + mShapeDerivativeDirection);
        }

        KRATOS_CATCH("");
    }

    ///@}

}; /* Class LeastSquaresShadowingBossakForwardScheme */

///@}

} /* namespace Kratos.*/

#endif /* KRATOS_LSS_BOSSAK_FORWARD_SCHEME_H_INCLUDED defined */
