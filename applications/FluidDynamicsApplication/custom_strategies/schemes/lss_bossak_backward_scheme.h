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

#if !defined(KRATOS_LSS_BOSSAK_BACKWARD_SCHEME_H_INCLUDED)
#define KRATOS_LSS_BOSSAK_BACKWARD_SCHEME_H_INCLUDED

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
#include "custom_utilities/fluid_lss_variable_utilities.h"
#include "fluid_dynamics_application_variables.h"

namespace Kratos
{
///@name Kratos Classes
///@{

template <class TSparseSpace, class TDenseSpace>
class LSSBossakBackwardScheme : public Scheme<TSparseSpace, TDenseSpace>
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(LSSBossakBackwardScheme);

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
    explicit LSSBossakBackwardScheme(
        FluidLSSVariableUtilities::Pointer pFluidLeastSquaresShadowingVariableUtilities,
        const Variable<Vector>& rAuxVariable,
        const double BossakAlpha,
        const IndexType Dimension,
        const IndexType BlockSize,
        const IndexType EchoLevel)
        : BaseType(),
          mpFluidLeastSquaresShadowingVariableUtilities(pFluidLeastSquaresShadowingVariableUtilities),
          mrAuxVariable(rAuxVariable),
          mDimension(Dimension),
          mEchoLevel(EchoLevel),
          mBossak(BossakAlpha),
          mBossakConstants(mBossak),
          mAdjointSlipUtilities(Dimension, BlockSize)
    {
        KRATOS_TRY

        KRATOS_ERROR_IF_NOT(mpFluidLeastSquaresShadowingVariableUtilities->GetPrimalIndirectVariablesList().size() == BlockSize)
                << "Provided block size does not match with the number of primal variables provided in least squares shadowing utilities. [ block size = "
                << BlockSize << ", number of primal variables in least squares shadowing utilities = " << mpFluidLeastSquaresShadowingVariableUtilities->GetPrimalIndirectVariablesList().size() << " ].\n";

        KRATOS_INFO_IF(this->Info(), mEchoLevel > 0) << "Created [ Dimensionality = " << Dimension << ", BlockSize = " << BlockSize << " ].\n";

        KRATOS_CATCH("");
    }

    /// Destructor.
    ~LSSBossakBackwardScheme() override = default;

    ///@}
    ///@name Operations
    ///@{

    int Check(const ModelPart& rModelPart) const override
    {
        KRATOS_TRY

        const IndexType buffer_size = rModelPart.GetBufferSize();
        KRATOS_ERROR_IF(buffer_size < 2) << "Buffer size needs to be greater than 1 in " << rModelPart.FullName() << " [ buffer size = " << buffer_size << " ].\n";

        const double delta_time = rModelPart.GetProcessInfo()[DELTA_TIME];
        KRATOS_ERROR_IF(delta_time > 0.0) << "LSSBossakBackwardScheme should be run in backward time with negative delta time. [ delta_time = " << delta_time << " ].\n";

        return BaseType::Check(rModelPart);

        KRATOS_CATCH("");
    }

    void Initialize(ModelPart& rModelPart) override
    {
        KRATOS_TRY

        BaseType::Initialize(rModelPart);

        block_for_each(rModelPart.Elements(), Vector(), [&](ModelPart::ElementType& rElement, Vector& rTLS) {
            rElement.GetValuesVector(rTLS);
            noalias(rTLS) = ZeroVector(rTLS.size());
            rElement.SetValue(mrAuxVariable, rTLS);
        });

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

        auto& r_process_info = rModelPart.GetProcessInfo();

        // set bossak constants
        mBossakConstants.SetConstants(-r_process_info[DELTA_TIME]);

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
        TLS& r_tls = mSMPStorage;

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
        TLS& r_tls = mSMPStorage;

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

        const auto& r_process_info = rModelPart.GetProcessInfo();
        const double coeff_1 = (mBossak.GetGamma() - 1) / mBossak.GetGamma();

        // update mu
        const auto& r_adjoint_variable_pointers_list = mpFluidLeastSquaresShadowingVariableUtilities->GetAdjointIndirectVariablesList();
        const auto& r_adjoint_first_derivative_variable_pointers_list = mpFluidLeastSquaresShadowingVariableUtilities->GetAdjointFirstDerivativeIndirectVariablesList();
        const IndexType number_of_variables = r_adjoint_variable_pointers_list.size();

        // clear current first derivative variables
        block_for_each(rModelPart.Nodes(), [&](ModelPart::NodeType& rNode) {
            for (IndexType i = 0; i < number_of_variables; ++i) {
                r_adjoint_first_derivative_variable_pointers_list[i](rNode) = 0.0;
            }
        });

        struct TLS
        {
            Vector mAdjointValues;
            Vector mCurrentAux;
            Vector mPreviousAux;
            Matrix mResidualSecondDerivatives;
            Matrix mRotatedResidualSecondDerivatives;
        };

        // add the elemental contributions first
        block_for_each(rModelPart.Elements(), TLS(), [&](ModelPart::ElementType& rElement, TLS& rTLS) {
            auto& r_geometry = rElement.GetGeometry();

            rElement.GetValuesVector(rTLS.mAdjointValues, r_process_info);
            rElement.CalculateSecondDerivativesLHS(rTLS.mResidualSecondDerivatives, r_process_info);
            mAdjointSlipUtilities.CalculateRotatedSlipConditionAppliedNonSlipNonShapeVariableDerivatives(rTLS.mRotatedResidualSecondDerivatives, rTLS.mResidualSecondDerivatives, r_geometry);

            if (rTLS.mPreviousAux.size() != rTLS.mAdjointValues.size()) {
                rTLS.mCurrentAux.resize(rTLS.mAdjointValues.size());
                rTLS.mPreviousAux.resize(rTLS.mAdjointValues.size());
            }

            noalias(rTLS.mPreviousAux) = rElement.GetValue(mrAuxVariable);
            noalias(rTLS.mCurrentAux) = prod(rTLS.mRotatedResidualSecondDerivatives, rTLS.mAdjointValues);
            noalias(rElement.GetValue(mrAuxVariable)) = rTLS.mCurrentAux * this->mBossakConstants.mPreviousStepSecondDerivativeCoefficient;
            noalias(rTLS.mCurrentAux) = rTLS.mCurrentAux * this->mBossakConstants.mCurrentStepSecondDerivativeCoefficient +  rTLS.mPreviousAux;

            IndexType local_index = 0;
            for (IndexType i = 0; i < r_geometry.PointsNumber(); ++i) {
                auto& r_node = r_geometry[i];

                r_node.SetLock();
                for (IndexType j = 0; j < number_of_variables; ++j) {
                    r_adjoint_first_derivative_variable_pointers_list[j](r_node) -= rTLS.mCurrentAux[local_index++];
                }
                r_node.UnSetLock();
            }
        });

        for (IndexType i = 0; i < number_of_variables; ++i) {
            r_adjoint_first_derivative_variable_pointers_list[i].ExecuteForVariable([&](const Variable<double>& rVariable){
                rModelPart.GetCommunicator().AssembleCurrentData(rVariable);
            });
        }

        // add nodal contributions
        block_for_each(rModelPart.Nodes(), [&](ModelPart::NodeType& rNode) {
            for (IndexType i = 0; i < number_of_variables; ++i){
                r_adjoint_first_derivative_variable_pointers_list[i](rNode) += r_adjoint_first_derivative_variable_pointers_list[i](rNode, 1) * coeff_1;
            }
        });

        KRATOS_INFO_IF(this->Info(), mEchoLevel > 0) << "Computed adjoint first derivatives in " << rModelPart.FullName() << ".\n";

        KRATOS_CATCH("")
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "LSSBossakBackwardScheme";
    }

    ///@}

private:
    ///@name Private Classes
    ///@{

    struct BossakConstants
    {
        double mC1;
        double mC2;
        const TimeDiscretization::Bossak& mrBossak;
        double mCurrentStepSecondDerivativeCoefficient;
        double mPreviousStepSecondDerivativeCoefficient;

        explicit BossakConstants(const TimeDiscretization::Bossak& rBossak) : mrBossak(rBossak)
        {
            mCurrentStepSecondDerivativeCoefficient = 1 - mrBossak.GetAlphaM();
            mPreviousStepSecondDerivativeCoefficient = mrBossak.GetAlphaM();
        }

        void SetConstants(const double DeltaTime)
        {
            mC1 = 1 / (DeltaTime * mrBossak.GetGamma());
            mC2 = 1 / (DeltaTime * mrBossak.GetGamma() * mrBossak.GetGamma());
        }
    };

    struct TLS
    {
        Vector mAdjointValues;

        // LHS TLS
        Matrix mResidualFirstDerivatives;
        Matrix mRotatedResidualFirstDerivatives;
        Matrix mResidualSecondDerivatives;
        Matrix mRotatedResidualSecondDerivatives;

        // RHS TLS
        Vector rPreviousSecondDerivativesValues;
    };

    ///@}
    ///@name Private Member Variables
    ///@{

    FluidLSSVariableUtilities::Pointer mpFluidLeastSquaresShadowingVariableUtilities;

    const Variable<Vector>& mrAuxVariable;

    const IndexType mDimension;

    const IndexType mEchoLevel;

    const TimeDiscretization::Bossak mBossak;

    BossakConstants mBossakConstants;

    const FluidAdjointSlipUtilities mAdjointSlipUtilities;

    static TLS mSMPStorage;
    #pragma omp threadprivate(mSMPStorage)

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

        TLS& r_tls = mSMPStorage;

        CalculateEntityLHSContribution<TEntityType>(
            rLHS, rEquationId, rEntity, r_tls.mResidualFirstDerivatives, r_tls.mRotatedResidualFirstDerivatives,
            r_tls.mResidualSecondDerivatives, r_tls.mRotatedResidualSecondDerivatives, rCurrentProcessInfo);

        CalculateEntityRHSContribution<TEntityType>(
            rRHS, rEntity, r_tls.rPreviousSecondDerivativesValues);

        // Calculate system contributions in residual form.
        if (rLHS.size1() != 0) {
            rEntity.GetValuesVector(r_tls.mAdjointValues);
            noalias(rRHS) -= prod(rLHS, r_tls.mAdjointValues);
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
        mAdjointSlipUtilities.CalculateRotatedSlipConditionAppliedNonSlipNonShapeVariableDerivatives(
            rEntityRotatedResidualSecondDerivatives, rEntityResidualSecondDerivatives, rEntity.GetGeometry());

        if (rLHS.size1() != rEntityRotatedResidualFirstDerivatives.size1() || rLHS.size2() != rEntityRotatedResidualFirstDerivatives.size2()) {
            rLHS.resize(rEntityRotatedResidualFirstDerivatives.size1(), rEntityRotatedResidualFirstDerivatives.size2(), false);
        }

        noalias(rLHS) = rEntityRotatedResidualFirstDerivatives + rEntityRotatedResidualSecondDerivatives * (mBossakConstants.mC1 * mBossakConstants.mCurrentStepSecondDerivativeCoefficient);

        KRATOS_CATCH("");
    }

    template<class TEntityType>
    void CalculateEntityRHSContribution(
        Vector& rRHS,
        TEntityType& rEntity,
        Vector& rPreviousSecondDerivativesValues)
    {
        KRATOS_TRY

        rEntity.GetSecondDerivativesVector(rPreviousSecondDerivativesValues, 1);

        mpFluidLeastSquaresShadowingVariableUtilities->GetLSSValues(rRHS, rEntity);
        rRHS *= -2.0;
        noalias(rRHS) -= rEntity.GetValue(mrAuxVariable) * mBossakConstants.mC1;
        noalias(rRHS) -= rPreviousSecondDerivativesValues * mBossakConstants.mC2;

        KRATOS_CATCH("");
    }

    ///@}

}; /* Class LSSBossakBackwardScheme */

///@}

} /* namespace Kratos.*/

#endif /* KRATOS_LSS_BOSSAK_BACKWARD_SCHEME_H_INCLUDED defined */
