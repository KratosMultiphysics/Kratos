// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Richard Faasse
//

#include "backward_euler_quasistatic_U_Pw_scheme.hpp"

namespace Kratos
{
template <class TSparseSpace, class TDenseSpace>
class PrototypeNullSteppingScheme : public BackwardEulerQuasistaticUPwScheme<TSparseSpace, TDenseSpace>
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(PrototypeNullSteppingScheme);

    void CalculateSystemContributions(
        Element& rCurrentElement,
        typename GeoMechanicsTimeIntegrationScheme<TSparseSpace, TDenseSpace>::LocalSystemMatrixType& LHS_Contribution,
        typename GeoMechanicsTimeIntegrationScheme<TSparseSpace, TDenseSpace>::LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        const ProcessInfo&             CurrentProcessInfo) override
    {
        rCurrentElement.CalculateLeftHandSide(LHS_Contribution, CurrentProcessInfo);

        rCurrentElement.Calculate(INTERNAL_FORCES_VECTOR, RHS_Contribution, CurrentProcessInfo);
        const auto fraction_of_unbalance = (CurrentProcessInfo[TIME] - CurrentProcessInfo[START_TIME]) /
                                           (CurrentProcessInfo[END_TIME] - CurrentProcessInfo[START_TIME]);

        const auto f_ext = -mInternalForcesAtStartbyElementId[rCurrentElement.GetId()] +
                           (mExternalForcesAtStartbyElementId[rCurrentElement.GetId()] +
                            mInternalForcesAtStartbyElementId[rCurrentElement.GetId()]) *
                               fraction_of_unbalance;

        RHS_Contribution += f_ext;

        rCurrentElement.EquationIdVector(EquationId, CurrentProcessInfo);
    }

    void CalculateRHSContribution(Element& rCurrentElement,
                                  typename GeoMechanicsTimeIntegrationScheme<TSparseSpace, TDenseSpace>::LocalSystemVectorType& RHS_Contribution,
                                  Element::EquationIdVectorType& EquationId,
                                  const ProcessInfo&             CurrentProcessInfo) override
    {
        rCurrentElement.Calculate(INTERNAL_FORCES_VECTOR, RHS_Contribution, CurrentProcessInfo);
        const auto fraction_of_unbalance = (CurrentProcessInfo[TIME] - CurrentProcessInfo[START_TIME]) /
                                           (CurrentProcessInfo[END_TIME] - CurrentProcessInfo[START_TIME]);

        const auto f_ext = -mInternalForcesAtStartbyElementId[rCurrentElement.GetId()] +
                           (mExternalForcesAtStartbyElementId[rCurrentElement.GetId()] +
                            mInternalForcesAtStartbyElementId[rCurrentElement.GetId()]) *
                               fraction_of_unbalance;

        RHS_Contribution += f_ext;
        rCurrentElement.EquationIdVector(EquationId, CurrentProcessInfo);
    }

    void CalculateSystemContributions(
        Condition& rCurrentCondition,
        typename GeoMechanicsTimeIntegrationScheme<TSparseSpace, TDenseSpace>::LocalSystemMatrixType& LHS_Contribution,
        typename GeoMechanicsTimeIntegrationScheme<TSparseSpace, TDenseSpace>::LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        const ProcessInfo&             CurrentProcessInfo) override
    {
        rCurrentCondition.CalculateLocalSystem(LHS_Contribution, RHS_Contribution, CurrentProcessInfo);

        const auto fraction_of_unbalance = (CurrentProcessInfo[TIME] - CurrentProcessInfo[START_TIME]) /
                                           (CurrentProcessInfo[END_TIME] - CurrentProcessInfo[START_TIME]);

        RHS_Contribution *= fraction_of_unbalance;

        rCurrentCondition.EquationIdVector(EquationId, CurrentProcessInfo);
    }

    void CalculateRHSContribution(Condition& rCurrentCondition,
                                  typename GeoMechanicsTimeIntegrationScheme<TSparseSpace, TDenseSpace>::LocalSystemVectorType& RHS_Contribution,
                                  Element::EquationIdVectorType& EquationId,
                                  const ProcessInfo&             CurrentProcessInfo) override
    {
        rCurrentCondition.CalculateRightHandSide(RHS_Contribution, CurrentProcessInfo);
        const auto fraction_of_unbalance = (CurrentProcessInfo[TIME] - CurrentProcessInfo[START_TIME]) /
                                           (CurrentProcessInfo[END_TIME] - CurrentProcessInfo[START_TIME]);
        RHS_Contribution *= fraction_of_unbalance;
        rCurrentCondition.EquationIdVector(EquationId, CurrentProcessInfo);
    }

    void InitializeSolutionStep(
        ModelPart& rModelPart,
        typename GeoMechanicsTimeIntegrationScheme<TSparseSpace, TDenseSpace>::TSystemMatrixType& A,
        typename GeoMechanicsTimeIntegrationScheme<TSparseSpace, TDenseSpace>::TSystemVectorType& dX,
        typename GeoMechanicsTimeIntegrationScheme<TSparseSpace, TDenseSpace>::TSystemVectorType& b) override
    {
        BackwardEulerQuasistaticUPwScheme<TSparseSpace, TDenseSpace>::InitializeSolutionStep(
            rModelPart, A, dX, b);

        if (!mIsInitialized) {
            for (auto& rElement : rModelPart.Elements()) {
                mExternalForcesAtStartbyElementId.insert({rElement.GetId(), Vector{}});
                rElement.Calculate(EXTERNAL_FORCES_VECTOR,
                                   mExternalForcesAtStartbyElementId[rElement.GetId()],
                                   rModelPart.GetProcessInfo());

                mInternalForcesAtStartbyElementId.insert({rElement.GetId(), Vector{}});
                rElement.Calculate(INTERNAL_FORCES_VECTOR,
                                   mInternalForcesAtStartbyElementId[rElement.GetId()],
                                   rModelPart.GetProcessInfo());
            }
            mIsInitialized = true;
        }
    }

private:
    std::map<std::size_t, Vector> mExternalForcesAtStartbyElementId;
    std::map<std::size_t, Vector> mInternalForcesAtStartbyElementId;

    bool mIsInitialized = false;
};
} // namespace Kratos