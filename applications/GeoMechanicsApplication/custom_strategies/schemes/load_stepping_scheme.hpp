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

#pragma once

#include "geomechanics_static_scheme.hpp"

namespace Kratos
{
template <class TSparseSpace, class TDenseSpace>
class LoadSteppingScheme : public GeoMechanicsStaticScheme<TSparseSpace, TDenseSpace>
{
public:
    void CalculateRHSContribution(Element& rCurrentElement,
                                  GeoMechanicsTimeIntegrationScheme<TSparseSpace, TDenseSpace>::LocalSystemVectorType& RHS_Contribution,
                                  Element::EquationIdVectorType& EquationId,
                                  const ProcessInfo&             CurrentProcessInfo) override
    {
        const auto fraction_of_unbalance = (CurrentProcessInfo[TIME] - CurrentProcessInfo[START_TIME]) /
                                     (CurrentProcessInfo[END_TIME] - CurrentProcessInfo[START_TIME]);
        rCurrentElement.Calculate(INTERNAL_FORCES_VECTOR, RHS_Contribution, CurrentProcessInfo);
        RHS_Contribution -= mInternalForcesAtStartByElementId.at(rCurrentElement.GetId());
        RHS_Contribution +=
            fraction_of_unbalance * (mInternalForcesAtStartByElementId.at(rCurrentElement.GetId()) +
                                     mExternalForcesAtStartByElementId.at(rCurrentElement.GetId()));
        rCurrentElement.EquationIdVector(EquationId, CurrentProcessInfo);
    }

    using GeoMechanicsStaticScheme<TSparseSpace, TDenseSpace>::CalculateRHSContribution;

    void CalculateSystemContributions(
        Element& rCurrentElement,
        GeoMechanicsTimeIntegrationScheme<TSparseSpace, TDenseSpace>::LocalSystemMatrixType& LHS_Contribution,
        GeoMechanicsTimeIntegrationScheme<TSparseSpace, TDenseSpace>::LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        const ProcessInfo&             CurrentProcessInfo) override
    {
        CalculateRHSContribution(rCurrentElement, RHS_Contribution, EquationId, CurrentProcessInfo);
        this->CalculateLHSContribution(rCurrentElement, LHS_Contribution, EquationId, CurrentProcessInfo);
    }

    void CalculateSystemContributions(
        Condition& rCurrentCondition,
        GeoMechanicsTimeIntegrationScheme<TSparseSpace, TDenseSpace>::LocalSystemMatrixType& LHS_Contribution,
        GeoMechanicsTimeIntegrationScheme<TSparseSpace, TDenseSpace>::LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        const ProcessInfo&             CurrentProcessInfo) override
    {
        CalculateRHSContribution(rCurrentCondition, RHS_Contribution, EquationId, CurrentProcessInfo);
        this->CalculateLHSContribution(rCurrentCondition, LHS_Contribution, EquationId, CurrentProcessInfo);
    }

    using GeoMechanicsStaticScheme<TSparseSpace, TDenseSpace>::CalculateSystemContributions;

    void CalculateRHSContribution(Condition& rCurrentCondition,
                                  GeoMechanicsTimeIntegrationScheme<TSparseSpace, TDenseSpace>::LocalSystemVectorType& RHS_Contribution,
                                  Element::EquationIdVectorType& EquationId,
                                  const ProcessInfo&             CurrentProcessInfo) override
    {
        GeoMechanicsStaticScheme<TSparseSpace, TDenseSpace>::CalculateRHSContribution(
            rCurrentCondition, RHS_Contribution, EquationId, CurrentProcessInfo);

        const auto fraction_of_unbalance = (CurrentProcessInfo[TIME] - CurrentProcessInfo[START_TIME]) /
                                     (CurrentProcessInfo[END_TIME] - CurrentProcessInfo[START_TIME]);
        RHS_Contribution *= fraction_of_unbalance;
    }

    void InitializeSolutionStep(
        ModelPart&                                                                       rModelPart,
        GeoMechanicsTimeIntegrationScheme<TSparseSpace, TDenseSpace>::TSystemMatrixType& rA,
        GeoMechanicsTimeIntegrationScheme<TSparseSpace, TDenseSpace>::TSystemVectorType& rDx,
        GeoMechanicsTimeIntegrationScheme<TSparseSpace, TDenseSpace>::TSystemVectorType& rB) override
    {
        GeoMechanicsStaticScheme<TSparseSpace, TDenseSpace>::InitializeSolutionStep(rModelPart, rA, rDx, rB);

        if (!mIsInitialized) {
            for (auto& rElement : rModelPart.Elements()) {
                mInternalForcesAtStartByElementId.insert({rElement.GetId(), Vector{}});
                rElement.Calculate(INTERNAL_FORCES_VECTOR,
                                   mInternalForcesAtStartByElementId[rElement.GetId()],
                                   rModelPart.GetProcessInfo());

                mExternalForcesAtStartByElementId.insert({rElement.GetId(), Vector{}});
                rElement.Calculate(EXTERNAL_FORCES_VECTOR,
                                   mExternalForcesAtStartByElementId[rElement.GetId()],
                                   rModelPart.GetProcessInfo());
            }
            mIsInitialized = true;
        }
    }

private:
    std::map<std::size_t, Vector> mInternalForcesAtStartByElementId;
    std::map<std::size_t, Vector> mExternalForcesAtStartByElementId;
    bool                          mIsInitialized = false;
};
} // namespace Kratos