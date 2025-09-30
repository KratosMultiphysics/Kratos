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
    using TSystemVectorType = GeoMechanicsTimeIntegrationScheme<TSparseSpace, TDenseSpace>::TSystemVectorType;
    using TSystemMatrixType = GeoMechanicsTimeIntegrationScheme<TSparseSpace, TDenseSpace>::TSystemMatrixType;
    using TLocalSystemMatrixType = GeoMechanicsTimeIntegrationScheme<TSparseSpace, TDenseSpace>::LocalSystemMatrixType;
    using TLocalSystemVectorType = GeoMechanicsTimeIntegrationScheme<TSparseSpace, TDenseSpace>::LocalSystemVectorType;

    void CalculateSystemContributions(Element&                       rCurrentElement,
                                      TLocalSystemMatrixType&        LHS_Contribution,
                                      TLocalSystemVectorType&        RHS_Contribution,
                                      Element::EquationIdVectorType& EquationId,
                                      const ProcessInfo&             CurrentProcessInfo) override
    {
        CalculateRHSContribution(rCurrentElement, RHS_Contribution, EquationId, CurrentProcessInfo);
        this->CalculateLHSContribution(rCurrentElement, LHS_Contribution, EquationId, CurrentProcessInfo);
    }

    void CalculateRHSContribution(Element&                       rCurrentElement,
                                  TLocalSystemVectorType&        RHS_Contribution,
                                  Element::EquationIdVectorType& EquationId,
                                  const ProcessInfo&             CurrentProcessInfo) override
    {
        Vector internal_forces;
        rCurrentElement.Calculate(INTERNAL_FORCES_VECTOR, internal_forces, CurrentProcessInfo);

        const auto fraction_of_unbalance = CalculateFractionOfUnbalance(CurrentProcessInfo);
        RHS_Contribution =
            -internal_forces +
            (1 - fraction_of_unbalance) * mInternalForcesAtStartByElementId.at(rCurrentElement.GetId()) +
            fraction_of_unbalance * mExternalForcesAtStartByElementId.at(rCurrentElement.GetId());

        rCurrentElement.EquationIdVector(EquationId, CurrentProcessInfo);
    }

    void CalculateSystemContributions(Condition&                     rCurrentCondition,
                                      TLocalSystemMatrixType&        LHS_Contribution,
                                      TLocalSystemVectorType&        RHS_Contribution,
                                      Element::EquationIdVectorType& EquationId,
                                      const ProcessInfo&             CurrentProcessInfo) override
    {
        CalculateRHSContribution(rCurrentCondition, RHS_Contribution, EquationId, CurrentProcessInfo);
        this->CalculateLHSContribution(rCurrentCondition, LHS_Contribution, EquationId, CurrentProcessInfo);
    }

    void CalculateRHSContribution(Condition&                     rCurrentCondition,
                                  TLocalSystemVectorType&        RHS_Contribution,
                                  Element::EquationIdVectorType& EquationId,
                                  const ProcessInfo&             CurrentProcessInfo) override
    {
        GeoMechanicsStaticScheme<TSparseSpace, TDenseSpace>::CalculateRHSContribution(
            rCurrentCondition, RHS_Contribution, EquationId, CurrentProcessInfo);

        RHS_Contribution *= CalculateFractionOfUnbalance(CurrentProcessInfo);
    }

    void InitializeSolutionStep(ModelPart& rModelPart, TSystemMatrixType& rA, TSystemVectorType& rDx, TSystemVectorType& rB) override
    {
        GeoMechanicsStaticScheme<TSparseSpace, TDenseSpace>::InitializeSolutionStep(rModelPart, rA, rDx, rB);

        // Saving the internal and external forces should be done only once. However, the Initialize
        // function (which would normally be the place for initializations like this), is called too
        // early to be able to calculate external and internal forces (which is why it's done in
        // InitializeSolutionStep).
        if (!mIsInitialized) {
            SaveInternalForces(rModelPart);
            SaveExternalForces(rModelPart);

            mIsInitialized = true;
        }
    }

private:
    std::map<std::size_t, TLocalSystemVectorType> mInternalForcesAtStartByElementId;
    std::map<std::size_t, TLocalSystemVectorType> mExternalForcesAtStartByElementId;
    bool                                          mIsInitialized = false;

    static double CalculateFractionOfUnbalance(const ProcessInfo& rProcessInfo)
    {
        return (rProcessInfo[TIME] - rProcessInfo[START_TIME]) /
               (rProcessInfo[END_TIME] - rProcessInfo[START_TIME]);
    }

    void SaveInternalForces(ModelPart& rModelPart)
    {
        std::ranges::transform(
            rModelPart.Elements(),
            std::inserter(mInternalForcesAtStartByElementId, mInternalForcesAtStartByElementId.end()),
            [&rModelPart](auto& rElement) {
            Vector forces;
            rElement.Calculate(INTERNAL_FORCES_VECTOR, forces, rModelPart.GetProcessInfo());
            return std::make_pair(rElement.GetId(), forces);
        });
    }

    void SaveExternalForces(ModelPart& rModelPart)
    {
        std::ranges::transform(
            rModelPart.Elements(),
            std::inserter(mExternalForcesAtStartByElementId, mExternalForcesAtStartByElementId.end()),
            [&rModelPart](auto& rElement) {
            Vector forces;
            rElement.Calculate(EXTERNAL_FORCES_VECTOR, forces, rModelPart.GetProcessInfo());
            return std::make_pair(rElement.GetId(), forces);
        });
    }
};
} // namespace Kratos