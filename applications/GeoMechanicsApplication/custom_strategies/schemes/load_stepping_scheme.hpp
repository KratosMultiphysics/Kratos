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
    KRATOS_CLASS_POINTER_DEFINITION(LoadSteppingScheme);

    using TSystemVectorType =
        typename GeoMechanicsTimeIntegrationScheme<TSparseSpace, TDenseSpace>::TSystemVectorType;
    using TSystemMatrixType =
        typename GeoMechanicsTimeIntegrationScheme<TSparseSpace, TDenseSpace>::TSystemMatrixType;
    using TLocalSystemMatrixType =
        typename GeoMechanicsTimeIntegrationScheme<TSparseSpace, TDenseSpace>::LocalSystemMatrixType;
    using TLocalSystemVectorType =
        typename GeoMechanicsTimeIntegrationScheme<TSparseSpace, TDenseSpace>::LocalSystemVectorType;

    void CalculateSystemContributions(Element&                       rCurrentElement,
                                      TLocalSystemMatrixType&        rLHS_Contribution,
                                      TLocalSystemVectorType&        rRHS_Contribution,
                                      Element::EquationIdVectorType& rEquationId,
                                      const ProcessInfo&             rCurrentProcessInfo) override
    {
        this->CalculateLHSContribution(rCurrentElement, rLHS_Contribution, rEquationId, rCurrentProcessInfo);
        CalculateRHSContribution(rCurrentElement, rRHS_Contribution, rEquationId, rCurrentProcessInfo);
    }

    void CalculateRHSContribution(Element&                       rCurrentElement,
                                  TLocalSystemVectorType&        rRHS_Contribution,
                                  Element::EquationIdVectorType& rEquationId,
                                  const ProcessInfo&             rCurrentProcessInfo) override
    {
        KRATOS_ERROR_IF(mExternalForcesAtStartByElementId.empty() ||
                        mInternalForcesAtStartByElementId.empty())
            << "The load stepping scheme is not initialized properly. Make sure all "
               "initializations are done before calculating system contributions.";

        Vector internal_forces;
        rCurrentElement.Calculate(INTERNAL_FORCES_VECTOR, internal_forces, rCurrentProcessInfo);

        const auto load_fraction = CalculateLoadFraction(rCurrentProcessInfo);
        rRHS_Contribution =
            mInternalForcesAtStartByElementId.at(rCurrentElement.GetId()) +
            load_fraction * (mExternalForcesAtStartByElementId.at(rCurrentElement.GetId()) -
                             mInternalForcesAtStartByElementId.at(rCurrentElement.GetId())) -
            internal_forces;

        rCurrentElement.EquationIdVector(rEquationId, rCurrentProcessInfo);
    }

    void CalculateSystemContributions(Condition&                     rCurrentCondition,
                                      TLocalSystemMatrixType&        rLHS_Contribution,
                                      TLocalSystemVectorType&        rRHS_Contribution,
                                      Element::EquationIdVectorType& rEquationId,
                                      const ProcessInfo&             rCurrentProcessInfo) override
    {
        GeoMechanicsTimeIntegrationScheme<TSparseSpace, TDenseSpace>::CalculateSystemContributions(
            rCurrentCondition, rLHS_Contribution, rRHS_Contribution, rEquationId, rCurrentProcessInfo);
        rRHS_Contribution *= CalculateLoadFraction(rCurrentProcessInfo);
    }

    void CalculateRHSContribution(Condition&                     rCurrentCondition,
                                  TLocalSystemVectorType&        rRHS_Contribution,
                                  Element::EquationIdVectorType& rEquationId,
                                  const ProcessInfo&             rCurrentProcessInfo) override
    {
        GeoMechanicsTimeIntegrationScheme<TSparseSpace, TDenseSpace>::CalculateRHSContribution(
            rCurrentCondition, rRHS_Contribution, rEquationId, rCurrentProcessInfo);

        rRHS_Contribution *= CalculateLoadFraction(rCurrentProcessInfo);
    }

    void InitializeSolutionStep(ModelPart& rModelPart, TSystemMatrixType& rA, TSystemVectorType& rDx, TSystemVectorType& rB) override
    {
        GeoMechanicsStaticScheme<TSparseSpace, TDenseSpace>::InitializeSolutionStep(rModelPart, rA, rDx, rB);

        // Saving the internal and external forces should be done only once. However, the Initialize
        // function (which would normally be the place for initializations like this), is called too
        // early to be able to calculate external and internal forces (which is why it's done in
        // InitializeSolutionStep).
        if (!mIsInitialized) {
            mInternalForcesAtStartByElementId = CalculateElementForces(rModelPart, INTERNAL_FORCES_VECTOR);
            mExternalForcesAtStartByElementId = CalculateElementForces(rModelPart, EXTERNAL_FORCES_VECTOR);

            mIsInitialized = true;
        }
    }

    void FinalizeSolutionStep(ModelPart& rModelPart, TSystemMatrixType& rA, TSystemVectorType& rDx, TSystemVectorType& rb) override
    {
        GeoMechanicsStaticScheme<TSparseSpace, TDenseSpace>::FinalizeSolutionStep(rModelPart, rA, rDx, rb);
        KRATOS_INFO("Load stepping")
            << "Fraction of unbalance: " << CalculateLoadFraction(rModelPart.GetProcessInfo()) << "\n";
    }

private:
    std::map<std::size_t, TLocalSystemVectorType> mInternalForcesAtStartByElementId;
    std::map<std::size_t, TLocalSystemVectorType> mExternalForcesAtStartByElementId;
    bool                                          mIsInitialized = false;

    static double CalculateLoadFraction(const ProcessInfo& rProcessInfo)
    {
        return (rProcessInfo[TIME] - rProcessInfo[START_TIME]) /
               (rProcessInfo[END_TIME] - rProcessInfo[START_TIME]);
    }

    std::map<std::size_t, TLocalSystemVectorType> CalculateElementForces(ModelPart& rModelPart,
                                                                         Variable<Vector>& rForces)
    {
        std::map<std::size_t, TLocalSystemVectorType> result;
        std::ranges::transform(rModelPart.Elements(), std::inserter(result, result.end()),
                               [&rModelPart, &rForces](auto& rElement) {
            Vector forces;
            rElement.Calculate(rForces, forces, rModelPart.GetProcessInfo());
            return std::make_pair(rElement.GetId(), forces);
        });

        return result;
    }
};
} // namespace Kratos