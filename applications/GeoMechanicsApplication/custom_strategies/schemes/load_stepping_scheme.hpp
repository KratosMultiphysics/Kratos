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
                                  const ProcessInfo& CurrentProcessInfo) override
    {
        rCurrentElement.Calculate(INTERNAL_FORCES_VECTOR, RHS_Contribution, CurrentProcessInfo);
        RHS_Contribution -= mInternalForcesAtStartByElementId.at(rCurrentElement.GetId());
    }

    void InitializeSolutionStep(ModelPart& rModelPart,
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
            }
            mIsInitialized = true;
        }
    }


private:
    std::map<std::size_t, Vector> mInternalForcesAtStartByElementId;
    bool mIsInitialized = false;
};
}