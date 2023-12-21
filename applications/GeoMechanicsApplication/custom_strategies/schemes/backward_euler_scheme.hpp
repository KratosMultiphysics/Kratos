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

#include "geomechanics_time_integration_scheme.hpp"

namespace Kratos
{

template <class TSparseSpace, class TDenseSpace>
class BackwardEulerScheme
    : public GeoMechanicsTimeIntegrationScheme<TSparseSpace, TDenseSpace>
{
public:
    BackwardEulerScheme(const Variable<double>& rVariable,
                        const Variable<double>& rDeltaTimeVariable,
                        const Variable<double>& rDeltaTimeVariableCoefficient)
        : GeoMechanicsTimeIntegrationScheme<TSparseSpace, TDenseSpace>(
              rVariable, rDeltaTimeVariable, rDeltaTimeVariableCoefficient)
    {
    }

protected:

    inline void SetTimeFactors(ModelPart& rModelPart) override
    {
        KRATOS_TRY

        const auto delta_time = rModelPart.GetProcessInfo()[DELTA_TIME];
        this->SetDeltaTime(delta_time);
        rModelPart.GetProcessInfo()[this->mDeltaTimeVariableCoefficient] = 1.0 / delta_time;

        KRATOS_CATCH("")
    }

    void UpdateScalarTimeDerivative(Node& rNode,
                                    const Variable<double>& variable,
                                    const Variable<double>& dt_variable) const override
    {
        const auto delta_variable = rNode.FastGetSolutionStepValue(variable, 0) -
                                    rNode.FastGetSolutionStepValue(variable, 1);
        rNode.FastGetSolutionStepValue(dt_variable, 0) =
            delta_variable / this->GetDeltaTime();
    }
};

} // namespace Kratos
