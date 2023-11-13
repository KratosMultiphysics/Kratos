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

namespace Kratos {

template <class TSparseSpace, class TDenseSpace>
class GeneralizedNewmarkScheme
    : public GeoMechanicsTimeIntegrationScheme<TSparseSpace, TDenseSpace> {
public:
    explicit GeneralizedNewmarkScheme(double theta)
        : GeoMechanicsTimeIntegrationScheme<TSparseSpace, TDenseSpace>(), mTheta(theta)
    {
    }

protected:
    void UpdateScalarTimeDerivative(Node& rNode,
                                    const Variable<double>& variable,
                                    const Variable<double>& dt_variable) const
    {
        const double delta_variable = rNode.FastGetSolutionStepValue(variable) -
                                      rNode.FastGetSolutionStepValue(variable, 1);
        const auto& previous_dt_variable =
            rNode.FastGetSolutionStepValue(dt_variable, 1);

        rNode.FastGetSolutionStepValue(dt_variable) =
            (delta_variable - (1.0 - mTheta) * this->GetDeltaTime() * previous_dt_variable) /
            (mTheta * this->GetDeltaTime());
    }

    double mTheta = 0.0;
};

} // namespace Kratos
