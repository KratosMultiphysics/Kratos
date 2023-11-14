// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Vahid Galavi
//

#pragma once

// Project includes
#include "geomechanics_time_integration_scheme.hpp"
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/schemes/scheme.h"
#include "utilities/parallel_utilities.h"

// Application includes
#include "generalized_backward_euler_scheme.hpp"
#include "geo_mechanics_application_variables.h"

namespace Kratos {

template <class TSparseSpace, class TDenseSpace>
class BackwardEulerQuasistaticPwScheme
    : public GeneralizedBackwardEulerScheme<TSparseSpace, TDenseSpace> {
public:
    KRATOS_CLASS_POINTER_DEFINITION(BackwardEulerQuasistaticPwScheme);

    BackwardEulerQuasistaticPwScheme()
        : GeneralizedBackwardEulerScheme<TSparseSpace, TDenseSpace>(
              WATER_PRESSURE, DT_WATER_PRESSURE, DT_PRESSURE_COEFFICIENT)
    {
    }
}; // Class BackwardEulerQuasistaticPwScheme

} // namespace Kratos