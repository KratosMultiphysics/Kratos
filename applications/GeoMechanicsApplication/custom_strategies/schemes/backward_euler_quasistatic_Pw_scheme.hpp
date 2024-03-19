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
#include "includes/define.h"

// Application includes
#include "backward_euler_scheme.hpp"
#include "geo_mechanics_application_variables.h"

namespace Kratos
{

template <class TSparseSpace, class TDenseSpace>
class BackwardEulerQuasistaticPwScheme : public BackwardEulerScheme<TSparseSpace, TDenseSpace>
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(BackwardEulerQuasistaticPwScheme);

    BackwardEulerQuasistaticPwScheme()
        : BackwardEulerScheme<TSparseSpace, TDenseSpace>(
              {FirstOrderScalarVariable(WATER_PRESSURE, DT_WATER_PRESSURE, DT_PRESSURE_COEFFICIENT)}, {})
    {
    }
}; // Class BackwardEulerQuasistaticPwScheme

} // namespace Kratos