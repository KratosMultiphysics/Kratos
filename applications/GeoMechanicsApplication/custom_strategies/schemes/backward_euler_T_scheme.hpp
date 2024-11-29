// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Mohamed Nabi
//                   John van Esch
//                   Richard Faasse
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
class BackwardEulerTScheme : public BackwardEulerScheme<TSparseSpace, TDenseSpace>
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(BackwardEulerTScheme);

    BackwardEulerTScheme()
        : BackwardEulerScheme<TSparseSpace, TDenseSpace>(
              {FirstOrderScalarVariable(TEMPERATURE, DT_TEMPERATURE, DT_TEMPERATURE_COEFFICIENT)}, {})
    {
    }
}; // Class BackwardEulerTScheme
} // namespace Kratos
