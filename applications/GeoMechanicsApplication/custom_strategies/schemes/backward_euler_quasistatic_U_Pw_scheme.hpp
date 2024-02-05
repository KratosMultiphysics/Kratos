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
#include "backward_euler_scheme.hpp"
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/schemes/scheme.h"
#include "utilities/parallel_utilities.h"

// Application includes
#include "geo_mechanics_application_variables.h"

namespace Kratos
{

template <class TSparseSpace, class TDenseSpace>
class BackwardEulerQuasistaticUPwScheme : public BackwardEulerScheme<TSparseSpace, TDenseSpace>
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(BackwardEulerQuasistaticUPwScheme);

    BackwardEulerQuasistaticUPwScheme()
        : BackwardEulerScheme<TSparseSpace, TDenseSpace>(
              {FirstOrderScalarVariable(WATER_PRESSURE, DT_WATER_PRESSURE, DT_PRESSURE_COEFFICIENT)},
              {SecondOrderVectorVariable(DISPLACEMENT), SecondOrderVectorVariable(ROTATION)})
    {
    }

protected:
    std::string Info() const override { return "BackwardEulerQuasistaticUPwScheme"; }
}; // Class BackwardEulerQuasistaticUPwScheme

} // namespace Kratos