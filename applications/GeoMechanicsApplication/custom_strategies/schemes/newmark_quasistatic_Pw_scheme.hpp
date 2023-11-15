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
#include "includes/model_part.h"
#include "solving_strategies/schemes/scheme.h"

// Application includes
#include "custom_strategies/schemes/newmark_quasistatic_U_Pw_scheme.hpp"
#include "geo_mechanics_application_variables.h"

namespace Kratos {

template <class TSparseSpace, class TDenseSpace>
class NewmarkQuasistaticPwScheme
    : public GeneralizedNewmarkScheme<TSparseSpace, TDenseSpace> {
public:
    KRATOS_CLASS_POINTER_DEFINITION(NewmarkQuasistaticPwScheme);

    using BaseType = Scheme<TSparseSpace, TDenseSpace>;
    using TSystemMatrixType = typename BaseType::TSystemMatrixType;
    using TSystemVectorType = typename BaseType::TSystemVectorType;

    explicit NewmarkQuasistaticPwScheme(double theta)
        : GeneralizedNewmarkScheme<TSparseSpace, TDenseSpace>(
              theta, WATER_PRESSURE, DT_WATER_PRESSURE, DT_PRESSURE_COEFFICIENT)
    {
    }

}; // Class NewmarkQuasistaticPwScheme

} // namespace Kratos