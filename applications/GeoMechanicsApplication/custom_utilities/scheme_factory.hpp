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

#include "solving_strategies/schemes/scheme.h"
#include "custom_strategies/schemes/backward_euler_quasistatic_U_Pw_scheme.hpp"
#include <memory>

namespace Kratos
{

template<class TSparseSpace, class TDenseSpace>
class SchemeFactory
{
public:
    using SchemeType = Scheme<TSparseSpace, TDenseSpace>;

    static std::shared_ptr<SchemeType> Create(const Parameters& rSolverSettings)
    {
        KRATOS_ERROR_IF_NOT(rSolverSettings.Has("scheme_type")) << "scheme_type is not defined, aborting";
        KRATOS_ERROR_IF_NOT(rSolverSettings.Has("solution_type")) << "solution_type is not defined, aborting";

        if (rSolverSettings["scheme_type"].GetString() == "Backward_Euler" &&
            rSolverSettings["solution_type"].GetString() == "Quasi-Static")
        {
            return std::make_shared<BackwardEulerQuasistaticUPwScheme<TSparseSpace, TDenseSpace>>();
        }

        KRATOS_ERROR << "Specified solution_type/scheme_type is not supported, aborting";
    }
};

}
