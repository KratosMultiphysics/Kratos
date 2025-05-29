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

#include "custom_strategies/schemes/backward_euler_quasistatic_U_Pw_scheme.hpp"
#include "custom_strategies/schemes/newmark_dynamic_U_Pw_scheme.hpp"
#include "solving_strategies/schemes/scheme.h"
#include <memory>

namespace Kratos
{

template <class TSparseSpace, class TDenseSpace>
class SchemeFactory
{
public:
    using SchemeType = Scheme<TSparseSpace, TDenseSpace>;

    static std::shared_ptr<SchemeType> Create(const Parameters& rSolverSettings)
    {
        KRATOS_ERROR_IF_NOT(rSolverSettings.Has("scheme_type"))
            << "scheme_type is not defined, aborting";
        KRATOS_ERROR_IF_NOT(rSolverSettings.Has("solution_type"))
            << "solution_type is not defined, aborting";

        if (rSolverSettings["scheme_type"].GetString() == "Backward_Euler" &&
            rSolverSettings["solution_type"].GetString() == "Quasi-Static") {
            return std::make_shared<BackwardEulerQuasistaticUPwScheme<TSparseSpace, TDenseSpace>>();
        }

        if (rSolverSettings["scheme_type"].GetString() == "Newmark" &&
            rSolverSettings["solution_type"].GetString() == "dynamic") {
            KRATOS_ERROR_IF_NOT(rSolverSettings.Has("newmark_beta"))
                << "'newmark_beta' is not defined, aborting";
            KRATOS_ERROR_IF_NOT(rSolverSettings.Has("newmark_gamma"))
                << "'newmark_gamma' is not defined, aborting";
            KRATOS_ERROR_IF_NOT(rSolverSettings.Has("newmark_theta"))
                << "'newmark_theta' is not defined, aborting";

            return std::make_shared<NewmarkDynamicUPwScheme<TSparseSpace, TDenseSpace>>(
                rSolverSettings["newmark_beta"].GetDouble(), rSolverSettings["newmark_gamma"].GetDouble(),
                rSolverSettings["newmark_theta"].GetDouble());
        }

        KRATOS_ERROR << "Specified combination of solution_type ("
                     << rSolverSettings["solution_type"].GetString() << ") and scheme_type ("
                     << rSolverSettings["scheme_type"].GetString() << ") is not supported, aborting";
    }
};

} // namespace Kratos
