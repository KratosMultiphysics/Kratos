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
#include "custom_strategies/schemes/backward_euler_quasistatic_Pw_scheme.hpp"
#include "custom_strategies/schemes/load_stepping_scheme.hpp"
#include "custom_strategies/schemes/newmark_dynamic_U_Pw_scheme.hpp"
#include "custom_strategies/schemes/newmark_quasistatic_damped_U_Pw_scheme.hpp"
#include "custom_strategies/schemes/newmark_quasistatic_Pw_scheme.hpp"
#include "custom_strategies/schemes/newmark_quasistatic_U_Pw_scheme.hpp"
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

        const auto scheme_type   = rSolverSettings["scheme_type"].GetString();
        const auto solution_type = rSolverSettings["solution_type"].GetString();
        
        // Determine the solver type (Pw or U_Pw)
        std::string solver_type = "U_Pw";
        if (rSolverSettings.Has("solver_type")) {
            solver_type = rSolverSettings["solver_type"].GetString();
        }

        // Backward Euler schemes
        if (scheme_type == "Backward_Euler") {
            if (solution_type == "Quasi-Static") {
                if (solver_type == "Pw") {
                    return std::make_shared<BackwardEulerQuasistaticPwScheme<TSparseSpace, TDenseSpace>>();
                }
                // Default to U-Pw scheme
                return std::make_shared<BackwardEulerQuasistaticUPwScheme<TSparseSpace, TDenseSpace>>();
            }
        }

        // Newmark schemes
        if (scheme_type == "Newmark") {
            KRATOS_ERROR_IF_NOT(rSolverSettings.Has("newmark_beta"))
                << "'newmark_beta' is not defined, aborting";
            KRATOS_ERROR_IF_NOT(rSolverSettings.Has("newmark_gamma"))
                << "'newmark_gamma' is not defined, aborting";
            KRATOS_ERROR_IF_NOT(rSolverSettings.Has("newmark_theta"))
                << "'newmark_theta' is not defined, aborting";

            const auto beta  = rSolverSettings["newmark_beta"].GetDouble();
            const auto gamma = rSolverSettings["newmark_gamma"].GetDouble();
            const auto theta = rSolverSettings["newmark_theta"].GetDouble();

            if (solution_type == "dynamic") {
                return std::make_shared<NewmarkDynamicUPwScheme<TSparseSpace, TDenseSpace>>(beta, gamma, theta);
            }
            
            if (solution_type == "Quasi-Static") {
                if (solver_type == "Pw") {
                    return std::make_shared<NewmarkQuasistaticPwScheme<TSparseSpace, TDenseSpace>>(theta);
                }

                const auto rayleigh_m = rSolverSettings.Has("rayleigh_m") ? rSolverSettings["rayleigh_m"].GetDouble() : 0.0;
                const auto rayleigh_k = rSolverSettings.Has("rayleigh_k") ? rSolverSettings["rayleigh_k"].GetDouble() : 0.0;

                if (rayleigh_m < 1.0e-20 && rayleigh_k < 1.0e-20) {
                    return std::make_shared<NewmarkQuasistaticUPwScheme<TSparseSpace, TDenseSpace>>(
                        beta, gamma, theta);
                } else {
                    return std::make_shared<NewmarkQuasistaticDampedUPwScheme<TSparseSpace, TDenseSpace>>(
                        beta, gamma, theta);
                }
            }
        }

        // Static schemes
        if (solution_type == "static") {
            if (scheme_type == "load_stepping") {
                return std::make_shared<LoadSteppingScheme<TSparseSpace, TDenseSpace>>();
            }
            return std::make_shared<GeoMechanicsStaticScheme<TSparseSpace, TDenseSpace>>();
        }

        KRATOS_ERROR << "Specified combination of solution_type (" << solution_type
                     << ") and scheme_type (" << scheme_type << ") is not supported, aborting";
    }
};

} // namespace Kratos
