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


#include "solving_strategies/convergencecriterias/mixed_generic_criteria.h"

namespace Kratos
{

template<class TSparseSpace, class TDenseSpace>
class ConvergenceCriteriaFactory
{
public:
    using MixedGenericCriteriaType = MixedGenericCriteria<TSparseSpace, TDenseSpace>;
    using ConvergenceVariableListType = typename MixedGenericCriteriaType::ConvergenceVariableListType;
    using ConvergenceCriteriaType = ConvergenceCriteria<TSparseSpace, TDenseSpace>;

    static std::shared_ptr<ConvergenceCriteriaType> Create(const Parameters& rSolverSettings)
    {
        KRATOS_ERROR_IF_NOT(rSolverSettings.Has("convergence_criterion"))
        << "No convergence_criterion is defined, aborting.";

        if (rSolverSettings["convergence_criterion"].GetString() == "displacement_criterion")
        {
            double rel_tol = 0;
            if (rSolverSettings.Has("displacement_relative_tolerance"))
            {
                rel_tol = rSolverSettings["displacement_relative_tolerance"].GetDouble();
            }

            double abs_tol = 0;
            if (rSolverSettings.Has("displacement_absolute_tolerance"))
            {
                abs_tol = rSolverSettings["displacement_absolute_tolerance"].GetDouble();
            }
            VariableData* p_displacement = &DISPLACEMENT;
            ConvergenceVariableListType convergence_settings{std::make_tuple(p_displacement, rel_tol, abs_tol)};
            return std::make_shared<MixedGenericCriteriaType>(MixedGenericCriteriaType(convergence_settings));
        }

        KRATOS_ERROR << "The convergence_criterion (" << rSolverSettings["convergence_criterion"].GetString() << ") is unknown, supported criteria are: 'displacement_criterion'";
    }
};

}
