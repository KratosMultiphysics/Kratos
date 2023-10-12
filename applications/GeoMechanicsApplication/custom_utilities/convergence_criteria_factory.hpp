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

#include "solving_strategies/convergencecriterias/displacement_criteria.h"

namespace Kratos
{

template<class TSparseSpace, class TDenseSpace>
class ConvergenceCriteriaFactory
{
public:
    using ConvergenceCriteriaType = ConvergenceCriteria<TSparseSpace, TDenseSpace>;

    static std::shared_ptr<ConvergenceCriteriaType> Create(const Parameters& rSolverSettings)
    {
        KRATOS_ERROR_IF_NOT(rSolverSettings.Has("convergence_criterion"))
        << "No convergence_criterion is defined, aborting.";

        if (rSolverSettings["convergence_criterion"].GetString() == "displacement_criterion")
        {
            Parameters convergenceInputs;

            std::vector<std::string> entries_to_copy = {"displacement_absolute_tolerance", "displacement_relative_tolerance"};
            for (const std::string& entry : entries_to_copy)
            {
                if (rSolverSettings.Has(entry))
                {
                    convergenceInputs.AddValue(entry, rSolverSettings[entry]);
                }
            }
            return std::make_shared<DisplacementCriteria<TSparseSpace, TDenseSpace>>(convergenceInputs);
        }

        KRATOS_ERROR << "The convergence_criterion (" << rSolverSettings["convergence_criterion"].GetString() << ") is unknown, supported criteria are: 'displacement_criterion'";
    }
};

}
