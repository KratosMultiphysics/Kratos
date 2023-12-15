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
#include "solving_strategies/convergencecriterias/residual_criteria.h"
#include "parameters_utilities.h"

namespace Kratos
{

template <class TSparseSpace, class TDenseSpace>
class ConvergenceCriteriaFactory
{
public:
    using ConvergenceCriteriaType = ConvergenceCriteria<TSparseSpace, TDenseSpace>;

    static std::shared_ptr<ConvergenceCriteriaType> Create(const Parameters& rSolverSettings)
    {
        KRATOS_ERROR_IF_NOT(rSolverSettings.Has("convergence_criterion"))
            << "No convergence_criterion is defined, aborting.";

        if (rSolverSettings["convergence_criterion"].GetString() ==
            "displacement_criterion")
        {
            const std::vector<std::string> entries_to_copy = {
                "displacement_absolute_tolerance",
                "displacement_relative_tolerance"};
            const Parameters convergence_inputs =
                ParametersUtilities::CopyOptionalParameters(rSolverSettings, entries_to_copy);
            return std::make_shared<DisplacementCriteria<TSparseSpace, TDenseSpace>>(convergence_inputs);
        }
        if (rSolverSettings["convergence_criterion"].GetString() ==
            "residual_criterion")
        {
            const std::vector<std::string> entries_to_copy = {
                "residual_absolute_tolerance", "residual_relative_tolerance"};
            const auto convergence_inputs = ParametersUtilities::CopyOptionalParameters(
                rSolverSettings, entries_to_copy);
            return std::make_shared<ResidualCriteria<TSparseSpace, TDenseSpace>>(convergence_inputs);
        }

        KRATOS_ERROR << "The convergence_criterion ("
                     << rSolverSettings["convergence_criterion"].GetString() << ") is unknown, "
                     << "supported criteria are: 'displacement_criterion', "
                        "'residual_criterion'."
                     << std::endl;
    }
};

} // namespace Kratos
