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

#include "custom_utilities/string_utilities.h"
#include "includes/variables.h" // for WATER_PRESSURE
#include "parameters_utilities.h"
#include "solving_strategies/convergencecriterias/and_criteria.h"
#include "solving_strategies/convergencecriterias/displacement_criteria.h"
#include "solving_strategies/convergencecriterias/mixed_generic_criteria.h"
#include "solving_strategies/convergencecriterias/or_criteria.h"
#include "solving_strategies/convergencecriterias/residual_criteria.h"

#include <string>

namespace Kratos
{

template <class TSparseSpace, class TDenseSpace>
class ConvergenceCriteriaFactory
{
public:
    using ConvergenceCriterionType      = ConvergenceCriteria<TSparseSpace, TDenseSpace>;
    using ConvergenceCriterionSharedPtr = std::shared_ptr<ConvergenceCriterionType>;
    using MixedGenericCriterionType     = MixedGenericCriteria<TSparseSpace, TDenseSpace>;

    static ConvergenceCriterionSharedPtr Create(const Parameters& rSolverSettings)
    {
        using namespace std::string_literals;

        KRATOS_ERROR_IF_NOT(rSolverSettings.Has("convergence_criterion"s))
            << "No convergence_criterion is defined, aborting.";

        const auto convergence_criterion_type =
            GeoStringUtilities::ToLower(rSolverSettings["convergence_criterion"s].GetString());

        // The following switch resembles the one in method `_ConstructConvergenceCriterion` (see
        // `geomechanics_U_Pw_solver.py`)
        if (convergence_criterion_type == "displacement_criterion") {
            return CreateDisplacementCriterion(rSolverSettings);
        }

        if (convergence_criterion_type == "residual_criterion") {
            return CreateResidualCriterion(rSolverSettings);
        }

        if (convergence_criterion_type == "and_criterion") {
            return std::make_shared<And_Criteria<TSparseSpace, TDenseSpace>>(
                CreateResidualCriterion(rSolverSettings), CreateDisplacementCriterion(rSolverSettings));
        }

        if (convergence_criterion_type == "or_criterion") {
            return std::make_shared<Or_Criteria<TSparseSpace, TDenseSpace>>(
                CreateResidualCriterion(rSolverSettings), CreateDisplacementCriterion(rSolverSettings));
        }

        if (convergence_criterion_type == "water_pressure_criterion") {
            return CreateWaterPressureCriterion(rSolverSettings);
        }

        if (convergence_criterion_type == "displacement_and_water_pressure_criterion") {
            return CreateAndCriterion(CreateDisplacementCriterion(rSolverSettings),
                                      CreateWaterPressureCriterion(rSolverSettings));
        }

        KRATOS_ERROR << "The convergence_criterion (" << convergence_criterion_type << ") is unknown, "
                     << "supported criteria are: 'displacement_criterion', "
                        "'residual_criterion'."
                     << std::endl;
    }

private:
    static ConvergenceCriterionSharedPtr CreateDisplacementCriterion(const Parameters& rSolverSettings)
    {
        using namespace std::string_literals;

        const auto entries_to_copy =
            std::vector{"displacement_absolute_tolerance"s, "displacement_relative_tolerance"s};
        const auto convergence_inputs =
            ParametersUtilities::CopyOptionalParameters(rSolverSettings, entries_to_copy);
        return std::make_shared<DisplacementCriteria<TSparseSpace, TDenseSpace>>(convergence_inputs);
    }

    static ConvergenceCriterionSharedPtr CreateResidualCriterion(const Parameters& rSolverSettings)
    {
        using namespace std::string_literals;

        const auto entries_to_copy = std::vector{"residual_absolute_tolerance"s, "residual_relative_tolerance"s};
        const auto convergence_inputs =
            ParametersUtilities::CopyOptionalParameters(rSolverSettings, entries_to_copy);
        return std::make_shared<ResidualCriteria<TSparseSpace, TDenseSpace>>(convergence_inputs);
    }

    static ConvergenceCriterionSharedPtr CreateWaterPressureCriterion(const Parameters& rSolverSettings)
    {
        using namespace std::string_literals;

        const auto convergence_variables = std::vector{
            std::make_tuple<const VariableData*, MixedGenericCriterionType::TDataType, MixedGenericCriterionType::TDataType>(
                &WATER_PRESSURE, rSolverSettings["water_pressure_relative_tolerance"s].GetDouble(),
                rSolverSettings["water_pressure_absolute_tolerance"s].GetDouble())};
        return std::make_shared<MixedGenericCriterionType>(convergence_variables);
    }

    static ConvergenceCriterionSharedPtr CreateAndCriterion(ConvergenceCriterionSharedPtr FirstCriterion,
                                                            ConvergenceCriterionSharedPtr SecondCriterion)
    {
        return std::make_shared<And_Criteria<TSparseSpace, TDenseSpace>>(
            std::move(FirstCriterion), std::move(SecondCriterion));
    }
};

} // namespace Kratos
