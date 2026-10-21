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

#include "custom_utilities/geomechanics_residual_criteria.hpp"
#include "custom_utilities/string_utilities.h"
#include "includes/variables.h" // for WATER_PRESSURE
#include "parameters_utilities.h"
#include "solving_strategies/convergencecriterias/and_criteria.h"
#include "solving_strategies/convergencecriterias/displacement_criteria.h"
#include "solving_strategies/convergencecriterias/mixed_generic_criteria.h"
#include "solving_strategies/convergencecriterias/or_criteria.h"

#include <algorithm>
#include <functional>
#include <map>
#include <string>

using namespace std::string_literals;

namespace Kratos
{

template <class TSparseSpace, class TDenseSpace>
class ConvergenceCriteriaFactory
{
public:
    using ConvergenceCriterionType      = ConvergenceCriteria<TSparseSpace, TDenseSpace>;
    using ConvergenceCriterionSharedPtr = std::shared_ptr<ConvergenceCriterionType>;
    using DisplacementCriterionType     = DisplacementCriteria<TSparseSpace, TDenseSpace>;
    using ResidualCriterionType         = GeoMechanicsResidualCriteria<TSparseSpace, TDenseSpace>;
    using AndCriterionType              = And_Criteria<TSparseSpace, TDenseSpace>;
    using OrCriterionType               = Or_Criteria<TSparseSpace, TDenseSpace>;
    using MixedGenericCriterionType     = MixedGenericCriteria<TSparseSpace, TDenseSpace>;

    static ConvergenceCriterionSharedPtr Create(const Parameters& rSolverSettings)
    {
        KRATOS_ERROR_IF_NOT(rSolverSettings.Has("convergence_criterion"s))
            << "No convergence_criterion is defined, aborting.";

        // The following map resembles the switch in method `_ConstructConvergenceCriterion` (see
        // `geomechanics_U_Pw_solver.py`)
        using Creator = std::function<ConvergenceCriterionSharedPtr(const Parameters&)>;
        static const auto creator_map = std::map<std::string, Creator, std::less<>>{
            {"displacement_criterion"s, CreateDisplacementCriterion},
            {"residual_criterion"s, CreateResidualCriterion},
            {"and_criterion"s,
             [](const Parameters& rSolverSettings) {
            return CreateAndCriterion(CreateResidualCriterion(rSolverSettings),
                                      CreateDisplacementCriterion(rSolverSettings));
        }},
            {"or_criterion"s,
             [](const Parameters& rSolverSettings) {
            return CreateOrCriterion(CreateResidualCriterion(rSolverSettings),
                                     CreateDisplacementCriterion(rSolverSettings));
        }},
            {"water_pressure_criterion"s, CreateWaterPressureCriterion},
            {"displacement_and_water_pressure_criterion"s, [](const Parameters& rSolverSettings) {
            return CreateAndCriterion(CreateDisplacementCriterion(rSolverSettings),
                                      CreateWaterPressureCriterion(rSolverSettings));
        }}};

        const auto convergence_criterion_type =
            GeoStringUtilities::ToLower(rSolverSettings["convergence_criterion"s].GetString());

        if (!creator_map.contains(convergence_criterion_type)) {
            auto supported_criteria = std::vector<std::string>{};
            std::ranges::transform(creator_map, std::back_inserter(supported_criteria),
                                   [](const auto& rKeyValuePair) { return rKeyValuePair.first; });
            KRATOS_ERROR << "The convergence_criterion (" << convergence_criterion_type << ") is unknown, "
                         << "supported criteria are: " << GeoStringUtilities::Join(supported_criteria, ", "s)
                         << std::endl;
        }

        return creator_map.at(convergence_criterion_type)(rSolverSettings);
    }

private:
    static ConvergenceCriterionSharedPtr CreateDisplacementCriterion(const Parameters& rSolverSettings)
    {
        const auto entries_to_copy =
            std::vector{"displacement_absolute_tolerance"s, "displacement_relative_tolerance"s};
        const auto convergence_inputs =
            ParametersUtilities::CopyOptionalParameters(rSolverSettings, entries_to_copy);
        return std::make_shared<DisplacementCriterionType>(convergence_inputs);
    }

    static ConvergenceCriterionSharedPtr CreateResidualCriterion(const Parameters& rSolverSettings)
    {
        const auto entries_to_copy = std::vector{"residual_absolute_tolerance"s, "residual_relative_tolerance"s};
        const auto convergence_inputs =
            ParametersUtilities::CopyOptionalParameters(rSolverSettings, entries_to_copy);
        auto residual_representation = ResidualCriterionType::ResidualRepresentation::Raw;
        if (rSolverSettings.Has("strategy_type"s) &&
            GeoStringUtilities::ToLower(rSolverSettings["strategy_type"s].GetString()) == "quasi_newton"s) {
            residual_representation = ResidualCriterionType::ResidualRepresentation::ConstraintReduced;
        }
        if (rSolverSettings.Has("residual_representation"s)) {
            const auto representation =
                GeoStringUtilities::ToLower(rSolverSettings["residual_representation"s].GetString());
            if (representation == "raw"s) {
                residual_representation = ResidualCriterionType::ResidualRepresentation::Raw;
            } else if (representation == "constraint_reduced"s) {
                residual_representation = ResidualCriterionType::ResidualRepresentation::ConstraintReduced;
            } else {
                KRATOS_ERROR << "Unknown residual_representation (" << representation
                             << "), supported values are: raw, constraint_reduced." << std::endl;
            }
        }
        return std::make_shared<ResidualCriterionType>(convergence_inputs, residual_representation);
    }

    static ConvergenceCriterionSharedPtr CreateWaterPressureCriterion(const Parameters& rSolverSettings)
    {
        const auto convergence_variables = std::vector{
            std::make_tuple<const VariableData*, typename MixedGenericCriterionType::TDataType, typename MixedGenericCriterionType::TDataType>(
                &WATER_PRESSURE, rSolverSettings["water_pressure_relative_tolerance"s].GetDouble(),
                rSolverSettings["water_pressure_absolute_tolerance"s].GetDouble())};
        return std::make_shared<MixedGenericCriterionType>(convergence_variables);
    }

    static ConvergenceCriterionSharedPtr CreateAndCriterion(ConvergenceCriterionSharedPtr pFirstCriterion,
                                                            ConvergenceCriterionSharedPtr pSecondCriterion)
    {
        return std::make_shared<AndCriterionType>(std::move(pFirstCriterion), std::move(pSecondCriterion));
    }

    static ConvergenceCriterionSharedPtr CreateOrCriterion(ConvergenceCriterionSharedPtr pFirstCriterion,
                                                           ConvergenceCriterionSharedPtr pSecondCriterion)
    {
        return std::make_shared<OrCriterionType>(std::move(pFirstCriterion), std::move(pSecondCriterion));
    }
};

} // namespace Kratos
