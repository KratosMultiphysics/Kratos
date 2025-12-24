// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Ignasi de Pouplana,
//                   Vahid Galavi
//

#include "custom_processes/apply_boundary_hydrostatic_pressure_table_process.h"
#include "geo_mechanics_application_variables.h"
#include "includes/model_part.h"

#include <string>

namespace Kratos
{
using namespace std::string_literals;

ApplyBoundaryHydrostaticPressureTableProcess::ApplyBoundaryHydrostaticPressureTableProcess(ModelPart& model_part,
                                                                                           Parameters rParameters)
    : ApplyConstantBoundaryHydrostaticPressureProcess(model_part, rParameters)
{
    KRATOS_TRY

    unsigned int TableId = rParameters["table"].GetInt();
    mpTable              = model_part.pGetTable(TableId);
    mTimeUnitConverter   = model_part.GetProcessInfo()[TIME_UNIT_CONVERTER];

    KRATOS_CATCH("")
}

void ApplyBoundaryHydrostaticPressureTableProcess::ExecuteInitializeSolutionStep()
{
    KRATOS_TRY

    const auto& r_variable = KratosComponents<Variable<double>>::Get(GetVariableName());
    const auto  time       = GetModelPart().GetProcessInfo()[TIME] / mTimeUnitConverter;
    const auto  delta_h    = mpTable->GetValue(time);

    block_for_each(GetModelPart().Nodes(), [&delta_h, &r_variable, this](Node& rNode) {
        const auto distance = GetReferenceCoordinate() - rNode.Coordinates()[GetGravityDirection()];
        const auto pressure = GetSpecificWeight() * (distance + delta_h);
        rNode.FastGetSolutionStepValue(r_variable) = std::max(pressure, 0.0);
    });
    KRATOS_CATCH("")
}

std::string ApplyBoundaryHydrostaticPressureTableProcess::Info() const
{
    return "ApplyBoundaryHydrostaticPressureTableProcess"s;
}

} // namespace Kratos