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

#include "custom_processes/apply_boundary_phreatic_surface_pressure_table_process.h"
#include "geo_mechanics_application_variables.h"
#include "includes/model_part.h"

#include <string>

namespace Kratos
{
using namespace std::string_literals;

ApplyBoundaryPhreaticSurfacePressureTableProcess::ApplyBoundaryPhreaticSurfacePressureTableProcess(
    ModelPart& model_part, Parameters rParameters)
    : ApplyConstantBoundaryPhreaticSurfacePressureProcess(model_part, rParameters)
{
    KRATOS_TRY

    unsigned int TableId = rParameters["table"].GetInt();
    mpTable              = model_part.pGetTable(TableId);
    mTimeUnitConverter   = model_part.GetProcessInfo()[TIME_UNIT_CONVERTER];

    KRATOS_CATCH("")
}

void ApplyBoundaryPhreaticSurfacePressureTableProcess::ExecuteInitializeSolutionStep()
{
    KRATOS_TRY

    const Variable<double>& var    = KratosComponents<Variable<double>>::Get(mVariableName);
    const double            Time   = mrModelPart.GetProcessInfo()[TIME] / mTimeUnitConverter;
    const double            deltaH = mpTable->GetValue(Time);

    Vector3 direction            = ZeroVector(3);
    direction[mGravityDirection] = 1.0;

    block_for_each(mrModelPart.Nodes(), [&deltaH, &var, &direction, this](Node& rNode) {
        double       distance               = inner_prod(mNormalVector, rNode.Coordinates());
        const double d                      = inner_prod(mNormalVector, direction);
        distance                            = -(distance - mEqRHS) / d;
        const double pressure               = mSpecificWeight * (distance + deltaH);
        rNode.FastGetSolutionStepValue(var) = std::max(pressure, 0.0);
    });

    KRATOS_CATCH("")
}

std::string ApplyBoundaryPhreaticSurfacePressureTableProcess::Info() const
{
    return "ApplyBoundaryPhreaticSurfacePressureTableProcess"s;
}

} // namespace Kratos