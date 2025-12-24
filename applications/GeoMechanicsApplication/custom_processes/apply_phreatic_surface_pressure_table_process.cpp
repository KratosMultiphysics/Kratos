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

#include "apply_phreatic_surface_pressure_table_process.h"
#include "geo_mechanics_application_variables.h"
#include "includes/model_part.h"

#include <string>

namespace Kratos
{
using namespace std::string_literals;

ApplyPhreaticSurfacePressureTableProcess::ApplyPhreaticSurfacePressureTableProcess(ModelPart& model_part,
                                                                                   Parameters rParameters)
    : ApplyConstantPhreaticSurfacePressureProcess(model_part, rParameters)
{
    KRATOS_TRY

    unsigned int TableId = rParameters["table"].GetInt();
    mpTable              = model_part.pGetTable(TableId);
    mTimeUnitConverter   = model_part.GetProcessInfo()[TIME_UNIT_CONVERTER];

    KRATOS_CATCH("")
}

void ApplyPhreaticSurfacePressureTableProcess::ExecuteInitializeSolutionStep()
{
    KRATOS_TRY

    const Variable<double>& var    = KratosComponents<Variable<double>>::Get(mVariableName);
    const double            Time   = mrModelPart.GetProcessInfo()[TIME] / mTimeUnitConverter;
    const double            deltaH = mpTable->GetValue(Time);

    Vector3 direction            = ZeroVector(3);
    direction[mGravityDirection] = 1.0;

    block_for_each(mrModelPart.Nodes(), [&var, &direction, &deltaH, this](Node& rNode) {
        double       distance = inner_prod(mNormalVector, rNode.Coordinates());
        const double d        = inner_prod(mNormalVector, direction);
        distance              = -(distance - mEqRHS) / d;
        distance += deltaH;
        const double pressure = -PORE_PRESSURE_SIGN_FACTOR * mSpecificWeight * distance;
        if (mIsSeepage) {
            if (pressure < PORE_PRESSURE_SIGN_FACTOR * mPressureTensionCutOff) { // Before 0. was used i.s.o. the tension cut off value -> no effect in any test.
                rNode.FastGetSolutionStepValue(var) = pressure;
                if (mIsFixed) rNode.Fix(var);
            } else {
                if (mIsFixedProvided) rNode.Free(var);
            }
        } else {
            rNode.FastGetSolutionStepValue(var) =
                std::min(pressure, PORE_PRESSURE_SIGN_FACTOR * mPressureTensionCutOff);
        }
    });

    KRATOS_CATCH("")
}

std::string ApplyPhreaticSurfacePressureTableProcess::Info() const
{
    return "ApplyPhreaticSurfacePressureTableProcess"s;
}

} // namespace Kratos