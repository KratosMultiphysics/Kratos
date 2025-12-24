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

#include "apply_hydrostatic_pressure_table_process.h"
#include "geo_mechanics_application_variables.h"
#include "includes/model_part.h"

#include <string>

namespace Kratos
{
using namespace std::string_literals;

ApplyHydrostaticPressureTableProcess::ApplyHydrostaticPressureTableProcess(ModelPart& model_part, Parameters rParameters)
    : ApplyConstantHydrostaticPressureProcess(model_part, rParameters)
{
    KRATOS_TRY

    unsigned int TableId = rParameters["table"].GetInt();
    mpTable              = model_part.pGetTable(TableId);
    mTimeUnitConverter   = model_part.GetProcessInfo()[TIME_UNIT_CONVERTER];

    KRATOS_CATCH("")
}

void ApplyHydrostaticPressureTableProcess::ExecuteInitializeSolutionStep()
{
    KRATOS_TRY

    const Variable<double>& var    = KratosComponents<Variable<double>>::Get(mVariableName);
    const double            Time   = mrModelPart.GetProcessInfo()[TIME] / mTimeUnitConverter;
    const double            deltaH = mpTable->GetValue(Time);

    block_for_each(mrModelPart.Nodes(), [&var, &deltaH, this](Node& rNode) {
        const double distance = mReferenceCoordinate - rNode.Coordinates()[mGravityDirection];
        const double pressure = -PORE_PRESSURE_SIGN_FACTOR * mSpecificWeight * (distance + deltaH);
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

std::string ApplyHydrostaticPressureTableProcess::Info() const
{
    return "ApplyHydrostaticPressureTableProcess"s;
}

} // namespace Kratos