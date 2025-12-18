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

#include "apply_phreatic_line_pressure_table_process.h"
#include "geo_mechanics_application_variables.h"
#include "includes/model_part.h"

namespace Kratos
{

ApplyPhreaticLinePressureTableProcess::ApplyPhreaticLinePressureTableProcess(ModelPart& model_part, Parameters rParameters)
    : ApplyConstantPhreaticLinePressureProcess(model_part, rParameters)
{
    KRATOS_TRY

    for (unsigned int i = 0; i < mpTable.size(); ++i) {
        unsigned int TableId = rParameters["table"][i].GetInt();
        if (TableId > 0) {
            mpTable[i] = model_part.pGetTable(TableId);
        } else {
            mpTable[i] = nullptr;
        }
    }

    mTimeUnitConverter = model_part.GetProcessInfo()[TIME_UNIT_CONVERTER];

    KRATOS_CATCH("")
}

void ApplyPhreaticLinePressureTableProcess::ExecuteInitializeSolutionStep()
{
    KRATOS_TRY

    const Variable<double>& var = KratosComponents<Variable<double>>::Get(mVariableName);

    const double        Time = mrModelPart.GetProcessInfo()[TIME] / mTimeUnitConverter;
    array_1d<double, 2> deltaH;
    for (unsigned int i = 0; i < mpTable.size(); ++i) {
        if (mpTable[i]) {
            deltaH[i] = mpTable[i]->GetValue(Time);
        } else {
            deltaH[i] = 0.0;
        }
    }

    array_1d<double, 2> y;
    y[0] = mFirstReferenceCoordinate[mGravityDirection];
    y[1] = mSecondReferenceCoordinate[mGravityDirection];
    y += deltaH;

    mSlope = (y[1] - y[0]) / (mSecondReferenceCoordinate[mHorizontalDirection] -
                              mFirstReferenceCoordinate[mHorizontalDirection]);

    block_for_each(mrModelPart.Nodes(), [&var, &y, this](Node& rNode) {
        const double pressure = PORE_PRESSURE_SIGN_FACTOR * CalculatePressurewithTable(rNode, y);
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

std::string ApplyPhreaticLinePressureTableProcess::Info() const
{
    return "ApplyPhreaticLinePressureTableProcess";
}

double ApplyPhreaticLinePressureTableProcess::CalculatePressurewithTable(const Node& rNode,
                                                                         const array_1d<double, 2>& y) const
{
    double horCoord = rNode.Coordinates()[mHorizontalDirection];
    horCoord        = std::max(horCoord, mMinHorizontalCoordinate);
    horCoord        = std::min(horCoord, mMaxHorizontalCoordinate);
    const double height = mSlope * (horCoord - mFirstReferenceCoordinate[mHorizontalDirection]) + y[0];
    const double distance = height - rNode.Coordinates()[mGravityDirection];
    return -mSpecificWeight * distance;
}

} // namespace Kratos