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

#pragma once

#include "includes/table.h"

#include "custom_processes/apply_constant_phreatic_line_pressure_process.hpp"
#include "geo_mechanics_application_variables.h"

namespace Kratos
{

class ApplyPhreaticLinePressureTableProcess : public ApplyConstantPhreaticLinePressureProcess
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(ApplyPhreaticLinePressureTableProcess);

    /// Defining a table with double argument and result type as table type.
    using TableType = Table<double, double>;

    ApplyPhreaticLinePressureTableProcess(ModelPart& model_part, Parameters rParameters)
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

    ApplyPhreaticLinePressureTableProcess(const ApplyPhreaticLinePressureTableProcess&) = delete;
    ApplyPhreaticLinePressureTableProcess& operator=(const ApplyPhreaticLinePressureTableProcess&) = delete;
    ~ApplyPhreaticLinePressureTableProcess() override = default;

    /// this function will be executed at every time step BEFORE performing the solve phase
    void ExecuteInitializeSolutionStep() override
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

    /// Turn back information as a string.
    std::string Info() const override { return "ApplyPhreaticLinePressureTableProcess"; }

protected:
    double CalculatePressurewithTable(const Node& rNode, const array_1d<double, 2>& y) const
    {
        double horCoord = rNode.Coordinates()[mHorizontalDirection];
        horCoord        = std::max(horCoord, mMinHorizontalCoordinate);
        horCoord        = std::min(horCoord, mMaxHorizontalCoordinate);
        const double height = mSlope * (horCoord - mFirstReferenceCoordinate[mHorizontalDirection]) + y[0];
        const double distance = height - rNode.Coordinates()[mGravityDirection];
        return -mSpecificWeight * distance;
    }

private:
    /// Member Variables
    array_1d<TableType::Pointer, 2> mpTable;
    double                          mTimeUnitConverter;
};

} // namespace Kratos