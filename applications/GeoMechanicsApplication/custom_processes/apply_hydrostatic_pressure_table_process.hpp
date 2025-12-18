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

#pragma once

#include "includes/table.h"

#include "custom_processes/apply_constant_hydrostatic_pressure_process.h"
#include "geo_mechanics_application_variables.h"

namespace Kratos
{

class ApplyHydrostaticPressureTableProcess : public ApplyConstantHydrostaticPressureProcess
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(ApplyHydrostaticPressureTableProcess);

    /// Defining a table with double argument and result type as table type.
    using TableType = Table<double, double>;

    ApplyHydrostaticPressureTableProcess(ModelPart& model_part, Parameters rParameters)
        : ApplyConstantHydrostaticPressureProcess(model_part, rParameters)
    {
        KRATOS_TRY

        unsigned int TableId = rParameters["table"].GetInt();
        mpTable              = model_part.pGetTable(TableId);
        mTimeUnitConverter   = model_part.GetProcessInfo()[TIME_UNIT_CONVERTER];

        KRATOS_CATCH("")
    }

    ApplyHydrostaticPressureTableProcess(const ApplyHydrostaticPressureTableProcess&) = delete;
    ApplyHydrostaticPressureTableProcess& operator=(const ApplyHydrostaticPressureTableProcess&) = delete;
    ~ApplyHydrostaticPressureTableProcess() override = default;

    /// this function will be executed at every time step BEFORE performing the solve phase
    void ExecuteInitializeSolutionStep() override
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

    /// Turn back information as a string.
    std::string Info() const override { return "ApplyHydrostaticPressureTableProcess"; }

private:
    /// Member Variables
    TableType::Pointer mpTable;
    double             mTimeUnitConverter;
};

} // namespace Kratos