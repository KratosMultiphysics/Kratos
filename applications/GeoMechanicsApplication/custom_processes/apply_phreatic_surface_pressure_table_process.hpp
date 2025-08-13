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

#include "custom_processes/apply_constant_phreatic_surface_pressure_process.hpp"
#include "geo_mechanics_application_variables.h"

namespace Kratos
{

class ApplyPhreaticSurfacePressureTableProcess : public ApplyConstantPhreaticSurfacePressureProcess
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(ApplyPhreaticSurfacePressureTableProcess);

    /// Defining a table with double argument and result type as table type.
    using TableType = Table<double, double>;

    ApplyPhreaticSurfacePressureTableProcess(ModelPart& model_part, Parameters rParameters)
        : ApplyConstantPhreaticSurfacePressureProcess(model_part, rParameters)
    {
        KRATOS_TRY

        unsigned int TableId = rParameters["table"].GetInt();
        mpTable              = model_part.pGetTable(TableId);
        mTimeUnitConverter   = model_part.GetProcessInfo()[TIME_UNIT_CONVERTER];

        KRATOS_CATCH("")
    }

    ApplyPhreaticSurfacePressureTableProcess(const ApplyPhreaticSurfacePressureTableProcess&) = delete;
    ApplyPhreaticSurfacePressureTableProcess& operator=(const ApplyPhreaticSurfacePressureTableProcess&) = delete;
    ~ApplyPhreaticSurfacePressureTableProcess() override = default;

    /// this function will be executed at every time step BEFORE performing the solve phase
    void ExecuteInitializeSolutionStep() override
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

    /// Turn back information as a string.
    std::string Info() const override { return "ApplyPhreaticSurfacePressureTableProcess"; }

private:
    TableType::Pointer mpTable;
    double             mTimeUnitConverter;
};

} // namespace Kratos