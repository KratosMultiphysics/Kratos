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

#include "custom_processes/apply_constant_boundary_hydrostatic_pressure_process.hpp"
#include "geo_mechanics_application_variables.h"

namespace Kratos
{

class ApplyBoundaryHydrostaticPressureTableProcess : public ApplyConstantBoundaryHydrostaticPressureProcess
{

public:

    KRATOS_CLASS_POINTER_DEFINITION(ApplyBoundaryHydrostaticPressureTableProcess);

    /// Defining a table with double argument and result type as table type.
    using TableType = Table<double,double>;

    ApplyBoundaryHydrostaticPressureTableProcess(ModelPart& model_part,
                                                 Parameters rParameters
                                                 ) : ApplyConstantBoundaryHydrostaticPressureProcess(model_part, rParameters)
    {
        KRATOS_TRY

        unsigned int TableId = rParameters["table"].GetInt();
        mpTable = model_part.pGetTable(TableId);
        mTimeUnitConverter = model_part.GetProcessInfo()[TIME_UNIT_CONVERTER];

        KRATOS_CATCH("")
    }

    ApplyBoundaryHydrostaticPressureTableProcess(const ApplyBoundaryHydrostaticPressureTableProcess&) = delete;
    ApplyBoundaryHydrostaticPressureTableProcess& operator=(const ApplyBoundaryHydrostaticPressureTableProcess&) = delete;
    ~ApplyBoundaryHydrostaticPressureTableProcess() override = default;

    /// this function will be executed at every time step BEFORE performing the solve phase
    void ExecuteInitializeSolutionStep() override
    {
        KRATOS_TRY

        const Variable<double> &var = KratosComponents< Variable<double> >::Get(mVariableName);
        const double Time   = mrModelPart.GetProcessInfo()[TIME]/mTimeUnitConverter;
        const double deltaH = mpTable->GetValue(Time);

        block_for_each(mrModelPart.Nodes(), [&deltaH, &var, this](Node& rNode){
            const double distance = mReferenceCoordinate - rNode.Coordinates()[mGravityDirection];
            const double pressure = mSpecificWeight * (distance + deltaH);
            rNode.FastGetSolutionStepValue(var) = std::max(pressure,0.0);
        });
        KRATOS_CATCH("")
    }

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "ApplyBoundaryHydrostaticPressureTableProcess";
    }

private:
    /// Member Variables
    TableType::Pointer mpTable;
    double mTimeUnitConverter;

};

}