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

#include "custom_processes/apply_constant_boundary_phreatic_line_pressure_process.hpp"
#include "geo_mechanics_application_variables.h"

namespace Kratos
{

class ApplyBoundaryPhreaticLinePressureTableProcess : public ApplyConstantBoundaryPhreaticLinePressureProcess
{

public:

    KRATOS_CLASS_POINTER_DEFINITION(ApplyBoundaryPhreaticLinePressureTableProcess);

    /// Defining a table with double argument and result type as table type.
    using TableType = Table<double,double>;

    ApplyBoundaryPhreaticLinePressureTableProcess(ModelPart& model_part,
                                                 Parameters rParameters
                                                 ) : ApplyConstantBoundaryPhreaticLinePressureProcess(model_part, rParameters)
    {
        KRATOS_TRY

        unsigned int TableId = rParameters["table"].GetInt();
        mpTable            = model_part.pGetTable(TableId);
        mTimeUnitConverter = model_part.GetProcessInfo()[TIME_UNIT_CONVERTER];

        KRATOS_CATCH("")
    }

    ApplyBoundaryPhreaticLinePressureTableProcess(const ApplyBoundaryPhreaticLinePressureTableProcess&) = delete;
    ApplyBoundaryPhreaticLinePressureTableProcess& operator=(const ApplyBoundaryPhreaticLinePressureTableProcess&) = delete;
    ~ApplyBoundaryPhreaticLinePressureTableProcess() override = default;

    /// this function will be executed at every time step BEFORE performing the solve phase
    void ExecuteInitializeSolutionStep() override
    {
        KRATOS_TRY

        const auto& r_variable = KratosComponents<Variable<double>>::Get(GetVariableName());
        const auto  time       = GetModelPart().GetProcessInfo()[TIME] / mTimeUnitConverter;
        const auto  delta_h    = mpTable->GetValue(time);

        block_for_each(GetModelPart().Nodes(), [&delta_h, &r_variable, this](Node& rNode) {
            auto xcoord = rNode.Coordinates()[GetHorizontalDirection()];
            xcoord      = std::max(xcoord, GetMinHorizontalCoordinate());
            xcoord      = std::min(xcoord, GetMaxHorizontalCoordinate());
            const auto height =
                GetSlope() * (xcoord - GetFirstReferenceCoordinate()[GetHorizontalDirection()]) +
                GetFirstReferenceCoordinate()[GetGravityDirection()];
            const auto distance = height - rNode.Coordinates()[GetGravityDirection()];
            const auto pressure = GetSpecificWeight() * (distance + delta_h);
            rNode.FastGetSolutionStepValue(r_variable) = std::max(pressure, 0.0);
        });

        KRATOS_CATCH("")
    }

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "ApplyBoundaryPhreaticLinePressureTableProcess";
    }

private:
    /// Member Variables
    TableType::Pointer mpTable;
    double mTimeUnitConverter;

};

}