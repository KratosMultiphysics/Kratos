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
    typedef Table<double,double> TableType;

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

        const Variable<double> &var = KratosComponents< Variable<double> >::Get(mVariableName);
        const double Time   = mrModelPart.GetProcessInfo()[TIME]/mTimeUnitConverter;
        const double deltaH = mpTable->GetValue(Time);

        block_for_each(mrModelPart.Nodes(), [&deltaH, &var, this](Node& rNode){
            double xcoord = rNode.Coordinates()[mHorizontalDirection];
            xcoord = std::max(xcoord,mMinHorizontalCoordinate);
            xcoord = std::min(xcoord,mMaxHorizontalCoordinate);
            double height = mSlope * ( xcoord - mFirstReferenceCoordinate[mHorizontalDirection]) + mFirstReferenceCoordinate[mGravityDirection];
            const double distance = height - rNode.Coordinates()[mGravityDirection];
            const double pressure = mSpecificWeight * (distance + deltaH);
            rNode.FastGetSolutionStepValue(var) = std::max(pressure,0.0);
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