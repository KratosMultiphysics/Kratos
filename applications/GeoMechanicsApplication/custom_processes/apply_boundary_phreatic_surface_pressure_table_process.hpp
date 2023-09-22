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

#include "custom_processes/apply_constant_boundary_phreatic_surface_pressure_process.hpp"
#include "geo_mechanics_application_variables.h"

namespace Kratos
{

class ApplyBoundaryPhreaticSurfacePressureTableProcess : public ApplyConstantBoundaryPhreaticSurfacePressureProcess
{

public:

    KRATOS_CLASS_POINTER_DEFINITION(ApplyBoundaryPhreaticSurfacePressureTableProcess);

    /// Defining a table with double argument and result type as table type.
    using TableType = Table<double,double>;

    ApplyBoundaryPhreaticSurfacePressureTableProcess(ModelPart& model_part,
                                                 Parameters rParameters
                                                 ) : ApplyConstantBoundaryPhreaticSurfacePressureProcess(model_part, rParameters)
    {
        KRATOS_TRY

        unsigned int TableId = rParameters["table"].GetInt();
        mpTable = model_part.pGetTable(TableId);
        mTimeUnitConverter = model_part.GetProcessInfo()[TIME_UNIT_CONVERTER];

        KRATOS_CATCH("")
    }

    ApplyBoundaryPhreaticSurfacePressureTableProcess(const ApplyBoundaryPhreaticSurfacePressureTableProcess&) = delete;
    ApplyBoundaryPhreaticSurfacePressureTableProcess& operator=(const ApplyBoundaryPhreaticSurfacePressureTableProcess&) = delete;
    ~ApplyBoundaryPhreaticSurfacePressureTableProcess() override = default;

    /// this function will be executed at every time step BEFORE performing the solve phase
    void ExecuteInitializeSolutionStep() override
    {
        KRATOS_TRY

        const Variable<double> &var = KratosComponents< Variable<double> >::Get(mVariableName);
        const double Time   = mrModelPart.GetProcessInfo()[TIME] / mTimeUnitConverter;
        const double deltaH = mpTable->GetValue(Time);

        Vector3 direction = ZeroVector(3);
        direction[mGravityDirection] = 1.0;

        block_for_each(mrModelPart.Nodes(), [&deltaH, &var, &direction, this](Node& rNode){
            double distance = inner_prod(mNormalVector, rNode.Coordinates());
            const double d  = inner_prod(mNormalVector, direction);
            distance = -(distance - mEqRHS) / d;
            const double pressure = mSpecificWeight * ( distance + deltaH );
            rNode.FastGetSolutionStepValue(var) = std::max(pressure,0.0);
        });

        KRATOS_CATCH("")
    }

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "ApplyBoundaryPhreaticSurfacePressureTableProcess";
    }

private:
    /// Member Variables
    TableType::Pointer mpTable;
    double mTimeUnitConverter;

};

}