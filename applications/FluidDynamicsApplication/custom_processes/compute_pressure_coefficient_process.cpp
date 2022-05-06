//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Eduard Gómez
//
//

// System includes

// External includes

// Project includes

// Application includes
#include "compute_pressure_coefficient_process.h"
#include "fluid_dynamics_application_variables.h"
#include "utilities/variable_utils.h"

namespace Kratos
{

ComputePressureCoefficientProcess::ComputePressureCoefficientProcess(
    ModelPart& rModelPart,
    Parameters Params)
    : Process(),
      mrModelPart(rModelPart)
{
    // Check default settings
    this->CheckDefaultsAndProcessSettings(Params);
}

ComputePressureCoefficientProcess::ComputePressureCoefficientProcess(
    Model &rModel,
    Parameters Params)
    : Process(),
      mrModelPart(rModel.GetModelPart(Params["model_part_name"].GetString()))
{
    // Check default settings
    this->CheckDefaultsAndProcessSettings(Params);
}


void ComputePressureCoefficientProcess::CheckDefaultsAndProcessSettings(Parameters Params)
{
    KRATOS_TRY

    Parameters default_parameters( R"(
    {
        "model_part_name"     : "PLEASE_PROVIDE_A_MODELPART_NAME",
        "freestream_density" : 0,
        "freestream_velocity" : 0,
        "freestream_pressure" : 0,
        "compute_only_as_postprocess" : true
    })" );

    Params.ValidateAndAssignDefaults(default_parameters);

    mComputeAsPostProcess = Params["compute_only_as_postprocess"].GetBool();

    mFreestreamStaticPressure = Params["freestream_pressure"].GetDouble();

    const double freestream_density = Params["freestream_density"].GetDouble();
    const double free_stream_velocity = Params["freestream_velocity"].GetDouble();
    mFreestreamDynamicPressure = 0.5 * freestream_density * free_stream_velocity * free_stream_velocity;

    KRATOS_CATCH("")
}


void ComputePressureCoefficientProcess::ExecuteInitialize()
{
    VariableUtils().SetNonHistoricalVariableToZero(PRESSURE_COEFFICIENT, mrModelPart.Nodes());
}

void ComputePressureCoefficientProcess::ExecuteFinalizeSolutionStep()
{
    if(!mComputeAsPostProcess)
    {
        Execute();
    }
}

void ComputePressureCoefficientProcess::ExecuteBeforeOutputStep()
{
    if(mComputeAsPostProcess)
    {
        Execute();
    }   
}


void ComputePressureCoefficientProcess::Execute()
{
    KRATOS_TRY

    block_for_each(mrModelPart.Nodes(),
        [&](NodeType& r_node)
        {
            const auto pressure = r_node.GetSolutionStepValue(PRESSURE);
            const double cp = (pressure - mFreestreamStaticPressure) / mFreestreamDynamicPressure;
            r_node.SetValue(PRESSURE_COEFFICIENT, cp);
        }
    );

    KRATOS_CATCH("")
}


}