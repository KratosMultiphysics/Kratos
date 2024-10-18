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
    Model &rModel,
    Parameters Params)
    : Process(),
      mrModelPart(rModel.GetModelPart(Params["model_part_name"].GetString()))
{
    // Check default settings
    KRATOS_TRY
    
    Params.ValidateAndAssignDefaults(GetDefaultParameters());

    SelectExecutionTime(Params);
    SelectPressureGetter(Params);
    ReadFreestreamValues(Params);

    KRATOS_CATCH("")
}


const Parameters ComputePressureCoefficientProcess::GetDefaultParameters() const
{
    return Parameters( R"(
    {
        "model_part_name"     : "PLEASE_PROVIDE_A_MODELPART_NAME",
        "freestream_density" : 0,
        "freestream_velocity" : 0,
        "freestream_pressure" : 0,
        "pressure_is_historical" : true,
        "execution_step" : "ExecuteBeforeOutputStep"
    })" );
}


void ComputePressureCoefficientProcess::SelectExecutionTime(Parameters Params)
{
    const std::string execution_step = Params["execution_step"].GetString();
    if(execution_step == "ExecuteFinalizeSolutionStep") {
        mComputeAsPostProcess = false;
    }
    else if(execution_step == "ExecuteBeforeOutputStep") {
        mComputeAsPostProcess = true;
    }
    else {
        KRATOS_ERROR << "Invalid value for 'execution_step'. Try any of 'ExecuteFinalizeSolutionStep', 'ExecuteBeforeOutputStep'.";
    }
}


void ComputePressureCoefficientProcess::SelectPressureGetter(Parameters Params)
{
    if(Params["pressure_is_historical"].GetBool()) {
        mGetPressure = [](NodeType const& r_node) { return r_node.FastGetSolutionStepValue(PRESSURE); };
    } else {
        mGetPressure = [](NodeType const& r_node) { return r_node.GetValue(PRESSURE); };
    }
}


void ComputePressureCoefficientProcess::ReadFreestreamValues(Parameters Params)
{
    constexpr double epsilon = 1e-12;

    mFreestreamStaticPressure = Params["freestream_pressure"].GetDouble();

    const double freestream_density = Params["freestream_density"].GetDouble();
    KRATOS_ERROR_IF(freestream_density < epsilon) << "Value of 'freestream_density' must be greater than zero";

    const double free_stream_velocity = Params["freestream_velocity"].GetDouble();
    KRATOS_ERROR_IF(std::abs(free_stream_velocity) < epsilon) << "Value of 'free_stream_velocity' must be non-zero";

    mFreestreamDynamicPressure = 0.5 * freestream_density * free_stream_velocity * free_stream_velocity;
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
            const auto pressure = mGetPressure(r_node);
            const double cp = (pressure - mFreestreamStaticPressure) / mFreestreamDynamicPressure;
            r_node.SetValue(PRESSURE_COEFFICIENT, cp);
        }
    );

    KRATOS_CATCH("")
}


}