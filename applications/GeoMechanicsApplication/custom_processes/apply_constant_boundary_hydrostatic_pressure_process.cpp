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

#include "apply_constant_boundary_hydrostatic_pressure_process.h"
#include "geo_mechanics_application_variables.h"
#include "includes/model_part.h"

namespace Kratos
{

ApplyConstantBoundaryHydrostaticPressureProcess::ApplyConstantBoundaryHydrostaticPressureProcess(ModelPart& model_part,
                                                                                                 Parameters rParameters)
    : Process(Flags()), mrModelPart(model_part)
{
    KRATOS_TRY

    // only include validation with c++11 since raw_literals do not exist in c++03
    Parameters default_parameters(R"(
            {
                "model_part_name":"PLEASE_CHOOSE_MODEL_PART_NAME",
                "variable_name": "PLEASE_PRESCRIBE_VARIABLE_NAME",
                "is_fixed": false,
                "gravity_direction" : 2,
                "reference_coordinate" : 0.0,
                "specific_weight" : 10000.0,
                "table" : 1
            }  )");

    // Some values have to be input by the user since no meaningful default value exist. For
    // this reason, we try to access them, so that an error is thrown if they don't exist.
    rParameters["reference_coordinate"];
    rParameters["variable_name"];
    rParameters["model_part_name"];

    mIsFixedProvided = rParameters.Has("is_fixed");

    // Now validate against defaults -- this also ensures no type mismatch
    rParameters.ValidateAndAssignDefaults(default_parameters);

    mVariableName        = rParameters["variable_name"].GetString();
    mIsFixed             = rParameters["is_fixed"].GetBool();
    mGravityDirection    = static_cast<unsigned int>(rParameters["gravity_direction"].GetInt());
    mReferenceCoordinate = rParameters["reference_coordinate"].GetDouble();
    mSpecificWeight      = rParameters["specific_weight"].GetDouble();

    KRATOS_CATCH("")
}

void ApplyConstantBoundaryHydrostaticPressureProcess::ExecuteInitialize()
{
    KRATOS_TRY

    const auto& r_variable = KratosComponents<Variable<double>>::Get(GetVariableName());

    block_for_each(GetModelPart().Nodes(), [&r_variable, this](Node& rNode) {
        if (mIsFixed) rNode.Fix(r_variable);
        else if (mIsFixedProvided) rNode.Free(r_variable);

        const auto pressure = GetSpecificWeight() *
                              (GetReferenceCoordinate() - rNode.Coordinates()[GetGravityDirection()]);
        rNode.FastGetSolutionStepValue(r_variable) = std::max(pressure, 0.);
    });

    KRATOS_CATCH("")
}

std::string ApplyConstantBoundaryHydrostaticPressureProcess::Info() const
{
    return "ApplyConstantBoundaryHydrostaticPressureProcess";
}

ModelPart& ApplyConstantBoundaryHydrostaticPressureProcess::GetModelPart() { return mrModelPart; }

const std::string& ApplyConstantBoundaryHydrostaticPressureProcess::GetVariableName() const
{
    return mVariableName;
}

unsigned int ApplyConstantBoundaryHydrostaticPressureProcess::GetGravityDirection() const
{
    return mGravityDirection;
}

double ApplyConstantBoundaryHydrostaticPressureProcess::GetReferenceCoordinate() const
{
    return mReferenceCoordinate;
}

double ApplyConstantBoundaryHydrostaticPressureProcess::GetSpecificWeight() const
{
    return mSpecificWeight;
}

} // namespace Kratos