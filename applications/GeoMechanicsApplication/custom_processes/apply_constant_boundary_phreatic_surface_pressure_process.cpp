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

#include "apply_constant_boundary_phreatic_surface_pressure_process.h"
#include "includes/model_part.h"
#include "utilities/math_utils.h"

namespace Kratos
{

ApplyConstantBoundaryPhreaticSurfacePressureProcess::ApplyConstantBoundaryPhreaticSurfacePressureProcess(
    ModelPart& model_part, Parameters rParameters)
    : Process(Flags()), mrModelPart(model_part)
{
    KRATOS_TRY

    // only include validation with c++11 since raw_literals do not exist in c++03
    Parameters default_parameters(R"(
            {
                "model_part_name":"PLEASE_CHOOSE_MODEL_PART_NAME",
                "variable_name": "PLEASE_PRESCRIBE_VARIABLE_NAME",
                "is_fixed": false,
                "gravity_direction": 1,
                "first_reference_coordinate":           [0.0,1.0,0.0],
                "second_reference_coordinate":          [1.0,0.5,0.0],
                "third_reference_coordinate":           [1.0,0.5,1.0],
                "specific_weight" : 10000.0,
                "table" : 1
            }  )");

    // Some values need to be mandatorily prescribed since no meaningful default value exist.
    // For this reason try accessing to them So that an error is thrown if they don't exist
    rParameters["first_reference_coordinate"];
    rParameters["second_reference_coordinate"];
    rParameters["third_reference_coordinate"];
    rParameters["variable_name"];
    rParameters["model_part_name"];
    mIsFixedProvided = rParameters.Has("is_fixed");

    // Now validate against defaults -- this also ensures no type mismatch
    rParameters.ValidateAndAssignDefaults(default_parameters);

    mVariableName              = rParameters["variable_name"].GetString();
    mIsFixed                   = rParameters["is_fixed"].GetBool();
    mGravityDirection          = rParameters["gravity_direction"].GetInt();
    mFirstReferenceCoordinate  = rParameters["first_reference_coordinate"].GetVector();
    mSecondReferenceCoordinate = rParameters["second_reference_coordinate"].GetVector();
    mThirdReferenceCoordinate  = rParameters["third_reference_coordinate"].GetVector();

    calculateEquationParameters();

    mSpecificWeight = rParameters["specific_weight"].GetDouble();

    KRATOS_CATCH("")
}

void ApplyConstantBoundaryPhreaticSurfacePressureProcess::ExecuteInitialize()
{
    KRATOS_TRY

    const Variable<double>& var = KratosComponents<Variable<double>>::Get(mVariableName);

    Vector3 direction            = ZeroVector(3);
    direction[mGravityDirection] = 1.0;

    block_for_each(mrModelPart.Nodes(), [&var, &direction, this](Node& rNode) {
        if (mIsFixed) rNode.Fix(var);
        else if (mIsFixedProvided) rNode.Free(var);

        double       distance               = inner_prod(mNormalVector, rNode.Coordinates());
        const double d                      = inner_prod(mNormalVector, direction);
        distance                            = -(distance - mEqRHS) / d;
        const double pressure               = mSpecificWeight * distance;
        rNode.FastGetSolutionStepValue(var) = std::max(pressure, 0.0);
    });

    KRATOS_CATCH("")
}

std::string ApplyConstantBoundaryPhreaticSurfacePressureProcess::Info() const
{
    return "ApplyConstantBoundaryPhreaticSurfacePressureProcess";
}

void ApplyConstantBoundaryPhreaticSurfacePressureProcess::calculateEquationParameters()
{
    const Vector3 v1 = mSecondReferenceCoordinate - mFirstReferenceCoordinate;
    const Vector3 v2 = mThirdReferenceCoordinate - mFirstReferenceCoordinate;
    mNormalVector    = MathUtils<>::CrossProduct(v1, v2);
    if (norm_2(mNormalVector) == 0.0)
        KRATOS_ERROR << "Normal vector to phreatic surface has zero size!" << std::endl;

    mEqRHS = inner_prod(mNormalVector, mFirstReferenceCoordinate);
}

} // namespace Kratos