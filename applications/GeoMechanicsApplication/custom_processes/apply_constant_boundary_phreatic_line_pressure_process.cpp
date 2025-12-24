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

#include "apply_constant_boundary_phreatic_line_pressure_process.h"
#include "includes/kratos_flags.h"
#include "includes/model_part.h"

namespace Kratos
{
using namespace std::string_literals;

ApplyConstantBoundaryPhreaticLinePressureProcess::ApplyConstantBoundaryPhreaticLinePressureProcess(
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
                "out_of_plane_direction": 2,
                "first_reference_coordinate":           [0.0,1.0,0.0],
                "second_reference_coordinate":          [1.0,0.5,0.0],
                "specific_weight" : 10000.0,
                "table" : 1
            }  )");

    // Some values need to be mandatory prescribed since no meaningful default value exist. For
    // this reason try accessing to them So that an error is thrown if they don't exist
    rParameters["first_reference_coordinate"];
    rParameters["second_reference_coordinate"];
    rParameters["variable_name"];
    rParameters["model_part_name"];

    mIsFixedProvided = rParameters.Has("is_fixed");

    // Now validate against defaults -- this also ensures no type mismatch
    rParameters.ValidateAndAssignDefaults(default_parameters);

    mVariableName        = rParameters["variable_name"].GetString();
    mIsFixed             = rParameters["is_fixed"].GetBool();
    mGravityDirection    = rParameters["gravity_direction"].GetInt();
    mOutOfPlaneDirection = rParameters["out_of_plane_direction"].GetInt();
    KRATOS_ERROR_IF(mGravityDirection == mOutOfPlaneDirection)
        << "Gravity direction cannot be the same as Out-of-Plane directions " << rParameters << std::endl;

    for (unsigned int i = 0; i < N_DIM_3D; ++i)
        if (i != mGravityDirection && i != mOutOfPlaneDirection) mHorizontalDirection = i;

    mFirstReferenceCoordinate  = rParameters["first_reference_coordinate"].GetVector();
    mSecondReferenceCoordinate = rParameters["second_reference_coordinate"].GetVector();

    mMinHorizontalCoordinate = std::min(mFirstReferenceCoordinate[mHorizontalDirection],
                                        mSecondReferenceCoordinate[mHorizontalDirection]);
    mMaxHorizontalCoordinate = std::max(mFirstReferenceCoordinate[mHorizontalDirection],
                                        mSecondReferenceCoordinate[mHorizontalDirection]);

    KRATOS_ERROR_IF_NOT(mMaxHorizontalCoordinate > mMinHorizontalCoordinate)
        << "First and second point on the phreatic line have the same horizontal coordinate"
        << rParameters << std::endl;

    mSlope =
        (mSecondReferenceCoordinate[mGravityDirection] - mFirstReferenceCoordinate[mGravityDirection]) /
        (mSecondReferenceCoordinate[mHorizontalDirection] - mFirstReferenceCoordinate[mHorizontalDirection]);

    mSpecificWeight = rParameters["specific_weight"].GetDouble();

    KRATOS_CATCH("")
}

void ApplyConstantBoundaryPhreaticLinePressureProcess::ExecuteInitialize()
{
    KRATOS_TRY

    const Variable<double>& var = KratosComponents<Variable<double>>::Get(mVariableName);

    block_for_each(mrModelPart.Nodes(), [&var, this](Node& rNode) {
        if (mIsFixed) rNode.Fix(var);
        else if (mIsFixedProvided) rNode.Free(var);

        double horCoord = rNode.Coordinates()[mHorizontalDirection];
        horCoord        = std::max(horCoord, mMinHorizontalCoordinate);
        horCoord        = std::min(horCoord, mMaxHorizontalCoordinate);
        const double height = mSlope * (horCoord - mFirstReferenceCoordinate[mHorizontalDirection]) +
                              mFirstReferenceCoordinate[mGravityDirection];
        const double distance               = height - rNode.Coordinates()[mGravityDirection];
        const double pressure               = mSpecificWeight * distance;
        rNode.FastGetSolutionStepValue(var) = std::max(pressure, 0.0);
    });

    KRATOS_CATCH("")
}

std::string ApplyConstantBoundaryPhreaticLinePressureProcess::Info() const
{
    return "ApplyConstantBoundaryPhreaticLinePressureProcess"s;
}

const ModelPart& ApplyConstantBoundaryPhreaticLinePressureProcess::GetModelPart() const
{
    return mrModelPart;
}

ModelPart& ApplyConstantBoundaryPhreaticLinePressureProcess::GetModelPart() { return mrModelPart; }

const std::string& ApplyConstantBoundaryPhreaticLinePressureProcess::GetVariableName() const
{
    return mVariableName;
}

unsigned int ApplyConstantBoundaryPhreaticLinePressureProcess::GetGravityDirection() const
{
    return mGravityDirection;
}

unsigned int ApplyConstantBoundaryPhreaticLinePressureProcess::GetHorizontalDirection() const
{
    return mHorizontalDirection;
}

double ApplyConstantBoundaryPhreaticLinePressureProcess::GetSpecificWeight() const
{
    return mSpecificWeight;
}

const Vector3& ApplyConstantBoundaryPhreaticLinePressureProcess::GetFirstReferenceCoordinate() const
{
    return mFirstReferenceCoordinate;
}

double ApplyConstantBoundaryPhreaticLinePressureProcess::GetSlope() const { return mSlope; }

double ApplyConstantBoundaryPhreaticLinePressureProcess::GetMinHorizontalCoordinate() const
{
    return mMinHorizontalCoordinate;
}

double ApplyConstantBoundaryPhreaticLinePressureProcess::GetMaxHorizontalCoordinate() const
{
    return mMaxHorizontalCoordinate;
}

} // namespace Kratos