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

#include "apply_constant_phreatic_line_pressure_process.h"
#include "includes/model_part.h"

#include <algorithm>

namespace Kratos
{

ApplyConstantPhreaticLinePressureProcess::ApplyConstantPhreaticLinePressureProcess(ModelPart& model_part,
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
                "is_seepage": false,
                "gravity_direction": 1,
                "out_of_plane_direction": 2,
                "first_reference_coordinate":           [0.0,1.0,0.0],
                "second_reference_coordinate":          [1.0,0.5,0.0],
                "specific_weight" : 10000.0,
                "pressure_tension_cut_off" : 0.0,
                "table" : [0,1]
            }  )");

    // Some values need to be mandatorily prescribed since no meaningful default value exist.
    // For this reason try accessing to them So that an error is thrown if they don't exist
    rParameters["first_reference_coordinate"];
    rParameters["second_reference_coordinate"];
    rParameters["variable_name"];
    rParameters["model_part_name"];

    mIsFixedProvided = rParameters.Has("is_fixed");

    // Now validate agains defaults -- this also ensures no type mismatch
    rParameters.ValidateAndAssignDefaults(default_parameters);

    mVariableName        = rParameters["variable_name"].GetString();
    mIsFixed             = rParameters["is_fixed"].GetBool();
    mIsSeepage           = rParameters["is_seepage"].GetBool();
    mGravityDirection    = rParameters["gravity_direction"].GetInt();
    mOutOfPlaneDirection = rParameters["out_of_plane_direction"].GetInt();
    if (mGravityDirection == mOutOfPlaneDirection)
        KRATOS_ERROR << "Gravity direction cannot be the same as Out-of-Plane directions"
                     << rParameters << std::endl;

    mHorizontalDirection = 0;
    for (unsigned int i = 0; i < N_DIM_3D; ++i)
        if (i != mGravityDirection && i != mOutOfPlaneDirection) mHorizontalDirection = i;

    mFirstReferenceCoordinate  = rParameters["first_reference_coordinate"].GetVector();
    mSecondReferenceCoordinate = rParameters["second_reference_coordinate"].GetVector();

    mMinHorizontalCoordinate = std::min(mFirstReferenceCoordinate[mHorizontalDirection],
                                        mSecondReferenceCoordinate[mHorizontalDirection]);
    mMaxHorizontalCoordinate = std::max(mFirstReferenceCoordinate[mHorizontalDirection],
                                        mSecondReferenceCoordinate[mHorizontalDirection]);

    if (mMaxHorizontalCoordinate <= mMinHorizontalCoordinate) {
        KRATOS_ERROR
            << "First and second point on the phreatic line have the same horizontal coordinate"
            << rParameters << std::endl;
    }

    mSlope =
        (mSecondReferenceCoordinate[mGravityDirection] - mFirstReferenceCoordinate[mGravityDirection]) /
        (mSecondReferenceCoordinate[mHorizontalDirection] - mFirstReferenceCoordinate[mHorizontalDirection]);

    mSpecificWeight        = rParameters["specific_weight"].GetDouble();
    mPressureTensionCutOff = rParameters["pressure_tension_cut_off"].GetDouble();

    KRATOS_CATCH("")
}

void ApplyConstantPhreaticLinePressureProcess::ExecuteInitialize()
{
    KRATOS_TRY

    const Variable<double>& var = KratosComponents<Variable<double>>::Get(mVariableName);
    block_for_each(mrModelPart.Nodes(), [&var, this](Node& rNode) {
        const double pressure = PORE_PRESSURE_SIGN_FACTOR * CalculatePressure(rNode);

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
            if (mIsFixed) rNode.Fix(var);
            else if (mIsFixedProvided) rNode.Free(var);
        }
    });

    KRATOS_CATCH("")
}

std::string ApplyConstantPhreaticLinePressureProcess::Info() const
{
    return "ApplyConstantPhreaticLinePressureProcess";
}

double ApplyConstantPhreaticLinePressureProcess::CalculatePressure(const Node& rNode) const
{
    double x            = rNode.Coordinates()[mHorizontalDirection];
    x                   = std::max(x, mMinHorizontalCoordinate);
    x                   = std::min(x, mMaxHorizontalCoordinate);
    const double height = mSlope * (x - mFirstReferenceCoordinate[mHorizontalDirection]) +
                          mFirstReferenceCoordinate[mGravityDirection];
    const double distance = height - rNode.Coordinates()[mGravityDirection];
    return -mSpecificWeight * distance;
}

} // namespace Kratos