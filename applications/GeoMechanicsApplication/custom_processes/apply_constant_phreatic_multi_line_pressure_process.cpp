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
//                   Jonathan Nuttall

#include "apply_constant_phreatic_multi_line_pressure_process.h"

namespace Kratos
{

ApplyConstantPhreaticMultiLinePressureProcess::ApplyConstantPhreaticMultiLinePressureProcess(ModelPart& model_part,
                                                                                             Parameters rParameters)
        : Process(Flags()) , mrModelPart(model_part)
{
    KRATOS_TRY

    InitializeParameters(rParameters);

    mVariableName = rParameters["variable_name"].GetString();
    mIsFixed = rParameters["is_fixed"].GetBool();
    mIsFixedProvided = rParameters.Has("is_fixed");
    mIsSeepage = rParameters["is_seepage"].GetBool();
    mSpecificWeight = rParameters["specific_weight"].GetDouble();
    mPressureTensionCutOff = rParameters["pressure_tension_cut_off"].GetDouble();
    mOutOfPlaneDirection = rParameters["out_of_plane_direction"].GetInt();

    InitializeCoordinates(rParameters);
    InitializeGravityDirection(rParameters);
    InitializeHorizontalDirection();

    ValidateCoordinates(rParameters);

    KRATOS_CATCH("")
}

void ApplyConstantPhreaticMultiLinePressureProcess::InitializeParameters(Parameters& rParameters) const
{
    Parameters default_parameters(R"(
            {
                "model_part_name":"PLEASE_CHOOSE_MODEL_PART_NAME",
                "variable_name": "PLEASE_PRESCRIBE_VARIABLE_NAME",
                "is_fixed": false,
                "is_seepage": false,
                "gravity_direction": 1,
                "out_of_plane_direction": 2,
                "x_coordinates":           [0.0,1.0],
                "y_coordinates":           [1.0,0.5],
				"z_coordinates":           [0.0,0.0],
                "specific_weight" : 10000.0,
                "pressure_tension_cut_off" : 0.0,
                "table" : [0,1]
            }  )" );

    // Some values are mandatory, since no meaningful default value exist. For this reason try accessing to them
    // So that an error is thrown if they don't exist
    rParameters["x_coordinates"];
    rParameters["y_coordinates"];
    rParameters["z_coordinates"];
    rParameters["variable_name"];
    rParameters["model_part_name"];

    // Now validate the defaults again -- this also ensures no type mismatch
    rParameters.ValidateAndAssignDefaults(default_parameters);
}

void ApplyConstantPhreaticMultiLinePressureProcess::ValidateCoordinates(const Parameters &rParameters) const
{
    if (GravityDirection() == OutOfPlaneDirection())
    {
        KRATOS_ERROR << "Gravity direction cannot be the same as Out-of-Plane directions"
                     << rParameters
                     << std::endl;
    }

    if (!std::is_sorted(HorizontalDirectionCoordinates().begin(), HorizontalDirectionCoordinates().end()))
    {
        KRATOS_ERROR << "The Horizontal Elements Coordinates are not ordered."
                     << rParameters
                     << std::endl;
    }
}

void ApplyConstantPhreaticMultiLinePressureProcess::InitializeCoordinates(const Parameters &rParameters)
{
    mXCoordinates = rParameters["x_coordinates"].GetVector();
    mYCoordinates = rParameters["y_coordinates"].GetVector();
    mZCoordinates = rParameters["z_coordinates"].GetVector();
}

void ApplyConstantPhreaticMultiLinePressureProcess::InitializeGravityDirection(const Parameters &rParameters)
{
    mGravityDirection = rParameters["gravity_direction"].GetInt();
    switch (GravityDirection()) {
        case 0:
            mGravityDirectionCoordinates = XCoordinates();
            break;
        case 1:
            mGravityDirectionCoordinates = YCoordinates();
            break;
        case 2:
            mGravityDirectionCoordinates = ZCoordinates();
            break;
        default:
            KRATOS_ERROR << "The Gravity direction is invalid";
    }
}

void ApplyConstantPhreaticMultiLinePressureProcess::InitializeHorizontalDirection()
{
    mHorizontalDirection = 0;
    for (unsigned int i=0; i<N_DIM_3D; ++i)
        if (i != GravityDirection() && i != OutOfPlaneDirection()) mHorizontalDirection = i;

    switch (HorizontalDirection()) {
        case 0:
            mHorizontalDirectionCoordinates = XCoordinates();
            break;
        case 1:
            mHorizontalDirectionCoordinates = YCoordinates();
            break;
        case 2:
            mHorizontalDirectionCoordinates = ZCoordinates();
            break;
        default:
            KRATOS_ERROR << "The Horizontal direction is invalid";
    }
}

void ApplyConstantPhreaticMultiLinePressureProcess::ExecuteInitialize()
{
    KRATOS_TRY

    KRATOS_ERROR_IF(mrModelPart.NumberOfNodes() <= 0) << "Number of nodes is smaller than or equal to zero";

    const Variable<double> &var = KratosComponents< Variable<double> >::Get(VariableName());

    block_for_each(mrModelPart.Nodes(), [&var, this](Node& rNode) {
        const double pressure = CalculatePressure(rNode);
        if (IsSeepage()) {
            if (pressure < PORE_PRESSURE_SIGN_FACTOR * PressureTensionCutOff()) {
                rNode.FastGetSolutionStepValue(var) = pressure;
                if (IsFixed()) rNode.Fix(var);
            } else {
                if (IsFixedProvided()) rNode.Free(var);
            }
        } else {
            if (IsFixed()) rNode.Fix(var);
            else if (IsFixedProvided()) rNode.Free(var);
            rNode.FastGetSolutionStepValue(var) = std::min(pressure, PORE_PRESSURE_SIGN_FACTOR * PressureTensionCutOff());
        }
    });

    KRATOS_CATCH("")
}

std::string ApplyConstantPhreaticMultiLinePressureProcess::Info() const
{
    return "ApplyConstantPhreaticMultiLinePressureProcess";
}

void ApplyConstantPhreaticMultiLinePressureProcess::PrintInfo(std::ostream &rOStream) const
{
    rOStream << "ApplyConstantPhreaticMultiLinePressureProcess";
}

const std::string &ApplyConstantPhreaticMultiLinePressureProcess::VariableName() const
{
    return mVariableName;
}

bool ApplyConstantPhreaticMultiLinePressureProcess::IsFixed() const
{
    return mIsFixed;
}

bool ApplyConstantPhreaticMultiLinePressureProcess::IsFixedProvided() const
{
    return mIsFixedProvided;
}

bool ApplyConstantPhreaticMultiLinePressureProcess::IsSeepage() const
{
    return mIsSeepage;
}

unsigned int ApplyConstantPhreaticMultiLinePressureProcess::GravityDirection() const
{
    return mGravityDirection;
}

double ApplyConstantPhreaticMultiLinePressureProcess::SpecificWeight() const
{
    return mSpecificWeight;
}

unsigned int ApplyConstantPhreaticMultiLinePressureProcess::OutOfPlaneDirection() const
{
    return mOutOfPlaneDirection;
}

unsigned int ApplyConstantPhreaticMultiLinePressureProcess::HorizontalDirection() const
{
    return mHorizontalDirection;
}

const Vector &ApplyConstantPhreaticMultiLinePressureProcess::HorizontalDirectionCoordinates() const
{
    return mHorizontalDirectionCoordinates;
}

const Vector &ApplyConstantPhreaticMultiLinePressureProcess::GravityDirectionCoordinates() const
{
    return mGravityDirectionCoordinates;
}

const Vector &ApplyConstantPhreaticMultiLinePressureProcess::XCoordinates() const
{
    return mXCoordinates;
}

const Vector &ApplyConstantPhreaticMultiLinePressureProcess::YCoordinates() const
{
    return mYCoordinates;
}

const Vector &ApplyConstantPhreaticMultiLinePressureProcess::ZCoordinates() const
{
    return mZCoordinates;
}

double ApplyConstantPhreaticMultiLinePressureProcess::PressureTensionCutOff() const
{
    return mPressureTensionCutOff;
}

int ApplyConstantPhreaticMultiLinePressureProcess::findIndex(const Node &rNode) const
{
    const auto& coords = rNode.Coordinates();
    const auto number_of_coordinates = static_cast<int>(HorizontalDirectionCoordinates().size());
    for (int index = 0; index < number_of_coordinates; ++index)
    {
        if (HorizontalDirectionCoordinates()[index] >= coords[HorizontalDirection()])
        {
            return index == 0 ? index : index - 1;
        }
    }

    return number_of_coordinates - 1;
}

double ApplyConstantPhreaticMultiLinePressureProcess::CalculatePressure(const Node &rNode,
                                                                        std::vector<double> deltaH) const
{
    // find nodes in horizontalDirectionCoordinates
    const int firstPointIndex = findIndex(rNode);
    const int secondPointIndex = firstPointIndex + 1;

    array_1d<double, 2> y;
    y[0] = GravityDirectionCoordinates()[firstPointIndex];
    y[1] = GravityDirectionCoordinates()[secondPointIndex];

    if (!deltaH.empty())
    {
        y[0] += deltaH[firstPointIndex];
        y[1] += deltaH[secondPointIndex];
    }


    const double slope = (y[1] - y[0])
                         / (HorizontalDirectionCoordinates()[secondPointIndex] - HorizontalDirectionCoordinates()[firstPointIndex]);

    const double height = slope * (rNode.Coordinates()[HorizontalDirection()] - HorizontalDirectionCoordinates()[firstPointIndex]) + y[0];
    const double distance = height - rNode.Coordinates()[GravityDirection()];
    return - PORE_PRESSURE_SIGN_FACTOR * SpecificWeight() * distance;
}

}