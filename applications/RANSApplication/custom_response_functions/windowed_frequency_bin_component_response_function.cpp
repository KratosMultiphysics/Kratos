//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

// System includes
#include <string>
#include <cmath>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "includes/ublas_interface.h"
#include "response_functions/adjoint_response_function.h"
#include "utilities/brute_force_point_locator.h"

// Application includes
#include "custom_utilities/fluid_calculation_utilities.h"

// Include base h
#include "windowed_frequency_bin_component_response_function.h"

namespace Kratos
{
template<unsigned int TDim>
WindowedFrequencyBinComponentResponseFunction<TDim>::WindowedFrequencyBinComponentResponseFunction(
    Parameters Settings,
    ModelPart& rModelPart)
    : mrModelPart(rModelPart)
{
    KRATOS_TRY;

    Parameters default_settings(R"(
        {
            "point_coordinates"         : [0.0, 0.0, 0.0],
            "velocity_direction"        : [0.0, 1.0, 0.0],
            "is_real_component"         : true,
            "frequency_bin_index"       : 1,
            "window_size"               : 300,
            "total_number_of_time_steps": 1000,
            "echo_level"                : 0
        })");

    Settings.ValidateAndAssignDefaults(default_settings);

    mEchoLevel = Settings["echo_level"].GetInt();

    KRATOS_ERROR_IF(Settings["point_coordinates"].GetVector().size() != 3)
        << "point coordinates should be of size 3. [ "
           "point_coordinates.size() = "
        << Settings["point_coordinates"].GetVector().size() << " ].\n";
    noalias(mPointCoordinates) = Settings["point_coordinates"].GetVector();

    KRATOS_ERROR_IF(Settings["velocity_direction"].GetVector().size() != 3)
        << "point coordinates should be of size 3. [ "
           "velocity_direction.size() = "
        << Settings["velocity_direction"].GetVector().size() << " ].\n";
    noalias(mVelocityDirection) = Settings["velocity_direction"].GetVector();

    mFrequencyBinIndex = Settings["frequency_bin_index"].GetInt();
    KRATOS_ERROR_IF(mFrequencyBinIndex < 0)
        << "Frequency bin index should be greater than or equal to zero. [ "
           "frequency_bin_index = "
        << mFrequencyBinIndex << " ].\n";

    mTotalNumberOfTimeSteps = Settings["total_number_of_time_steps"].GetInt();
    KRATOS_ERROR_IF(mTotalNumberOfTimeSteps <= 0)
        << "Total number of time steps should be greater than to zero. [ "
           "total_number_of_time_steps = "
        << mTotalNumberOfTimeSteps << " ].\n";

    mIsRealComponentRequested = Settings["is_real_component"].GetBool();

    mWindowSize = Settings["window_size"].GetInt();
    KRATOS_ERROR_IF(mWindowSize <= 0)
        << "Window size should be greater than to zero. [ "
        << "window_size = " << mWindowSize << " ].\n";

    KRATOS_CATCH("");
}

template <unsigned int TDim>
void WindowedFrequencyBinComponentResponseFunction<TDim>::Initialize()
{
    KRATOS_TRY;

    FindPointElement();

    KRATOS_ERROR_IF(mPointElementId == 0)
        << "No element or node can be found for the point at "
        << mPointCoordinates << " in " << mrModelPart.Name() << ".\n";

    const double delta_time = mrModelPart.GetProcessInfo()[DELTA_TIME];

    KRATOS_INFO_IF("WindowedFrequencyBinComponentResponseFunction", mEchoLevel > 0) << "Summary: \n"
        << "\tPoint coordinates         : " << mPointCoordinates << "\n"
        << "\tElement id for given point: " << mPointElementId << "\n"
        << "\tVelocity direction        : " << mVelocityDirection << "\n"
        << "\tWindow size               : " << mWindowSize << "\n"
        << "\tWindow duration           : " << mWindowSize * delta_time << " s\n"
        << "\tEvaluated frequency bin   : " << mFrequencyBinIndex << "\n"
        << "\tEvaluated frequency       : " << mFrequencyBinIndex / (delta_time * mTotalNumberOfTimeSteps) << " Hz \n"
        << "\tMaximum possible frequency: " << 1.0 / (delta_time * 2.0) << " Hz \n"
        << "\tTotal number of time steps: " << mTotalNumberOfTimeSteps << "\n"
        << "\tComponent type            : " << (mIsRealComponentRequested ? "real\n" : "imaginary\n");

    KRATOS_CATCH("");
}

template <unsigned int TDim>
void WindowedFrequencyBinComponentResponseFunction<TDim>::CalculateGradient(
    const Element& rAdjointElement,
    const Matrix& rResidualGradient,
    Vector& rResponseGradient,
    const ProcessInfo& rProcessInfo)
{
    rResponseGradient.clear();
}

template <unsigned int TDim>
void WindowedFrequencyBinComponentResponseFunction<TDim>::CalculateGradient(
    const Condition& rAdjointCondition,
    const Matrix& rResidualGradient,
    Vector& rResponseGradient,
    const ProcessInfo& rProcessInfo)
{
    rResponseGradient.clear();
}

template <unsigned int TDim>
void WindowedFrequencyBinComponentResponseFunction<TDim>::CalculateFirstDerivativesGradient(
    const Element& rAdjointElement,
    const Matrix& rResidualGradient,
    Vector& rResponseGradient,
    const ProcessInfo& rProcessInfo)
{
    if (rResponseGradient.size() != rResidualGradient.size1()) {
        rResponseGradient.resize(rResidualGradient.size1());
    }

    const int step = mrModelPart.GetProcessInfo()[STEP];
    const int offsetted_step = step - mTotalNumberOfTimeSteps + mWindowSize;

    if (rAdjointElement.Id() == mPointElementId && offsetted_step >= 0) {
        const auto& r_geometry = rAdjointElement.GetGeometry();
        const double hann_window_value = CalculateHannWindowValue();

        IndexType local_index = 0;
        if (mIsRealComponentRequested) {
            for (IndexType c = 0; c < r_geometry.PointsNumber(); ++c) {
                for (IndexType k = 0; k < TDim; ++k) {
                    rResponseGradient[local_index++] =
                        hann_window_value * mPointShapeFunctionValues[c] * mVelocityDirection[k] *
                        std::cos(2.0 * M_PI * mFrequencyBinIndex * step / mTotalNumberOfTimeSteps);
                }
            }
        } else {
            for (IndexType c = 0; c < r_geometry.PointsNumber(); ++c) {
                for (IndexType k = 0; k < TDim; ++k) {
                    rResponseGradient[local_index++] =
                        hann_window_value * mPointShapeFunctionValues[c] * mVelocityDirection[k] *
                        std::sin(2.0 * M_PI * mFrequencyBinIndex * step / mTotalNumberOfTimeSteps);
                }
            }
        }
    } else {
        rResponseGradient.clear();
    }
}

template <unsigned int TDim>
void WindowedFrequencyBinComponentResponseFunction<TDim>::CalculateFirstDerivativesGradient(
    const Condition& rAdjointCondition,
    const Matrix& rResidualGradient,
    Vector& rResponseGradient,
    const ProcessInfo& rProcessInfo)
{
    rResponseGradient.clear();
}

template <unsigned int TDim>
void WindowedFrequencyBinComponentResponseFunction<TDim>::CalculateSecondDerivativesGradient(
    const Element& rAdjointElement,
    const Matrix& rResidualGradient,
    Vector& rResponseGradient,
    const ProcessInfo& rProcessInfo)
{
    if (rResponseGradient.size() != rResidualGradient.size1()) {
        rResponseGradient.resize(rResidualGradient.size1());
    }

    rResponseGradient.clear();
}

template <unsigned int TDim>
void WindowedFrequencyBinComponentResponseFunction<TDim>::CalculateSecondDerivativesGradient(
    const Condition& rAdjointCondition,
    const Matrix& rResidualGradient,
    Vector& rResponseGradient,
    const ProcessInfo& rProcessInfo)
{
    rResponseGradient.clear();
}

template <unsigned int TDim>
void WindowedFrequencyBinComponentResponseFunction<TDim>::CalculatePartialSensitivity(
    Element& rAdjointElement,
    const Variable<array_1d<double, 3>>& rVariable,
    const Matrix& rSensitivityMatrix,
    Vector& rSensitivityGradient,
    const ProcessInfo& rProcessInfo)
{
    if (rSensitivityGradient.size() != rSensitivityMatrix.size1())
        rSensitivityGradient.resize(rSensitivityMatrix.size1(), false);

    rSensitivityGradient.clear();
}

template <unsigned int TDim>
void WindowedFrequencyBinComponentResponseFunction<TDim>::CalculatePartialSensitivity(
    Condition& rAdjointCondition,
    const Variable<array_1d<double, 3>>& rVariable,
    const Matrix& rSensitivityMatrix,
    Vector& rSensitivityGradient,
    const ProcessInfo& rProcessInfo)
{
    if (rSensitivityGradient.size() != rSensitivityMatrix.size1())
        rSensitivityGradient.resize(rSensitivityMatrix.size1(), false);

    rSensitivityGradient.clear();
}

template <unsigned int TDim>
double WindowedFrequencyBinComponentResponseFunction<TDim>::CalculateValue(ModelPart& rModelPart)
{
    KRATOS_TRY

    const auto& element = mrModelPart.GetElement(mPointElementId);

    array_1d<double, 3> current_velocity;
    FluidCalculationUtilities::EvaluateInPoint(element.GetGeometry(), mPointShapeFunctionValues,
                                               std::tie(current_velocity, VELOCITY));

    const double directional_velocity = inner_prod(current_velocity, mVelocityDirection);

    const int step = mrModelPart.GetProcessInfo()[STEP];

    KRATOS_WARNING_IF("WindowedFrequencyBinComponentResponseFunction", step > mTotalNumberOfTimeSteps)
        << "Current time step is greater than the total number of time steps. [ current time step = "
        << step << ", total number of time steps = " << mTotalNumberOfTimeSteps << " ].\n";

    const int offsetted_step = step - mTotalNumberOfTimeSteps + mWindowSize;

    if (offsetted_step >= 0) {
        const double hann_window_value = CalculateHannWindowValue();
        if (mIsRealComponentRequested) {
            return directional_velocity * hann_window_value *
                   std::cos(2.0 * M_PI * mFrequencyBinIndex * step / mTotalNumberOfTimeSteps);
        } else {
            return directional_velocity * hann_window_value *
                   std::sin(2.0 * M_PI * mFrequencyBinIndex * step / mTotalNumberOfTimeSteps);
        }
    } else {
        return 0.0;
    }

    KRATOS_CATCH("");
}

template <unsigned int TDim>
void WindowedFrequencyBinComponentResponseFunction<TDim>::FindPointElement()
{
    KRATOS_TRY

    BruteForcePointLocator brute_force_point_locator(mrModelPart);
    Point current_point(mPointCoordinates[0], mPointCoordinates[1], mPointCoordinates[2]);

    mPointElementId = brute_force_point_locator.FindElement(
        current_point, mPointShapeFunctionValues);

    KRATOS_CATCH("");
}

template <unsigned int TDim>
double WindowedFrequencyBinComponentResponseFunction<TDim>::CalculateHannWindowValue()
{
    KRATOS_TRY

    const int step = mrModelPart.GetProcessInfo()[STEP];
    const int offsetted_step = step - mTotalNumberOfTimeSteps + mWindowSize;
    return 0.5 * (1.0 - std::cos(2.0 * M_PI * offsetted_step / mWindowSize));

    KRATOS_CATCH("");
}

// template instantiations

template class WindowedFrequencyBinComponentResponseFunction<2>;
template class WindowedFrequencyBinComponentResponseFunction<3>;

} /* namespace Kratos.*/


