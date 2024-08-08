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
#include <cmath>

// External includes

// Project includes
#include "includes/variables.h"
#include "utilities/parallel_utilities.h"

// Application includes

// Include h
#include "fluid_fft_utilities.h"


namespace Kratos
{

FluidFFTUtilities::FluidFFTUtilities(
    const double TotalTime,
    const double WindowingLength,
    const double DeltaTime)
    : mTotalTime(TotalTime),
      mWindowingLength(WindowingLength),
      mDeltaTime(DeltaTime)
{
    KRATOS_TRY

    KRATOS_ERROR_IF(mWindowingLength <= 0.0)
        << "Invalid windowing length provided. Windowing length should be positive. [ windowing_length = "
        << mWindowingLength << " ].\n";

    KRATOS_ERROR_IF(mTotalTime <= 0.0)
        << "Invalid total time provided. Total time should be positive. [ total_time = "
        << mTotalTime << " ].\n";

    KRATOS_ERROR_IF(mDeltaTime <= 0.0)
        << "Invalid delta time provided. Delta time should be positive. [ delta_time = "
        << mDeltaTime << " ].\n";

    KRATOS_ERROR_IF(mTotalTime -  mWindowingLength < 0)
        << "Total time should be equal or grater than windowing time. Try reducing windowing length or increasing total time. [ total_time = "
        << mTotalTime << ", windowing_length = " << mWindowingLength << " ].\n";

    // 0.5 is added to avoid rounding off errors
    mTotalSteps = static_cast<IndexType>(mTotalTime / mDeltaTime + 0.5);
    mWindowingSteps = static_cast<IndexType>(mWindowingLength / mDeltaTime + 0.5);
    mNumberOfFrequencies = static_cast<IndexType>(mTotalSteps / 2);
    mFrequencyResolution = 1.0 / (mDeltaTime * mTotalSteps);
    mWindowingStepsCoefficient = 16.0 / std::pow(mWindowingSteps, 2);
    mWindowingStepsDerivativeCoefficient = 32.0 / std::pow(mWindowingSteps, 2);

    KRATOS_CATCH("")
}

std::tuple<std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>> FluidFFTUtilities::CalculateFFTFrequencyDistribution(
    const VectorType& rValues) const
{
    KRATOS_TRY

    KRATOS_ERROR_IF(rValues.size() != mTotalSteps)
        << "Number of steps in values does not match with number of steps from total time. [ number of steps in values = "
        << rValues.size() << ", number of steps from total time = " << mTotalSteps << ", total_time = "
        << mTotalTime << ", delta_time = " << mDeltaTime << " ].\n";

    VectorType frequency_real_components(mNumberOfFrequencies, 0.0);
    VectorType frequency_imag_components(mNumberOfFrequencies, 0.0);

    VectorType hann_windowing_coefficients(mWindowingSteps);
    IndexPartition<IndexType>(mWindowingSteps).for_each([&](const IndexType WindowingStepIndex){
        const double current_time = (mTotalSteps - mWindowingSteps + WindowingStepIndex + 1) * mDeltaTime;
        hann_windowing_coefficients[WindowingStepIndex] = CalculateHannWindowCoefficient(current_time);
    });

    IndexPartition<IndexType>(mNumberOfFrequencies).for_each([&](const IndexType& FrequencyBinIndex) {
        for (IndexType n = mTotalSteps - mWindowingSteps + 1; n <= mTotalSteps; ++n){
            const double drag = rValues[n-1];
            const double window_coefficient = hann_windowing_coefficients[n - (mTotalSteps - mWindowingSteps + 1)];

            const double current_time = n * mDeltaTime;
            const double real_coefficient = CalculateFFTRealCoefficient(FrequencyBinIndex, current_time);
            const double imag_coefficient = CalculateFFTImagCoefficient(FrequencyBinIndex, current_time);
            frequency_real_components[FrequencyBinIndex] += window_coefficient * drag * real_coefficient;
            frequency_imag_components[FrequencyBinIndex] += window_coefficient * drag * imag_coefficient;
        }
    });

    VectorType frequency_amplitudes_square(mNumberOfFrequencies);
    VectorType frequency_list(mNumberOfFrequencies);

    IndexPartition<IndexType>(mNumberOfFrequencies).for_each([&](const IndexType& FrequencyBinIndex) {
        frequency_list[FrequencyBinIndex] = FrequencyBinIndex * mFrequencyResolution;
        frequency_amplitudes_square[FrequencyBinIndex] = CalculateFFTAmplitudeSquare(frequency_real_components[FrequencyBinIndex], frequency_imag_components[FrequencyBinIndex]);
    });

    return std::make_tuple<VectorType, VectorType, VectorType, VectorType>(
        std::forward<VectorType>(frequency_list),
        std::forward<VectorType>(frequency_real_components),
        std::forward<VectorType>(frequency_imag_components),
        std::forward<VectorType>(frequency_amplitudes_square));

    KRATOS_CATCH("");
}

bool FluidFFTUtilities::IsWithinWindowingRange(const double Time) const
{
    return Time - mTotalTime + mWindowingLength >= -0.5 * mDeltaTime;
}

double FluidFFTUtilities::CalculateHannWindowCoefficient(const double Time) const
{
    KRATOS_TRY

    const double offsetted_time = std::max(std::min(Time, mTotalTime) - mTotalTime + mWindowingLength, 0.0);
    return 0.5 * (1.0 - std::cos(2.0 * M_PI * offsetted_time / mWindowingLength));

    KRATOS_CATCH("");
}

double FluidFFTUtilities::CalculateFFTRealCoefficient(
    const IndexType FrequencyBinIndex,
    const double Time) const
{
    KRATOS_TRY

    KRATOS_ERROR_IF(Time < 0.0) << "Time should be greater than or equal to zero. [ Time = " << Time << " ].\n";

    return std::cos(2.0 * M_PI * FrequencyBinIndex * Time / mTotalTime);

    KRATOS_CATCH("");
}

double FluidFFTUtilities::CalculateFFTImagCoefficient(
    const IndexType FrequencyBinIndex,
    const double Time) const
{
    KRATOS_TRY

    KRATOS_ERROR_IF(Time < 0.0) << "Time should be greater than or equal to zero. [ Time = " << Time << " ].\n";

    return std::sin(2.0 * M_PI * FrequencyBinIndex * Time / mTotalTime);

    KRATOS_CATCH("");
}

double FluidFFTUtilities::CalculateFFTAmplitudeSquare(
    const double RealValue,
    const double ImagValue) const
{
    KRATOS_TRY

    return (std::pow(RealValue, 2) + std::pow(ImagValue, 2)) * mWindowingStepsCoefficient;

    KRATOS_CATCH("");
}

double FluidFFTUtilities::CalculateFFTAmplitudeSquareDerivative(
    const double RealValue,
    const double RealValueDerivative,
    const double ImagValue,
    const double ImagValueDerivative) const
{
    KRATOS_TRY

    return (RealValue * RealValueDerivative + ImagValue * ImagValueDerivative) * mWindowingStepsDerivativeCoefficient;

    KRATOS_CATCH("");
}

double FluidFFTUtilities::GetFrequencyResolution() const
{
    KRATOS_TRY

    return 1.0 / mTotalTime;

    KRATOS_CATCH("");
}

double FluidFFTUtilities::GetFrequency(const double FrequencyBinIndex) const
{
    KRATOS_TRY

    return FrequencyBinIndex / mTotalTime;

    KRATOS_CATCH("");
}

double FluidFFTUtilities::GetMaximumFrequency() const
{
    KRATOS_TRY

    return 1.0 / (mDeltaTime * 2.0);

    KRATOS_CATCH("");
}

std::size_t FluidFFTUtilities::GetTotalNumberOfSteps() const
{
    KRATOS_TRY

    return mTotalSteps;

    KRATOS_CATCH("");
}

std::size_t FluidFFTUtilities::GetNumberOfWindowingSteps() const
{
    KRATOS_TRY

    return mWindowingSteps;

    KRATOS_CATCH("");
}

} // namespace Kratos
