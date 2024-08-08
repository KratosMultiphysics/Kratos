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

#if !defined(KRATOS_FLUID_FFT_UTILITIES_H)
#define KRATOS_FLUID_FFT_UTILITIES_H

// System includes
#include <vector>
#include <tuple>

// External includes

// Project includes
#include "includes/kratos_parameters.h"
#include "includes/process_info.h"

// Application includes

namespace Kratos
{
///@addtogroup FluidDynamicsApplication
///@{

///@name Kratos classes
///@{

class KRATOS_API(FLUID_DYNAMICS_APPLICATION) FluidFFTUtilities
{
public:
    ///@name Type Definitions
    ///@{

    using IndexType = std::size_t;

    using VectorType = std::vector<double>;

    ///@}
    ///@name Life Cycle
    ///@{

    FluidFFTUtilities(
        const double TotalTime,
        const double WindowingLength,
        const double DeltaTime);

    ///@}
    ///@name Operations
    ///@{

    std::tuple<VectorType, VectorType, VectorType, VectorType> CalculateFFTFrequencyDistribution(
        const VectorType& rValues) const;

    bool IsWithinWindowingRange(const double Time) const;

    double CalculateHannWindowCoefficient(const double Time) const;

    double CalculateFFTRealCoefficient(
        const IndexType FrequencyBinIndex,
        const double Time) const;

    double CalculateFFTImagCoefficient(
        const IndexType FrequencyBinIndex,
        const double Time) const;

    double CalculateFFTAmplitudeSquare(
        const double RealValue,
        const double ImagValue) const;

    double CalculateFFTAmplitudeSquareDerivative(
        const double RealValue,
        const double RealValueDerivative,
        const double ImagValue,
        const double ImagValueDerivative) const;

    double GetFrequencyResolution() const;

    double GetFrequency(const double FrequencyBinIndex) const;

    double GetMaximumFrequency() const;

    IndexType GetTotalNumberOfSteps() const;

    IndexType GetNumberOfWindowingSteps() const;

    ///@}
private:
    ///@name Private Members
    ///@{

    const double mTotalTime;
    const double mWindowingLength;
    const double mDeltaTime;

    double mFrequencyResolution;
    double mWindowingStepsCoefficient;
    double mWindowingStepsDerivativeCoefficient;
    IndexType mTotalSteps;
    IndexType mWindowingSteps;
    IndexType mNumberOfFrequencies;

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
};

///@}

///@}

} // namespace Kratos

#endif // KRATOS_FLUID_FFT_UTILITIES_H
