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

#pragma once

// System includes
#include <vector>
#include <tuple>
#include <string>

// External includes

// Project includes

// Application includes

namespace Kratos
{
///@addtogroup OptimizationApplication
///@{

///@name Kratos classes
///@{

class KRATOS_API(FLUID_DYNAMICS_APPLICATION) FFTUtils
{
public:
    ///@name Type Definitions
    ///@{

    using IndexType = std::size_t;

    using VectorType = std::vector<double>;

    KRATOS_CLASS_POINTER_DEFINITION(FFTUtils);

    ///@}
    ///@name Life Cycle
    ///@{

    FFTUtils(
        const double TotalTime,
        const double WindowingLength,
        const double DeltaTime);

    ///@}
    ///@name Public operations
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

    std::string Info() const;

    ///@}
private:
    ///@name Private member variables
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
};

///@}

///@}

} // namespace Kratos

