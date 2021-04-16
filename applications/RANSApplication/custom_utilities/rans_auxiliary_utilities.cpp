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
#include <vector>

// External includes

// Project includes
#include "input_output/logger.h"
#include "includes/define.h"
#include "includes/ublas_interface.h"
#include "utilities/parallel_utilities.h"

// Application includes

// Include base h
#include "rans_auxiliary_utilities.h"

namespace Kratos
{
void RansAuxiliaryUtilities::CalculateFrequencyBinValues(
    Vector& rRealValues,
    Vector& rImaginaryValues,
    const double Value,
    const IndexType CurrentStep,
    const IndexType TotalSteps,
    const std::vector<int>& rFrequencyBinIndexList)
{
    KRATOS_TRY

    KRATOS_WARNING_IF("CalculateFrequencyBinValues", rFrequencyBinIndexList.size() >= TotalSteps / 2)
        << "Frequency bin list is larger than allowed maximum number of frequencies according to Nyquist theorem. [ FrequencyBinIndexList.size() = "
        << rFrequencyBinIndexList.size() << ", allowed maximum frequency bin indices count = " << TotalSteps / 2 << " ].\n";

    KRATOS_WARNING_IF("CalculateFrequencyBinValues", CurrentStep > TotalSteps)
        << "Current time step index is larger than the total number of steps. [ Current times step index = "
        << CurrentStep << ", total number of steps = " << TotalSteps << " ].\n";

    if (rRealValues.size() != rFrequencyBinIndexList.size()) {
        rRealValues.resize(rFrequencyBinIndexList.size());
    }

    if (rImaginaryValues.size() != rFrequencyBinIndexList.size()) {
        rImaginaryValues.resize(rFrequencyBinIndexList.size());
    }

    block_for_each(rFrequencyBinIndexList, [&](const int FrequencyBinIndex){
        rRealValues[FrequencyBinIndex] =
            Value * std::cos(2.0 * M_PI * FrequencyBinIndex * CurrentStep / TotalSteps);
        rImaginaryValues[FrequencyBinIndex] =
            Value * std::sin(2.0 * M_PI * FrequencyBinIndex * CurrentStep / TotalSteps);
    });

    KRATOS_CATCH("");
}

void RansAuxiliaryUtilities::ApplyHannWindow(
    Vector& rOutput,
    const Vector& rInput,
    const IndexType WindowStartIndex,
    const IndexType WindowEndIndex)
{
    KRATOS_TRY

    KRATOS_ERROR_IF(WindowStartIndex >= rInput.size())
        << "Provided startining index is out of bounds. [ WindowStartIndex = " << WindowStartIndex
        << ", Input.size() = " << rInput.size() << " ].\n";
    KRATOS_ERROR_IF(WindowEndIndex >= rInput.size())
        << "Provided ending index is out of bounds. [ WindowEndIndex = " << WindowEndIndex
        << ", Input.size() = " << rInput.size() << " ].\n";
    KRATOS_ERROR_IF(WindowEndIndex - WindowStartIndex <= 0)
        << "Provided start/end index is not in correct order. [ "
           "WindowStartIndex = "
        << WindowStartIndex << ", WindowEndIndex = " << WindowEndIndex << " ].\n";

    if (rOutput.size() != rInput.size()) {
        rOutput.resize(rInput.size());
    }

    noalias(rOutput) = ZeroVector(rInput.size());

    const IndexType number_of_values = WindowEndIndex - WindowStartIndex + 1;

    KRATOS_ERROR_IF(number_of_values > rOutput.size())
        << "Number of values calculated from provided window start/end indices "
           "out of bounds in the given input vector. [ Number of values = "
        << number_of_values << ", rInput.size() = " << rInput.size() << " ].\n";

    IndexPartition<int>(number_of_values).for_each([&](const int Index){
        const IndexType current_offsetted_index = Index + WindowStartIndex;
        rOutput[current_offsetted_index] = rInput[current_offsetted_index] * 0.5 * (1 - std::cos(2.0 * M_PI * Index / number_of_values));
    });

    KRATOS_CATCH("");
}

double RansAuxiliaryUtilities::VectorSummation(
    const Vector& rValues)
{
    KRATOS_TRY

    return IndexPartition<int>(rValues.size()).for_each<SumReduction<double>>([&](const int Index) {
        return rValues[Index];
    });

    KRATOS_CATCH("");
}

} // namespace Kratos.
