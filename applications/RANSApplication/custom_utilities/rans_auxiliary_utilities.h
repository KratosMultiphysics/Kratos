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

#if !defined(KRATOS_RANS_AUXILIARY_UTILITIES_H_INCLUDED)
#define KRATOS_RANS_AUXILIARY_UTILITIES_H_INCLUDED

// System includes
#include <vector>
#include <string>

// External includes

// Project includes
#include "includes/ublas_interface.h"

// Application includes

namespace Kratos
{
///@name Kratos Classes
///@{

class KRATOS_API(RANS_APPLICATION) RansAuxiliaryUtilities
{
public:
    ///@name Type Definitions

    using IndexType = std::size_t;

    ///@}
    ///@name Static Operations
    ///@{

    static void CalculateFrequencyBinValues(
        Vector& rRealValues,
        Vector& rImaginaryValues,
        const double Value,
        const IndexType CurrentStep,
        const IndexType TotalSteps,
        const std::vector<int>& rFrequencyBinIndexList);

    static void ApplyHannWindow(
        Vector& rOutput,
        const Vector& rInput,
        const IndexType WindowStartIndex,
        const IndexType WindowEndIndex);

    static double VectorSummation(
        const Vector& rValues);

    ///@}

}; // Class RansAuxiliaryUtilities

///@}

} // namespace Kratos.

#endif // KRATOS_RANS_AUXILIARY_UTILITIES_H_INCLUDED defined
