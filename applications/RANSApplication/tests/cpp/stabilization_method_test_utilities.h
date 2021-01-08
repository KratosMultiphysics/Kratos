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

#if !defined(KRATOS_STABILIZATION_METHOD_TEST_UTILITIES_H_INCLUDED)
#define KRATOS_STABILIZATION_METHOD_TEST_UTILITIES_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/process_info.h"

// Application includes

namespace Kratos
{
namespace StabilizationMethodTestUtilities
{
void InitializeResidualBasedFluxCorrectedConstants(
    ProcessInfo& rProcessInfo);

void InitializeCrossWindStabilizationConstants(
    ProcessInfo& rProcessInfo);
} // namespace StabilizationMethodTestUtilities
} // namespace Kratos

#endif // KRATOS_STABILIZATION_METHOD_TEST_UTILITIES_H_INCLUDED