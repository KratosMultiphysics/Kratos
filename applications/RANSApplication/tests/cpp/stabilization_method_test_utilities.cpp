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

// External includes

// Project includes
#include "includes/cfd_variables.h"
#include "includes/process_info.h"
#include "includes/variables.h"

// Application includes
#include "rans_application_variables.h"

// Include base h
#include "stabilization_method_test_utilities.h"

namespace Kratos
{
namespace StabilizationMethodTestUtilities
{
void InitializeResidualBasedFluxCorrectedConstants(
    ProcessInfo& rProcessInfo)
{
    rProcessInfo.SetValue(DELTA_TIME, 2.6);
    rProcessInfo.SetValue(BOSSAK_ALPHA, -0.3);
    rProcessInfo.SetValue(DYNAMIC_TAU, 0.8);
    rProcessInfo.SetValue(RANS_STABILIZATION_DISCRETE_UPWIND_OPERATOR_COEFFICIENT, 1.8);
    rProcessInfo.SetValue(RANS_STABILIZATION_DIAGONAL_POSITIVITY_PRESERVING_COEFFICIENT, 2.8);
}
void InitializeCrossWindStabilizationConstants(
    ProcessInfo& rProcessInfo)
{
    rProcessInfo.SetValue(DELTA_TIME, 2.6);
    rProcessInfo.SetValue(BOSSAK_ALPHA, -0.3);
    rProcessInfo.SetValue(DYNAMIC_TAU, 0.8);
}
} // namespace StabilizationMethodTestUtilities
} // namespace Kratos
