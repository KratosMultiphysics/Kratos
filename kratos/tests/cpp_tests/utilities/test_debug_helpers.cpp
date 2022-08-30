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
//

// Project includes
#include "includes/debug_helpers.h"
#include "includes/ublas_interface.h"
#include "testing/testing.h"

namespace Kratos
{
namespace Testing
{
KRATOS_TEST_CASE_IN_SUITE(DebugHelpersKRATOS_WATCH_VECTOR_WITH_PRECISION, KratosCoreFastSuite)
{
    // This test is added merely to check debug_helpers.h is compilable since this header is not supposed to be
    // included in any of the production codes.
    Vector temp(1, 10.0);
    KRATOS_WATCH_VECTOR_WITH_PRECISION(temp, 12);
}

KRATOS_TEST_CASE_IN_SUITE(DebugHelpersKRATOS_WATCH_MATRIX_WITH_PRECISION, KratosCoreFastSuite)
{
    // This test is added merely to check debug_helpers.h is compilable since this header is not supposed to be
    // included in any of the production codes.    
    Matrix temp(1, 1, 10.0);
    KRATOS_WATCH_MATRIX_WITH_PRECISION(temp, 12);
}
} // namespace Testing
} // namespace Kratos
