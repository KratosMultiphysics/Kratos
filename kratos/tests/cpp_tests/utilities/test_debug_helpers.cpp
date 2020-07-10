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
    Vector temp(4, 10.0);
    KRATOS_WATCH_VECTOR_WITH_PRECISION(temp, 12);
}

KRATOS_TEST_CASE_IN_SUITE(DebugHelpersKRATOS_WATCH_MATRIX_WITH_PRECISION, KratosCoreFastSuite)
{
    Matrix temp(4, 4, 10.0);
    KRATOS_WATCH_MATRIX_WITH_PRECISION(temp, 12);
}
} // namespace Testing
} // namespace Kratos
