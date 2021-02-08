//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//
//

// System includes

// External includes

// Project includes
#include "includes/checks.h"
#include "testing/testing.h"

// Application includes
#include "custom_utilities/method_utilities.h"

namespace Kratos
{
namespace Testing
{
KRATOS_TEST_CASE_IN_SUITE(RaiseToPower_Int, KratosStatisticsFastSuite)
{
    KRATOS_TRY

    const int test_value = 10;
    const int return_value = MethodUtilities::RaiseToPower(test_value, 2);
    KRATOS_CHECK_EQUAL(std::pow(test_value, 2), return_value);

    KRATOS_CATCH("");
}

KRATOS_TEST_CASE_IN_SUITE(RaiseToPower_Double, KratosStatisticsFastSuite)
{
    KRATOS_TRY

    const double test_value = 10.2;
    const double return_value = MethodUtilities::RaiseToPower(test_value, 2);
    KRATOS_CHECK_NEAR(std::pow(test_value, 2), return_value, 1e-12);

    KRATOS_CATCH("");
}

KRATOS_TEST_CASE_IN_SUITE(RaiseToPower_Array, KratosStatisticsFastSuite)
{
    KRATOS_TRY

    array_1d<double, 3> test_value;
    test_value[0] = 2.0;
    test_value[1] = 3.1;
    test_value[2] = 6.4;
    const array_1d<double, 3>& return_value =
        MethodUtilities::RaiseToPower(test_value, 2);

    array_1d<double, 3> analytical_value;
    analytical_value[0] = std::pow(test_value[0], 2);
    analytical_value[1] = std::pow(test_value[1], 2);
    analytical_value[2] = std::pow(test_value[2], 2);

    KRATOS_CHECK_VECTOR_NEAR(analytical_value, return_value, 1e-12);

    KRATOS_CATCH("");
}

KRATOS_TEST_CASE_IN_SUITE(RaiseToPower_Vector, KratosStatisticsFastSuite)
{
    KRATOS_TRY

    Vector test_value(6);
    test_value[0] = 2.0;
    test_value[1] = 3.1;
    test_value[2] = 6.4;
    test_value[3] = 2.0;
    test_value[4] = 3.1;
    test_value[5] = 6.4;
    const Vector& return_value = MethodUtilities::RaiseToPower(test_value, 2);

    Vector analytical_value(6);
    analytical_value[0] = std::pow(test_value[0], 2);
    analytical_value[1] = std::pow(test_value[1], 2);
    analytical_value[2] = std::pow(test_value[2], 2);
    analytical_value[3] = std::pow(test_value[3], 2);
    analytical_value[4] = std::pow(test_value[4], 2);
    analytical_value[5] = std::pow(test_value[5], 2);

    KRATOS_CHECK_VECTOR_NEAR(analytical_value, return_value, 1e-12);

    KRATOS_CATCH("");
}

KRATOS_TEST_CASE_IN_SUITE(RaiseToPower_Matrix, KratosStatisticsFastSuite)
{
    KRATOS_TRY

    Matrix test_value{2, 2};
    test_value(0, 0) = 2.0;
    test_value(0, 1) = 3.0;
    test_value(1, 0) = 4.0;
    test_value(1, 1) = 6.0;
    const Matrix& return_value = MethodUtilities::RaiseToPower(test_value, 2);

    Matrix analytical_value{2, 2};
    analytical_value(0, 0) = std::pow(test_value(0, 0), 2);
    analytical_value(0, 1) = std::pow(test_value(0, 1), 2);
    analytical_value(1, 0) = std::pow(test_value(1, 0), 2);
    analytical_value(1, 1) = std::pow(test_value(1, 1), 2);

    KRATOS_CHECK_MATRIX_NEAR(analytical_value, return_value, 1e-12);

    KRATOS_CATCH("");
}
} // namespace Testing
} // namespace Kratos