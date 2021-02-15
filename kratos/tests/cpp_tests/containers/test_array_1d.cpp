//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license:
//                   kratos/license.txt
//
//  Main authors:    Philipp Bucher (https://github.com/philbucher)
//
//

// System includes

// External includes

// Project includes
#include "containers/array_1d.h"
#include "testing/testing.h"

namespace Kratos {
namespace Testing {

KRATOS_TEST_CASE_IN_SUITE(Array1DInitialization, KratosCoreFastSuite) {
    array_1d<double, 3> arr;
    KRATOS_CHECK_VECTOR_EQUAL(arr, ZeroVector(3));
}

KRATOS_TEST_CASE_IN_SUITE(Array1DInitializationSize, KratosCoreFastSuite) {
    array_1d<double, 3> arr(2);
    KRATOS_CHECK_VECTOR_EQUAL(arr, ZeroVector(3)); // still initializes all entries
}

KRATOS_TEST_CASE_IN_SUITE(Array1DInitializationValue, KratosCoreFastSuite) {
    array_1d<double, 3> arr(2, 2.2);
    Vector ref(3);
    ref[0] = 2.2;
    ref[1] = 2.2;
    ref[2] = 0.0;
    KRATOS_CHECK_VECTOR_EQUAL(arr, ref);
}

KRATOS_TEST_CASE_IN_SUITE(Array1DInitializerList, KratosCoreFastSuite) {
    array_1d<double, 3> arr {1.1, -2.3, 3.4};
    Vector ref(3);
    ref[0] = 1.1;
    ref[1] = -2.3;
    ref[2] = 3.4;
    KRATOS_CHECK_VECTOR_EQUAL(arr, ref);
}

KRATOS_TEST_CASE_IN_SUITE(Array1DInitializerListShort, KratosCoreFastSuite) {
    array_1d<double, 5> arr {1.1, -2.3, 3.4}; // shorter list, other values should be default initialized
    Vector ref(5);
    ref[0] = 1.1;
    ref[1] = -2.3;
    ref[2] = 3.4;
    ref[3] = 0.0;
    ref[4] = 0.0;
    KRATOS_CHECK_VECTOR_EQUAL(arr, ref);
}

} // namespace Testing.
} // namespace Kratos.
