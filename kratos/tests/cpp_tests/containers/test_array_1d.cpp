//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
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

KRATOS_TEST_CASE_IN_SUITE(Array1DInitializationValue, KratosCoreFastSuite) {
    array_1d<double, 3> arr(2, 2.2);
    Vector ref(3);
    ref[0] = 2.2;
    ref[1] = 2.2;
    KRATOS_CHECK_DOUBLE_EQUAL(arr[0], ref[0]);
    KRATOS_CHECK_DOUBLE_EQUAL(arr[1], ref[1]);
}

KRATOS_TEST_CASE_IN_SUITE(Array1DInitializerList, KratosCoreFastSuite) {
    array_1d<double, 3> arr {1.1, -2.3, 3.4};
    Vector ref(3);
    ref[0] = 1.1;
    ref[1] = -2.3;
    ref[2] = 3.4;
    KRATOS_CHECK_VECTOR_EQUAL(arr, ref);
}

} // namespace Testing.
} // namespace Kratos.
