//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vincent La (This an adaptation of https://github.com/vincentlaucsb/csv-parser)
//                   Vicente Mataix Ferrandiz
//
//

// System includes
#include <cmath>
#include <iostream>

// External includes

// Project includes
#include "testing/testing.h"
#include "csv_parser/include/csv.hpp"

namespace Kratos {

namespace Testing {

using namespace csv;

// TEMPLATE_TEST_CASE(CSVFieldget<>-StringValue, KratosExternalLibrariesFastSuite,
//     signed char, short int, int, long long int, double, long double) {
//     CSVField field("applesauce");
//     KRATOS_CHECK(field.get<>() == "applesauce");
//
//     // Assert that improper conversions attempts are thwarted
//     bool ex_caught = false;
//     try {
//         field.get<TestType>();
//     }
//     catch (std::runtime_error& err) {
//         KRATOS_CHECK(err.what() == csv::internals::ERROR_NAN);
//         ex_caught = true;
//     }
//
//     KRATOS_CHECK(ex_caught);
// }
//
// KRATOS_TEST_CASE_IN_SUITE(CSVFieldget<>-ErrorMessages, KratosExternalLibrariesFastSuite) {
//     CSVField field("applesauce");
//
//     bool ex_caught = false;
//     try {
//         field.get<double>();
//     }
//     catch (std::runtime_error& err) {
//         KRATOS_CHECK(err.what() == csv::internals::ERROR_NAN);
//         ex_caught = true;
//     }
//
//     KRATOS_CHECK(ex_caught);
// }
//
// KRATOS_TEST_CASE_IN_SUITE(CSVFieldget<>()-IntegralValue, KratosExternalLibrariesFastSuite) {
//     CSVField this_year("2019");
//     KRATOS_CHECK(this_year.get<>() == "2019");
//     KRATOS_CHECK(this_year.get<csv::string_view>() == "2019");
//     KRATOS_CHECK(this_year.get<int>() == 2019);
//     KRATOS_CHECK(this_year.get<long long int>() == 2019);
//     KRATOS_CHECK(this_year.get<float>() == 2019.0f);
//     KRATOS_CHECK(this_year.get<double>() == 2019.0);
//     KRATOS_CHECK(this_year.get<long double>() == 2019l);
//
//     bool ex_caught = false;
//     try {
//         this_year.get<signed char>();
//     }
//     catch (std::runtime_error& err) {
//         KRATOS_CHECK(err.what() == csv::internals::ERROR_OVERFLOW);
//         ex_caught = true;
//     }
//
//     KRATOS_CHECK(ex_caught);
// }
//
// KRATOS_TEST_CASE_IN_SUITE(CSVFieldget<>()-IntegerBoundaryValue, KratosExternalLibrariesFastSuite) {
//     // Note: Tests may fail if compiler defines typenames differently than
//     // Microsoft/GCC/clang
//     KRATOS_CHECK(CSVField("127").get<signed char>() == 127);
//     KRATOS_CHECK(CSVField("32767").get<short>() == 32767);
//     KRATOS_CHECK(CSVField("2147483647").get<int>() == 2147483647);
//
//     KRATOS_CHECK(CSVField("255").get<unsigned char>() == 255);
//     KRATOS_CHECK(CSVField("65535").get<unsigned short>() == 65535);
//     KRATOS_CHECK(CSVField("4294967295").get<unsigned>() == 4294967295);
// }
//
// // Test converting a small integer to unsigned and signed integer types
// TEMPLATE_TEST_CASE(CSVFieldget<>()-IntegralValuetoInt, KratosExternalLibrariesFastSuite,
//     unsigned char, unsigned short, unsigned int, unsigned long long,
//     char, short, int, long long int) {
//     CSVField savage("21");
//     KRATOS_CHECK(savage.get<TestType>() == 21);
// }
//
// KRATOS_TEST_CASE_IN_SUITE(CSVFieldget<>()-FloatingPointValue, KratosExternalLibrariesFastSuite) {
//     CSVField euler("2.718");
//     KRATOS_CHECK(euler.get<>() == "2.718");
//     KRATOS_CHECK(euler.get<csv::string_view>() == "2.718");
//     KRATOS_CHECK(euler.get<float>() == 2.718f);
//     KRATOS_CHECK(euler.get<double>() == 2.718);
//     KRATOS_CHECK(euler.get<long double>() == 2.718l);
// }
//
// TEMPLATE_TEST_CASE(CSVFieldget<>()-DisallowFloattoInt, KratosExternalLibrariesFastSuite) {
//     CSVField euler("2.718");
//     bool ex_caught = false;
//
//     try {
//         euler.get<TestType>();
//     }
//     catch (std::runtime_error& err) {
//         KRATOS_CHECK(err.what() == csv::internals::ERROR_FLOAT_TO_INT);
//         ex_caught = true;
//     }
//
//     KRATOS_CHECK(ex_caught);
// }
//
// TEMPLATE_TEST_CASE(CSVFieldget<>()-DisallowNegativetoUnsigned, KratosExternalLibrariesFastSuite) {
//     CSVField neg("-1337");
//     bool ex_caught = false;
//
//     try {
//         neg.get<TestType>();
//     }
//     catch (std::runtime_error& err) {
//         KRATOS_CHECK(err.what() == csv::internals::ERROR_NEG_TO_UNSIGNED);
//         ex_caught = true;
//     }
//
//     KRATOS_CHECK(ex_caught);
// }
//
// KRATOS_TEST_CASE_IN_SUITE("CSVFieldEqualityOperator, KratosExternalLibrariesFastSuite) {
//     CSVField field("3.14");
//     KRATOS_CHECK(field == "3.14");
//     KRATOS_CHECK(field == 3.14f);
//     KRATOS_CHECK(field == 3.14);
// }

} // namespace Testing.
} // namespace Kratos.
