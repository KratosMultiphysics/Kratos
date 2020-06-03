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
#include "csv.hpp"

namespace Kratos {

namespace Testing {

using namespace csv;

// TEMPLATE_TEST_CASE("CSVField get<> - String Value", "[test_csv_field_get_string]",
//     signed char, short int, int, long long int, double, long double) {
//     CSVField field("applesauce");
//     REQUIRE(field.get<>() == "applesauce");
//
//     // Assert that improper conversions attempts are thwarted
//     bool ex_caught = false;
//     try {
//         field.get<TestType>();
//     }
//     catch (std::runtime_error& err) {
//         REQUIRE(err.what() == csv::internals::ERROR_NAN);
//         ex_caught = true;
//     }
//
//     REQUIRE(ex_caught);
// }
//
// TEST_CASE("CSVField get<> - Error Messages", "[test_csv_field_get_error]") {
//     CSVField field("applesauce");
//
//     bool ex_caught = false;
//     try {
//         field.get<double>();
//     }
//     catch (std::runtime_error& err) {
//         REQUIRE(err.what() == csv::internals::ERROR_NAN);
//         ex_caught = true;
//     }
//
//     REQUIRE(ex_caught);
// }
//
// TEST_CASE("CSVField get<>() - Integral Value", "[test_csv_field_get_int]") {
//     CSVField this_year("2019");
//     REQUIRE(this_year.get<>() == "2019");
//     REQUIRE(this_year.get<csv::string_view>() == "2019");
//     REQUIRE(this_year.get<int>() == 2019);
//     REQUIRE(this_year.get<long long int>() == 2019);
//     REQUIRE(this_year.get<float>() == 2019.0f);
//     REQUIRE(this_year.get<double>() == 2019.0);
//     REQUIRE(this_year.get<long double>() == 2019l);
//
//     bool ex_caught = false;
//     try {
//         this_year.get<signed char>();
//     }
//     catch (std::runtime_error& err) {
//         REQUIRE(err.what() == csv::internals::ERROR_OVERFLOW);
//         ex_caught = true;
//     }
//
//     REQUIRE(ex_caught);
// }
//
// TEST_CASE("CSVField get<>() - Integer Boundary Value", "[test_csv_field_get_boundary]") {
//     // Note: Tests may fail if compiler defines typenames differently than
//     // Microsoft/GCC/clang
//     REQUIRE(CSVField("127").get<signed char>() == 127);
//     REQUIRE(CSVField("32767").get<short>() == 32767);
//     REQUIRE(CSVField("2147483647").get<int>() == 2147483647);
//
//     REQUIRE(CSVField("255").get<unsigned char>() == 255);
//     REQUIRE(CSVField("65535").get<unsigned short>() == 65535);
//     REQUIRE(CSVField("4294967295").get<unsigned>() == 4294967295);
// }
//
// // Test converting a small integer to unsigned and signed integer types
// TEMPLATE_TEST_CASE("CSVField get<>() - Integral Value to Int", "[test_csv_field_convert_int]",
//     unsigned char, unsigned short, unsigned int, unsigned long long,
//     char, short, int, long long int) {
//     CSVField savage("21");
//     REQUIRE(savage.get<TestType>() == 21);
// }
//
// TEST_CASE("CSVField get<>() - Floating Point Value", "[test_csv_field_get_float]") {
//     CSVField euler("2.718");
//     REQUIRE(euler.get<>() == "2.718");
//     REQUIRE(euler.get<csv::string_view>() == "2.718");
//     REQUIRE(euler.get<float>() == 2.718f);
//     REQUIRE(euler.get<double>() == 2.718);
//     REQUIRE(euler.get<long double>() == 2.718l);
// }
//
// TEMPLATE_TEST_CASE("CSVField get<>() - Disallow Float to Int", "[test_csv_field_get_float_as_int]",
//     unsigned char, unsigned short, unsigned int, unsigned long long int,
//     signed char, short, int, long long int) {
//     CSVField euler("2.718");
//     bool ex_caught = false;
//
//     try {
//         euler.get<TestType>();
//     }
//     catch (std::runtime_error& err) {
//         REQUIRE(err.what() == csv::internals::ERROR_FLOAT_TO_INT);
//         ex_caught = true;
//     }
//
//     REQUIRE(ex_caught);
// }
//
// TEMPLATE_TEST_CASE("CSVField get<>() - Disallow Negative to Unsigned", "[test_csv_field_no_unsigned_neg]",
//     unsigned char, unsigned short, unsigned int, unsigned long long int) {
//     CSVField neg("-1337");
//     bool ex_caught = false;
//
//     try {
//         neg.get<TestType>();
//     }
//     catch (std::runtime_error& err) {
//         REQUIRE(err.what() == csv::internals::ERROR_NEG_TO_UNSIGNED);
//         ex_caught = true;
//     }
//
//     REQUIRE(ex_caught);
// }
//
// TEST_CASE("CSVField Equality Operator", "[test_csv_field_operator==]") {
//     CSVField field("3.14");
//     REQUIRE(field == "3.14");
//     REQUIRE(field == 3.14f);
//     REQUIRE(field == 3.14);
// }

} // namespace Testing.
} // namespace Kratos.
