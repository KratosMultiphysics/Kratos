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
#include <string>

// External includes

// Project includes
#include "testing/testing.h"
#include "csv_parser/include/csv.hpp"

namespace Kratos {

namespace Testing {

using namespace csv;
using namespace csv::internals;

// TEST_CASE( "Recognize Integers Properly", "[dtype_int]" ) {
//     std::string a("1"), b(" 2018   "), c(" -69 ");
//     long double out;
//
//     REQUIRE(data_type(a, &out) ==  CSV_INT8);
//     REQUIRE(out == 1);
//
//     REQUIRE(data_type(b, &out) == CSV_INT16);
//     REQUIRE(out == 2018);
//
//     REQUIRE(data_type(c, &out) == CSV_INT8);
//     REQUIRE(out == -69);
// }
//
// TEST_CASE( "Recognize Strings Properly", "[dtype_str]" ) {
//     auto str = GENERATE(as<std::string> {}, "test", "999.999.9999", "510-123-4567", "510 123", "510 123 4567");
//
//     SECTION("String Recognition") {
//         REQUIRE(data_type(str) == CSV_STRING);
//     }
// }
//
// TEST_CASE( "Recognize Null Properly", "[dtype_null]" ) {
//     std::string null_str("");
//     REQUIRE( data_type(null_str) ==  CSV_NULL );
// }
//
// TEST_CASE( "Recognize Floats Properly", "[dtype_float]" ) {
//     std::string float_a("3.14"),
//         float_b("       -3.14            "),
//         e("2.71828");
//
//     long double out;
//
//     REQUIRE(data_type(float_a, &out) == CSV_DOUBLE);
//     REQUIRE(is_equal(out, 3.14L));
//
//     REQUIRE(data_type(float_b, &out) ==  CSV_DOUBLE);
//     REQUIRE(is_equal(out, -3.14L));
//
//     REQUIRE(data_type(e, &out) == CSV_DOUBLE);
//     REQUIRE(is_equal(out, 2.71828L));
// }
//
// TEST_CASE("Integer Size Recognition", "[int_sizes]") {
//     std::string s;
//     long double out;
//
//     SECTION("Boundary Values") {
//         s = std::to_string((long long)csv::internals::CSV_INT8_MAX);
//         REQUIRE(data_type(s, &out) == CSV_INT8);
//         REQUIRE(out == (long long)CSV_INT8_MAX);
//
//         s = std::to_string((long long)csv::internals::CSV_INT16_MAX);
//         REQUIRE(data_type(s, &out) == CSV_INT16);
//         REQUIRE(out == (long long)CSV_INT16_MAX);
//
//         s = std::to_string((long long)csv::internals::CSV_INT32_MAX);
//         REQUIRE(data_type(s, &out) == CSV_INT32);
//         REQUIRE(out == (long long)CSV_INT32_MAX);
//
//         // Note: data_type() doesn't have enough precision for CSV_INT64
//     }
//
//     SECTION("Integer Overflow") {
//         s = std::to_string((long long)csv::internals::CSV_INT16_MAX + 1);
//         REQUIRE(data_type(s, &out) == CSV_INT32);
//         REQUIRE(out == (long long)CSV_INT16_MAX + 1);
//
//         s = std::to_string((long long)csv::internals::CSV_INT32_MAX + 1);
//         REQUIRE(data_type(s, &out) == CSV_INT64);
//         REQUIRE(out == (long long)CSV_INT32_MAX + 1);
//
//         // Case: Integer too large to fit in int64 --> store in long double
//         s = std::to_string((long long)csv::internals::CSV_INT64_MAX);
//         s.append("1");
//         REQUIRE(data_type(s, &out) == CSV_DOUBLE);
//     }
// }
//
// TEST_CASE( "Recognize Sub-Unit Double Values", "[regression_double]" ) {
//     std::string s("0.15");
//     long double out;
//     REQUIRE(data_type(s, &out) == CSV_DOUBLE);
//     REQUIRE(is_equal(out, 0.15L));
// }
//
// TEST_CASE( "Recognize Double Values", "[regression_double2]" ) {
//     // Test converting double values back and forth
//     long double out = -1.0;
//     std::string s;
//
//     for (long double i = 0; i <= 2.0; i += 0.01) {
//         s = std::to_string(i);
//         REQUIRE(data_type(s, &out) == CSV_DOUBLE);
//         REQUIRE(is_equal(out, i));
//     }
// }
//
// //! [Parse Scientific Notation]
// TEST_CASE("Parse Scientific Notation", "[e_notation]") {
//     // Test parsing e notation
//     long double out;
//
//     REQUIRE(data_type("1E-06", &out) == CSV_DOUBLE);
//     REQUIRE(is_equal(out, 0.000001L));
//
//     REQUIRE(data_type("1e-06", &out) == CSV_DOUBLE);
//     REQUIRE(is_equal(out, 0.000001L));
//
//     REQUIRE(data_type("2.17222E+02", &out) == CSV_DOUBLE);
//     REQUIRE(is_equal(out, 217.222L));
//
//     REQUIRE(data_type("4.55E+10", &out) == CSV_DOUBLE);
//     REQUIRE(is_equal(out, 45500000000.0L));
//
//     REQUIRE(data_type("4.55E+11", &out) == CSV_DOUBLE);
//     REQUIRE(is_equal(out, 455000000000.0L));
//
//     REQUIRE(data_type("4.55E-1", &out) == CSV_DOUBLE);
//     REQUIRE(is_equal(out, 0.455L));
//
//     REQUIRE(data_type("4.55E-5", &out) == CSV_DOUBLE);
//     REQUIRE(is_equal(out, 0.0000455L));
//
//     REQUIRE(data_type("4.55E-000000000005", &out) == CSV_DOUBLE);
//     REQUIRE(is_equal(out, 0.0000455L));
// }
// //! [Parse Scientific Notation]
//
// //! [Scientific Notation Flavors]
// TEST_CASE("Parse Different Flavors of Scientific Notation", "[sci_notation_diversity]") {
//     auto number = GENERATE(as<std::string> {},
//         "4.55e5", "4.55E5",
//         "4.55E+5", "4.55e+5",
//         "4.55E+05",
//         "4.55e0000005", "4.55E0000005",
//         "4.55e+0000005", "4.55E+0000005");
//
//     SECTION("Recognize 455 thousand") {
//         long double out;
//         REQUIRE(data_type(number, &out) == CSV_DOUBLE);
//         REQUIRE(is_equal(out, 455000.0L));
//     }
// }
// //! [Scientific Notation Flavors]
//
// TEST_CASE("Parse Scientific Notation Malformed", "[sci_notation]") {
//     // Assert parsing butchered scientific notation won't cause a
//     // crash or any other weird side effects
//     auto butchered = GENERATE(as<std::string>{},
//         "4.55E000a",
//         "4.55000x40",
//         "4.55000E40E40");
//
//     SECTION("Butchered Parsing Attempt") {
//         REQUIRE(data_type(butchered) == CSV_STRING);
//     }
// }

} // namespace Testing.
} // namespace Kratos.
