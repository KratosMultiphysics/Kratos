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
#include "csv-parser/include/csv.hpp"

namespace Kratos {

namespace Testing {

using namespace csv;
using namespace csv::internals;

KRATOS_TEST_CASE_IN_SUITE(RecognizeIntegersProperly, KratosExternalLibrariesFastSuite )
{
    std::string a("1"), b(" 2018   "), c(" -69 ");
    long double out;

    KRATOS_CHECK_EQUAL(data_type(a, &out),  CSV_INT8);
    KRATOS_CHECK_EQUAL(out, 1);

    KRATOS_CHECK_EQUAL(data_type(b, &out), CSV_INT16);
    KRATOS_CHECK_EQUAL(out, 2018);

    KRATOS_CHECK_EQUAL(data_type(c, &out), CSV_INT8);
    KRATOS_CHECK_EQUAL(out, -69);
}

KRATOS_TEST_CASE_IN_SUITE(RecognizeNullProperly, KratosExternalLibrariesFastSuite )
{
    std::string null_str("");
    KRATOS_CHECK_EQUAL( data_type(null_str),  CSV_NULL );
}

KRATOS_TEST_CASE_IN_SUITE(RecognizeFloatsProperly, KratosExternalLibrariesFastSuite )
{
    std::string float_a("3.14"),
        float_b("       -3.14            "),
        e("2.71828");

    long double out;

    KRATOS_CHECK_EQUAL(data_type(float_a, &out), CSV_DOUBLE);
    KRATOS_CHECK(is_equal(out, 3.14L));

    KRATOS_CHECK_EQUAL(data_type(float_b, &out),  CSV_DOUBLE);
    KRATOS_CHECK(is_equal(out, -3.14L));

    KRATOS_CHECK_EQUAL(data_type(e, &out), CSV_DOUBLE);
    KRATOS_CHECK(is_equal(out, 2.71828L));
}

KRATOS_TEST_CASE_IN_SUITE(IntegerSizeRecognition, KratosExternalLibrariesFastSuite)
{
    std::string s;
    long double out;

    // Boundary Values
    {
        s = std::to_string((long long)csv::internals::CSV_INT8_MAX);
        KRATOS_CHECK_EQUAL(data_type(s, &out), CSV_INT8);
        KRATOS_CHECK_EQUAL(out, (long long)CSV_INT8_MAX);

        s = std::to_string((long long)csv::internals::CSV_INT16_MAX);
        KRATOS_CHECK_EQUAL(data_type(s, &out), CSV_INT16);
        KRATOS_CHECK_EQUAL(out, (long long)CSV_INT16_MAX);

        s = std::to_string((long long)csv::internals::CSV_INT32_MAX);
        KRATOS_CHECK_EQUAL(data_type(s, &out), CSV_INT32);
        KRATOS_CHECK_EQUAL(out, (long long)CSV_INT32_MAX);

        // Note: data_type() doesn't have enough precision for CSV_INT64
    }

    // Integer Overflow
    {
        s = std::to_string((long long)csv::internals::CSV_INT16_MAX + 1);
        KRATOS_CHECK_EQUAL(data_type(s, &out), CSV_INT32);
        KRATOS_CHECK_EQUAL(out, (long long)CSV_INT16_MAX + 1);

        s = std::to_string((long long)csv::internals::CSV_INT32_MAX + 1);
        KRATOS_CHECK_EQUAL(data_type(s, &out), CSV_INT64);
        KRATOS_CHECK_EQUAL(out, (long long)CSV_INT32_MAX + 1);

        // Case: Integer too large to fit in int64 --> store in long double
        s = std::to_string((long long)csv::internals::CSV_INT64_MAX);
        s.append("1");
        KRATOS_CHECK_EQUAL(data_type(s, &out), CSV_DOUBLE);
    }
}

KRATOS_TEST_CASE_IN_SUITE(RecognizeSubUnitDoubleValues, KratosExternalLibrariesFastSuite )
{
    std::string s("0.15");
    long double out;
    KRATOS_CHECK_EQUAL(data_type(s, &out), CSV_DOUBLE);
    KRATOS_CHECK(is_equal(out, 0.15L));
}

KRATOS_TEST_CASE_IN_SUITE(RecognizeDoubleValues, KratosExternalLibrariesFastSuite )
{
    // Test converting double values back and forth
    long double out = -1.0;
    std::string s;

    for (long double i = 0; i <= 2.0; i += 0.01) {
        s = std::to_string(i);
        KRATOS_CHECK_EQUAL(data_type(s, &out), CSV_DOUBLE);
        KRATOS_CHECK(is_equal(out, i));
    }
}

KRATOS_TEST_CASE_IN_SUITE(ParseScientificNotation, KratosExternalLibrariesFastSuite)
{
    // Test parsing e notation
    long double out;

    KRATOS_CHECK_EQUAL(data_type("1E-06", &out), CSV_DOUBLE);
    KRATOS_CHECK(is_equal(out, 0.000001L));

    KRATOS_CHECK_EQUAL(data_type("1e-06", &out), CSV_DOUBLE);
    KRATOS_CHECK(is_equal(out, 0.000001L));

    KRATOS_CHECK_EQUAL(data_type("2.17222E+02", &out), CSV_DOUBLE);
    KRATOS_CHECK(is_equal(out, 217.222L));

    KRATOS_CHECK_EQUAL(data_type("4.55E+10", &out), CSV_DOUBLE);
    KRATOS_CHECK(is_equal(out, 45500000000.0L));

    KRATOS_CHECK_EQUAL(data_type("4.55E+11", &out), CSV_DOUBLE);
    KRATOS_CHECK(is_equal(out, 455000000000.0L));

    KRATOS_CHECK_EQUAL(data_type("4.55E-1", &out), CSV_DOUBLE);
    KRATOS_CHECK(is_equal(out, 0.455L));

    KRATOS_CHECK_EQUAL(data_type("4.55E-5", &out), CSV_DOUBLE);
    KRATOS_CHECK(is_equal(out, 0.0000455L));

    KRATOS_CHECK_EQUAL(data_type("4.55E-000000000005", &out), CSV_DOUBLE);
    KRATOS_CHECK(is_equal(out, 0.0000455L));
}

} // namespace Testing.
} // namespace Kratos.
