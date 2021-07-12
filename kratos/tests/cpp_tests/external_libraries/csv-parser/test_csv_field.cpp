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
#include "csv-parser/csv.hpp"

namespace Kratos {

namespace Testing {

KRATOS_TEST_CASE_IN_SUITE(CSVFieldgetStringValue, KratosExternalLibrariesFastSuite)
{
    csv::CSVField field("applesauce");
    KRATOS_CHECK_EQUAL(field.get<>(), "applesauce");

    // Assert that improper conversions attempts are thwarted
    bool ex_caught = false;
    try {
        field.get<int>();
    }
    catch (std::runtime_error& err) {
        KRATOS_CHECK_EQUAL(err.what(), csv::internals::ERROR_NAN);
        ex_caught = true;
    }

    KRATOS_CHECK(ex_caught);
}

KRATOS_TEST_CASE_IN_SUITE(CSVFieldgetErrorMessages, KratosExternalLibrariesFastSuite)
{
    csv::CSVField field("applesauce");

    bool ex_caught = false;
    try {
        field.get<double>();
    }
    catch (std::runtime_error& err) {
        KRATOS_CHECK_EQUAL(err.what(), csv::internals::ERROR_NAN);
        ex_caught = true;
    }

    KRATOS_CHECK(ex_caught);
}

KRATOS_TEST_CASE_IN_SUITE(CSVFieldgetIntegralValue, KratosExternalLibrariesFastSuite)
{
    csv::CSVField this_year("2019");
    KRATOS_CHECK_EQUAL(this_year.get<>(), "2019");
    KRATOS_CHECK_EQUAL(this_year.get<csv::string_view>(), "2019");
    KRATOS_CHECK_EQUAL(this_year.get<int>(), 2019);
    KRATOS_CHECK_EQUAL(this_year.get<long long int>(), 2019);
    KRATOS_CHECK_EQUAL(this_year.get<float>(), 2019.0f);
    KRATOS_CHECK_EQUAL(this_year.get<double>(), 2019.0);
    KRATOS_CHECK_EQUAL(this_year.get<long double>(), 2019l);

    bool ex_caught = false;
    try {
        this_year.get<signed char>();
    }
    catch (std::runtime_error& err) {
        KRATOS_CHECK_EQUAL(err.what(), csv::internals::ERROR_OVERFLOW);
        ex_caught = true;
    }

    KRATOS_CHECK(ex_caught);
}

KRATOS_TEST_CASE_IN_SUITE(CSVFieldgetIntegerBoundaryValue, KratosExternalLibrariesFastSuite)
{
    // Note: Tests may fail if compiler defines typenames differently than
    // Microsoft/GCC/clang
    KRATOS_CHECK_EQUAL(csv::CSVField("127").get<signed char>(), 127);
    KRATOS_CHECK_EQUAL(csv::CSVField("32767").get<short>(), 32767);
    KRATOS_CHECK_EQUAL(csv::CSVField("2147483647").get<int>(), 2147483647);

    KRATOS_CHECK_EQUAL(csv::CSVField("255").get<unsigned char>(), 255);
    KRATOS_CHECK_EQUAL(csv::CSVField("65535").get<unsigned short>(), 65535);
    KRATOS_CHECK_EQUAL(csv::CSVField("4294967295").get<unsigned>(), 4294967295);
}

// Test converting a small integer to unsigned and signed integer types
KRATOS_TEST_CASE_IN_SUITE(CSVFieldgetIntegralValuetoInt, KratosExternalLibrariesFastSuite)
{
    csv::CSVField savage("21");
    KRATOS_CHECK_EQUAL(savage.get<int>(), 21);
}

KRATOS_TEST_CASE_IN_SUITE(CSVFieldgetFloatingPointValue, KratosExternalLibrariesFastSuite)
{
    csv::CSVField euler("2.718");
    KRATOS_CHECK_EQUAL(euler.get<>(), "2.718");
    KRATOS_CHECK_EQUAL(euler.get<csv::string_view>(), "2.718");
    KRATOS_CHECK_EQUAL(euler.get<float>(), 2.718f);
    KRATOS_CHECK_EQUAL(euler.get<double>(), 2.718);
    KRATOS_CHECK_EQUAL(euler.get<long double>(), 2.718l);
}

KRATOS_TEST_CASE_IN_SUITE(CSVFieldgetDisallowFloattoInt, KratosExternalLibrariesFastSuite)
{
    csv::CSVField euler("2.718");
    bool ex_caught = false;

    try {
        euler.get<int>();
    }
    catch (std::runtime_error& err) {
        KRATOS_CHECK_EQUAL(err.what(), csv::internals::ERROR_FLOAT_TO_INT);
        ex_caught = true;
    }

    KRATOS_CHECK(ex_caught);
}

KRATOS_TEST_CASE_IN_SUITE(CSVFieldgetDisallowNegativetoUnsigned, KratosExternalLibrariesFastSuite) {
    csv::CSVField neg("-1337");
    bool ex_caught = false;

    try {
        neg.get<std::size_t>();
    }
    catch (std::runtime_error& err) {
        KRATOS_CHECK_EQUAL(err.what(), csv::internals::ERROR_NEG_TO_UNSIGNED);
        ex_caught = true;
    }

    KRATOS_CHECK(ex_caught);
}

KRATOS_TEST_CASE_IN_SUITE(CSVFieldEqualityOperator, KratosExternalLibrariesFastSuite)
{
    csv::CSVField field("3.14");
    KRATOS_CHECK_EQUAL(field, "3.14");
    KRATOS_CHECK_EQUAL(field, 3.14f);
    KRATOS_CHECK_EQUAL(field, 3.14);
}

} // namespace Testing.
} // namespace Kratos.
