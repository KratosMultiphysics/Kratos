//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//                   Philipp Bucher (https://github.com/philbucher)
//

// System includes
#include <tuple>

// External includes

// Project includes
#include "testing/testing.h"
#include "utilities/string_utilities.h"

namespace Kratos::Testing {

KRATOS_TEST_CASE_IN_SUITE(ConvertCamelCaseToSnakeCase, KratosCoreFastSuite)
{
    const std::string CamelCase = "TestInCamelCase";
    const std::string snake_case = StringUtilities::ConvertCamelCaseToSnakeCase(CamelCase);
    KRATOS_EXPECT_EQ(snake_case, "test_in_camel_case");
}

KRATOS_TEST_CASE_IN_SUITE(ConvertSnakeCaseToCamelCase, KratosCoreFastSuite)
{
    KRATOS_EXPECT_EQ(StringUtilities::ConvertSnakeCaseToCamelCase("test_snake_case"),
                              "TestSnakeCase");
    KRATOS_EXPECT_EQ(StringUtilities::ConvertSnakeCaseToCamelCase("test"),
                              "Test");
    KRATOS_EXPECT_EQ(StringUtilities::ConvertSnakeCaseToCamelCase("t_e_s_t"),
                              "TEST");
    KRATOS_EXPECT_EQ(StringUtilities::ConvertSnakeCaseToCamelCase("_test_"),
                              "Test");
    KRATOS_EXPECT_EQ(StringUtilities::ConvertSnakeCaseToCamelCase("_te_st_"),
                              "TeSt");
    KRATOS_EXPECT_EQ(StringUtilities::ConvertSnakeCaseToCamelCase("_"),
                              "");
    KRATOS_EXPECT_EQ(StringUtilities::ConvertSnakeCaseToCamelCase(""),
                              "");


    KRATOS_EXPECT_EQ(StringUtilities::ConvertSnakeCaseToCamelCase("num3r1c41"),
                              "Num3r1c41");
    KRATOS_EXPECT_EQ(StringUtilities::ConvertSnakeCaseToCamelCase("num3r_1c41"),
                              "Num3r1c41");
    KRATOS_EXPECT_EQ(StringUtilities::ConvertSnakeCaseToCamelCase("3141"),
                              "3141");

    #define KRATOS_EXPECT_THROWS(...)                            \
        try {                                                   \
            __VA_ARGS__;                                        \
            KRATOS_ERROR  << "The expression " << #__VA_ARGS__  \
            << " was supposed to throw, but it didn't";         \
        }                                                       \
        catch (...) {}

    // deliberately ignoring [[nodiscard]] as it is not relevant for this test
    KRATOS_EXPECT_THROWS(std::ignore = StringUtilities::ConvertSnakeCaseToCamelCase("Test"))
    KRATOS_EXPECT_THROWS(std::ignore = StringUtilities::ConvertSnakeCaseToCamelCase("tesT"))
    KRATOS_EXPECT_THROWS(std::ignore = StringUtilities::ConvertSnakeCaseToCamelCase("te__st"))
    KRATOS_EXPECT_THROWS(std::ignore = StringUtilities::ConvertSnakeCaseToCamelCase("__te_st"))
    KRATOS_EXPECT_THROWS(std::ignore = StringUtilities::ConvertSnakeCaseToCamelCase("test__"))
    KRATOS_EXPECT_THROWS(std::ignore = StringUtilities::ConvertSnakeCaseToCamelCase("te st"))
    KRATOS_EXPECT_THROWS(std::ignore = StringUtilities::ConvertSnakeCaseToCamelCase(" test"))
    KRATOS_EXPECT_THROWS(std::ignore = StringUtilities::ConvertSnakeCaseToCamelCase("test "))
    KRATOS_EXPECT_THROWS(std::ignore = StringUtilities::ConvertSnakeCaseToCamelCase("*nullptr"))
    KRATOS_EXPECT_THROWS(std::ignore = StringUtilities::ConvertSnakeCaseToCamelCase("core/stringutils"))
    KRATOS_EXPECT_THROWS(std::ignore = StringUtilities::ConvertSnakeCaseToCamelCase("c-s"))
    KRATOS_EXPECT_THROWS(std::ignore = StringUtilities::ConvertSnakeCaseToCamelCase("-cs"))
    KRATOS_EXPECT_THROWS(std::ignore = StringUtilities::ConvertSnakeCaseToCamelCase("cs-"))
    KRATOS_EXPECT_THROWS(std::ignore = StringUtilities::ConvertSnakeCaseToCamelCase("ph@"))
    KRATOS_EXPECT_THROWS(std::ignore = StringUtilities::ConvertSnakeCaseToCamelCase("#include"))

    #undef KRATOS_EXPECT_THROWS
}

KRATOS_TEST_CASE_IN_SUITE(ErasePartialString, KratosCoreFastSuite)
{
    const std::string text_with_stuff = "textwithstuff";
    std::string text = StringUtilities::ErasePartialString(text_with_stuff, "withstuff");
    KRATOS_EXPECT_EQ(text, "text");
    text = StringUtilities::ErasePartialString(text_with_stuff, "with");
    KRATOS_EXPECT_EQ(text, "textstuff");
}

KRATOS_TEST_CASE_IN_SUITE(ContainsPartialString, KratosCoreFastSuite)
{
    const std::string text_with_stuff = "textwithstuff";
    bool found = StringUtilities::ContainsPartialString(text_with_stuff, "withstuff");
    KRATOS_EXPECT_EQ(found, true);
    found = StringUtilities::ContainsPartialString(text_with_stuff, "daleatucuerpoalegriamacarena");
    KRATOS_EXPECT_EQ(found, false);
}

KRATOS_TEST_CASE_IN_SUITE(RemoveWhiteSpaces, KratosCoreFastSuite)
{
    std::string text_with_spaces = "textwithspaces                ";
    std::string text = StringUtilities::RemoveWhiteSpaces(text_with_spaces);
    KRATOS_EXPECT_EQ(text, "textwithspaces");
    text_with_spaces = "textwithspaces                jojobizarreadventure";
    text = StringUtilities::RemoveWhiteSpaces(text_with_spaces);
    KRATOS_EXPECT_EQ(text, "textwithspacesjojobizarreadventure");
}

KRATOS_TEST_CASE_IN_SUITE(SplitStringByDelimiter, KratosCoreFastSuite)
{
    const std::string string_to_split = "abc, eeee,;esdsdsd,   ";
    const std::vector<std::string> splitted_string = StringUtilities::SplitStringByDelimiter(string_to_split, ',');
    KRATOS_EXPECT_EQ(splitted_string.size(), 4);

    KRATOS_EXPECT_EQ(splitted_string[0], "abc");
    KRATOS_EXPECT_EQ(splitted_string[1], " eeee");
    KRATOS_EXPECT_EQ(splitted_string[2], ";esdsdsd");
    KRATOS_EXPECT_EQ(splitted_string[3], "   ");
}

KRATOS_TEST_CASE_IN_SUITE(ReplaceAllSubstrings, KratosCoreFastSuite)
{
    const std::string string_to_replace = "Pokemon is an awesome show. Pokemon is the best!";
    const std::string correct_string = StringUtilities::ReplaceAllSubstrings(string_to_replace, "Pokemon", "JoJo Bizarre Adventure");

    KRATOS_EXPECT_EQ(correct_string, "JoJo Bizarre Adventure is an awesome show. JoJo Bizarre Adventure is the best!");
}

KRATOS_TEST_CASE_IN_SUITE(Trim, KratosCoreFastSuite)
{
    using namespace std::string_literals; // required for strings with null-chars => "asdf\0"s

    KRATOS_EXPECT_EQ(StringUtilities::Trim(""), "");
    KRATOS_EXPECT_EQ(StringUtilities::Trim(" "), "");
    KRATOS_EXPECT_EQ(StringUtilities::Trim("\n"), "");
    KRATOS_EXPECT_EQ(StringUtilities::Trim("\t"), "");
    KRATOS_EXPECT_EQ(StringUtilities::Trim("\0"), "\0");
    KRATOS_EXPECT_EQ(StringUtilities::Trim("\0", true), "");

    KRATOS_EXPECT_EQ(StringUtilities::Trim("Kratos"), "Kratos");
    KRATOS_EXPECT_EQ(StringUtilities::Trim(" Kratos "), "Kratos");
    KRATOS_EXPECT_EQ(StringUtilities::Trim(" Kra tos "), "Kra tos");
    KRATOS_EXPECT_EQ(StringUtilities::Trim(" Kra\0tos "), "Kra\0tos");
    KRATOS_EXPECT_EQ(StringUtilities::Trim(" Kra\ntos "), "Kra\ntos");
    KRATOS_EXPECT_EQ(StringUtilities::Trim(" \0 Kra\ntos \0 "s), "\0 Kra\ntos \0"s);
    KRATOS_EXPECT_EQ(StringUtilities::Trim(" \0 Kra\ntos \0 "s, true), "Kra\ntos");

    KRATOS_EXPECT_EQ(StringUtilities::Trim("Kratos MP"), "Kratos MP");
    KRATOS_EXPECT_EQ(StringUtilities::Trim(" Kratos MP "), "Kratos MP");
    KRATOS_EXPECT_EQ(StringUtilities::Trim("\tKratos MP\n"), "Kratos MP");
    KRATOS_EXPECT_EQ(StringUtilities::Trim(" Kra tos MP "), "Kra tos MP");
    KRATOS_EXPECT_EQ(StringUtilities::Trim(" Kra\0tos MP "), "Kra\0tos MP");
    KRATOS_EXPECT_EQ(StringUtilities::Trim(" Kra\ntos MP "), "Kra\ntos MP");
}

KRATOS_TEST_CASE_IN_SUITE(TrimLeft, KratosCoreFastSuite)
{
    using namespace std::string_literals; // required for strings with null-chars => "asdf\0"s

    KRATOS_EXPECT_EQ(StringUtilities::TrimLeft(""), "");
    KRATOS_EXPECT_EQ(StringUtilities::TrimLeft(" "), "");
    KRATOS_EXPECT_EQ(StringUtilities::TrimLeft("\n"), "");
    KRATOS_EXPECT_EQ(StringUtilities::TrimLeft("\t"), "");
    KRATOS_EXPECT_EQ(StringUtilities::TrimLeft("\0"), "\0");
    KRATOS_EXPECT_EQ(StringUtilities::TrimLeft("\0", true), "");

    KRATOS_EXPECT_EQ(StringUtilities::TrimLeft("Kratos"), "Kratos");
    KRATOS_EXPECT_EQ(StringUtilities::TrimLeft(" Kratos "), "Kratos ");
    KRATOS_EXPECT_EQ(StringUtilities::TrimLeft(" Kra tos "), "Kra tos ");
    KRATOS_EXPECT_EQ(StringUtilities::TrimLeft(" Kra\0tos "), "Kra\0tos ");
    KRATOS_EXPECT_EQ(StringUtilities::TrimLeft(" Kra\ntos "), "Kra\ntos ");
    KRATOS_EXPECT_EQ(StringUtilities::TrimLeft(" Kra\ntos\0"s), "Kra\ntos\0"s);
    KRATOS_EXPECT_EQ(StringUtilities::TrimLeft(" \0 Kra\ntos \0 "s), "\0 Kra\ntos \0 "s);
    KRATOS_EXPECT_EQ(StringUtilities::TrimLeft(" \0 Kra\ntos \0 "s, true), "Kra\ntos \0 "s);

    KRATOS_EXPECT_EQ(StringUtilities::TrimLeft("Kratos MP"), "Kratos MP");
    KRATOS_EXPECT_EQ(StringUtilities::TrimLeft(" Kratos MP "), "Kratos MP ");
    KRATOS_EXPECT_EQ(StringUtilities::TrimLeft("\tKratos MP\n"), "Kratos MP\n");
    KRATOS_EXPECT_EQ(StringUtilities::TrimLeft(" Kra tos MP "), "Kra tos MP ");
    KRATOS_EXPECT_EQ(StringUtilities::TrimLeft(" Kra\0tos MP "), "Kra\0tos MP ");
    KRATOS_EXPECT_EQ(StringUtilities::TrimLeft(" Kra\ntos MP "), "Kra\ntos MP ");
}

KRATOS_TEST_CASE_IN_SUITE(TrimRight, KratosCoreFastSuite)
{
    using namespace std::string_literals; // required for strings with null-chars => "asdf\0"s

    KRATOS_EXPECT_EQ(StringUtilities::TrimRight(""), "");
    KRATOS_EXPECT_EQ(StringUtilities::TrimRight(" "), "");
    KRATOS_EXPECT_EQ(StringUtilities::TrimRight("\n"), "");
    KRATOS_EXPECT_EQ(StringUtilities::TrimRight("\t"), "");
    KRATOS_EXPECT_EQ(StringUtilities::TrimRight("\0"), "\0");
    KRATOS_EXPECT_EQ(StringUtilities::TrimRight("\0", true), "");

    KRATOS_EXPECT_EQ(StringUtilities::TrimRight("Kratos"), "Kratos");
    KRATOS_EXPECT_EQ(StringUtilities::TrimRight(" Kratos "), " Kratos");
    KRATOS_EXPECT_EQ(StringUtilities::TrimRight(" Kra tos "), " Kra tos");
    KRATOS_EXPECT_EQ(StringUtilities::TrimRight(" Kra\0tos "), " Kra\0tos");
    KRATOS_EXPECT_EQ(StringUtilities::TrimRight(" Kra\ntos "), " Kra\ntos");
    KRATOS_EXPECT_EQ(StringUtilities::TrimRight(" Kra\ntos\0"s), " Kra\ntos\0"s);
    KRATOS_EXPECT_EQ(StringUtilities::TrimRight(" Kra\ntos\0"s, true), " Kra\ntos");
    KRATOS_EXPECT_EQ(StringUtilities::TrimRight(" \0 Kra\ntos \0 "s), " \0 Kra\ntos \0"s);
    KRATOS_EXPECT_EQ(StringUtilities::TrimRight(" \0 Kra\ntos \0 "s, true), " \0 Kra\ntos"s);

    KRATOS_EXPECT_EQ(StringUtilities::TrimRight("Kratos MP"), "Kratos MP");
    KRATOS_EXPECT_EQ(StringUtilities::TrimRight(" Kratos MP "), " Kratos MP");
    KRATOS_EXPECT_EQ(StringUtilities::TrimRight("\tKratos MP\n"), "\tKratos MP");
    KRATOS_EXPECT_EQ(StringUtilities::TrimRight(" Kra tos MP "), " Kra tos MP");
    KRATOS_EXPECT_EQ(StringUtilities::TrimRight(" Kra\0tos MP "), " Kra\0tos MP");
    KRATOS_EXPECT_EQ(StringUtilities::TrimRight(" Kra\ntos MP "), " Kra\ntos MP");
}

KRATOS_TEST_CASE_IN_SUITE(StringToVector, KratosCoreFastSuite)
{
    // BasicIntegerParsing
    {
        std::string input = "[1.0, 2.5, 3.14]";
        std::vector<double> expected = {1.0, 2.5, 3.14};
        std::vector<double> actual = StringUtilities::StringToVector<double>(input);
        for (size_t i = 0; i < expected.size(); ++i) {
            KRATOS_EXPECT_EQ(actual[i], expected[i]);
        }
    }

    // ParsingWithVariousSpaces
    {
        std::string input = "[10, 20, -5]";
        std::vector<int> expected = {10, 20, -5};
        std::vector<int> actual = StringUtilities::StringToVector<int>(input);
        for (size_t i = 0; i < expected.size(); ++i) {
            KRATOS_EXPECT_EQ(actual[i], expected[i]);
        }
    }

    // SingleElementVector
    {
        std::string input = "  [ 1.1 ,  2.2,3.3  , 4.4 ]  ";
        std::vector<double> expected = {1.1, 2.2, 3.3, 4.4};
        std::vector<double> actual = StringUtilities::StringToVector<double>(input);
        for (size_t i = 0; i < expected.size(); ++i) {
            KRATOS_EXPECT_EQ(actual[i], expected[i]);
        }

        std::string input_int = " [ 5,  15,25 ] ";
        std::vector<int> expected_int = {5, 15, 25};
        std::vector<int> actual_int = StringUtilities::StringToVector<int>(input_int);
        for (size_t i = 0; i < expected_int.size(); ++i) {
            KRATOS_EXPECT_EQ(actual_int[i], expected_int[i]);
        }
    }

    // EmptyVectorInput
    {
        std::string input = "[42.0]";
        std::vector<double> expected = {42.0};
        std::vector<double> actual = StringUtilities::StringToVector<double>(input);
        KRATOS_EXPECT_EQ(actual.size(), expected.size());
        KRATOS_EXPECT_EQ(actual.size(), 1);
        KRATOS_EXPECT_EQ(actual[0], expected[0]);

        std::string input_int = "[7]";
        std::vector<int> expected_int = {7};
        std::vector<int> actual_int = StringUtilities::StringToVector<int>(input_int);
        KRATOS_EXPECT_EQ(actual_int.size(), expected_int.size());
        KRATOS_EXPECT_EQ(actual_int.size(), 1);
        KRATOS_EXPECT_EQ(actual_int[0], expected_int[0]);
    }

    // IntegerTruncationFromDoubleString
    {
        std::string input1 = "[]";
        std::vector<double> expected_empty = {};
        std::vector<double> actual1 = StringUtilities::StringToVector<double>(input1);
        KRATOS_EXPECT_EQ(actual1.size(), 0);
        KRATOS_EXPECT_EQ(actual1.size(), expected_empty.size());

        std::string input2 = "[ ]"; // Empty with space
        std::vector<int> actual2 = StringUtilities::StringToVector<int>(input2);
        KRATOS_EXPECT_EQ(actual2.size(), 0);
        KRATOS_EXPECT_EQ(actual2.size(), expected_empty.size());

        std::string input3 = ""; // Completely empty string
        std::vector<double> actual3 = StringUtilities::StringToVector<double>(input3);
        KRATOS_EXPECT_EQ(actual3.size(), 0);
        KRATOS_EXPECT_EQ(actual3.size(), expected_empty.size());
    }

    // InputWithoutBrackets
    {
        // The function is designed to remove brackets if present.
        // If not present, it should still parse the comma-separated values.
        std::string input = "1.5, 2.5, 3.5";
        std::vector<double> expected = {1.5, 2.5, 3.5};
        std::vector<double> actual = StringUtilities::StringToVector<double>(input);
        for (size_t i = 0; i < expected.size(); ++i) {
            KRATOS_EXPECT_EQ(actual[i], expected[i]);
        }
    }

    // InvalidElementMixedWithValid
    {
        // "two" cannot be parsed as double. The function currently logs a warning
        // and continues, so "two" should be skipped.
        std::string input = "[1.0, two, 3.0, four, 5.5]";
        std::vector<double> expected = {1.0, 3.0, 5.5}; // "two" and "four" should be skipped
        std::vector<double> actual = StringUtilities::StringToVector<double>(input);
        for (size_t i = 0; i < expected.size(); ++i) {
            KRATOS_EXPECT_EQ(actual[i], expected[i]);
        }
    }

    // InvalidNumberFormat
    {
        std::string input_double = "[1.2.3, 4.5]"; // "1.2.3" is invalid for double
        std::vector<double> expected_double = {4.5};
        std::vector<double> actual_double = StringUtilities::StringToVector<double>(input_double);
        for (size_t i = 0; i < expected_double.size(); ++i) {
            KRATOS_EXPECT_EQ(actual_double[i], expected_double[i]);
        }

        std::string input_int = "[10, 20x, 30]"; // "20x" is invalid for int
        std::vector<int> expected_int = {10, 30};
        std::vector<int> actual_int = StringUtilities::StringToVector<int>(input_int);
        for (size_t i = 0; i < expected_int.size(); ++i) {
            KRATOS_EXPECT_EQ(actual_int[i], expected_int[i]);
        }
    }
}

KRATOS_TEST_CASE_IN_SUITE(CountValuesUntilPrefix, KratosCoreFastSuite)
{
    // Test case 1: Basic functionality
    {
        std::string input = "value1 value2 value3";
        std::string prefix = "value2";
        std::size_t expected_count = 1;
        std::size_t actual_count = StringUtilities::CountValuesUntilPrefix(input, prefix);
        KRATOS_EXPECT_EQ(actual_count, expected_count);
    }

    // Test case 2: Prefix not found
    {
        std::string input = "value1 value2 value3";
        std::string prefix = "value4";
        std::size_t expected_count = 3;
        std::size_t actual_count = StringUtilities::CountValuesUntilPrefix(input, prefix);
        KRATOS_EXPECT_EQ(actual_count, expected_count);
    }

    // Test case 3: Empty input
    {
        std::string input = "";
        std::string prefix = "value1";
        std::size_t expected_count = 0;
        std::size_t actual_count = StringUtilities::CountValuesUntilPrefix(input, prefix);
        KRATOS_EXPECT_EQ(actual_count, expected_count);
    }

    // Test case 4: Empty prefix
    {
        std::string input = "value1 value2 value3";
        std::string prefix = "";
        std::size_t expected_count = 0;
        std::size_t actual_count = StringUtilities::CountValuesUntilPrefix(input, prefix);
        KRATOS_EXPECT_EQ(actual_count, expected_count);
    }
}

} // namespace Kratos::Testing
