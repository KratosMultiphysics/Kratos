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

} // namespace Kratos::Testing
